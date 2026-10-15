"""Systematic GEO search for a disease's expression series, via NCBI E-utilities.
Writes a candidate inventory for manual triage.

Disease-agnostic: the only thing that changes between projects is the list of
query spellings. Getting that list wrong is the main way a search misses data --
AS cohorts hide under "axial spondyloarthritis" and "spondyloarthropathy", and
IBD cohorts under "Crohn", "ulcerative colitis" and "inflammatory bowel".

    python scripts/geo_search.py                         # ankylosing spondylitis
    python scripts/geo_search.py --disease ibd
    python scripts/geo_search.py --query "sarcoidosis" --query "Lofgren" \
        --out data/sarcoid_candidates.tsv
"""

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "michael@getwaverly.ai"
TOOL = "as-experiment"

# Cast wide: a disease is rarely filed under one name. Each entry is one GEO
# query; results are pooled and de-duplicated by accession.
DISEASES = {
    "as": [
        "ankylosing spondylitis",
        "axial spondyloarthritis",
        "spondyloarthritis",
        "spondyloarthropathy",
        "radiographic axial spondyloarthritis",
    ],
    "ibd": [
        "inflammatory bowel disease",
        "Crohn disease",
        "Crohn's disease",
        "ulcerative colitis",
        "ileal Crohn",
        "colitis",
        "IBD biopsy",
    ],
}
DEFAULT_DISEASE = "as"

# Restrict to series (GSE); the entry-type filter avoids GDS/GPL/GSM clutter.
# Organism is filtered afterwards on the `taxon` field rather than in the query,
# because the [Organism] index drops series whose taxon string is a mix.
FILTER = ' AND "gse"[Entry Type]'


def _get(endpoint: str, **params) -> str:
    params.setdefault("tool", TOOL)
    params.setdefault("email", EMAIL)
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as fh:
        return fh.read().decode("utf-8", "replace")


def esearch(term: str, retmax: int = 500) -> list[str]:
    raw = _get("esearch.fcgi", db="gds", term=term, retmax=retmax, retmode="json")
    return json.loads(raw)["esearchresult"].get("idlist", [])


def esummary(uids: list[str]) -> list[dict]:
    out = []
    for i in range(0, len(uids), 100):
        chunk = uids[i : i + 100]
        raw = _get("esummary.fcgi", db="gds", id=",".join(chunk), retmode="json")
        res = json.loads(raw)["result"]
        out.extend(res[u] for u in res.get("uids", []))
        time.sleep(0.4)  # NCBI: <=3 req/s without an API key
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--disease", choices=sorted(DISEASES), default=DEFAULT_DISEASE,
                   help="preset list of query spellings")
    p.add_argument("--query", action="append", default=[],
                   help="extra search term; repeatable. Overrides --disease.")
    p.add_argument("--out", type=Path, default=None,
                   help="output TSV (default data/<disease>_candidates.tsv)")
    p.add_argument("--min-samples", type=int, default=10)
    return p.parse_args()


def main() -> None:
    # Series titles contain non-latin-1 characters; the Windows console default
    # codec chokes on them.
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = parse_args()

    terms = args.query or DISEASES[args.disease]
    queries = [f'"{t}"[All Fields]' for t in terms]
    out_path = args.out or (
        Path(__file__).resolve().parents[1] / "data" /
        ("geo_candidates.tsv" if args.disease == "as"
         else f"{args.disease}_geo_candidates.tsv")
    )
    records: dict[str, dict] = {}

    for q in queries:
        term = q + FILTER
        uids = esearch(term)
        print(f"{len(uids):4d} hits  {q}")
        for rec in esummary(uids):
            acc = rec.get("accession", "")
            if not acc.startswith("GSE"):
                continue
            prev = records.get(acc)
            if prev is not None:
                prev["matched_queries"].append(q)
                continue
            records[acc] = {
                "accession": acc,
                "n_samples": int(rec.get("n_samples") or 0),
                "title": (rec.get("title") or "").strip(),
                "gdstype": rec.get("gdstype", ""),
                "taxon": rec.get("taxon", ""),
                "pdat": rec.get("pdat", ""),
                "summary": (rec.get("summary") or "").replace("\n", " ").strip(),
                "platform": ";".join(
                    f"GPL{p}" for p in str(rec.get("gpl", "")).split(";") if p
                ),
                "matched_queries": [q],
            }

    df = pd.DataFrame(records.values())
    df["matched_queries"] = df["matched_queries"].map(lambda x: " | ".join(x))

    # Keep expression / methylation profiling; drop pure genotyping, SNP, ChIP.
    keep = df["gdstype"].str.contains(
        "Expression profiling|Methylation profiling|Non-coding RNA profiling",
        case=False,
        na=False,
    )
    keep &= df["taxon"].str.contains("Homo sapiens", na=False)
    df = df[keep].copy()

    # Flag the series whose text suggests a case/control design -- the thing we
    # actually need. This is a triage hint only; labels get verified per study.
    text = (df["title"] + " " + df["summary"]).str.lower()
    df["hints_control"] = text.str.contains(
        r"healthy|control|normal donor|hc\b|versus|compared", regex=True, na=False
    )
    # Build the disease hint from the query terms themselves, so this works for
    # whatever disease was searched.
    stems = sorted({t.lower().split()[0].rstrip("'s") for t in terms})
    df["hints_disease"] = text.str.contains("|".join(stems), regex=True, na=False)
    df["hints_bulk"] = ~text.str.contains(
        r"single[- ]cell|scrna|single[- ]nucle|citeseq|cite-seq", regex=True, na=False
    )

    df = df.sort_values(
        ["hints_disease", "hints_control", "n_samples"], ascending=False
    )

    out = out_path
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"\n{len(df)} expression/methylation series -> {out}")

    strong = df[
        df["hints_disease"] & df["hints_control"] & (df["n_samples"] >= args.min_samples)
    ]
    print(f"\n{len(strong)} with disease + control hints and "
          f"n>={args.min_samples}:\n")
    with pd.option_context("display.max_colwidth", 95, "display.width", 250):
        print(strong[["accession", "n_samples", "gdstype", "title"]].to_string(index=False))


if __name__ == "__main__":
    main()
