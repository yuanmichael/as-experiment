"""Systematic GEO search for ankylosing spondylitis / axial spondyloarthritis
expression series, via NCBI E-utilities. Writes a candidate inventory to
data/geo_candidates.tsv for manual triage.

Usage:  python scripts/geo_search.py
"""

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

# Cast wide: AS has several names in the literature, and the newer cohorts are
# usually labelled "axial spondyloarthritis" rather than "ankylosing spondylitis".
QUERIES = [
    '"ankylosing spondylitis"[All Fields]',
    '"axial spondyloarthritis"[All Fields]',
    '"spondyloarthritis"[All Fields]',
    '"spondyloarthropathy"[All Fields]',
    '"radiographic axial spondyloarthritis"[All Fields]',
]

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


def main() -> None:
    # Series titles contain non-latin-1 characters; the Windows console default
    # codec chokes on them.
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    records: dict[str, dict] = {}

    for q in QUERIES:
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
    df["hints_as"] = text.str.contains(
        r"ankylosing|spondyloarthrit|spondylitis|\bas\b|axspa", regex=True, na=False
    )
    df["hints_bulk"] = ~text.str.contains(
        r"single[- ]cell|scrna|single[- ]nucle|citeseq|cite-seq", regex=True, na=False
    )

    df = df.sort_values(["hints_as", "hints_control", "n_samples"], ascending=False)

    out = Path(__file__).resolve().parents[1] / "data" / "geo_candidates.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\n{len(df)} expression/methylation series -> {out}")

    strong = df[df["hints_as"] & df["hints_control"] & (df["n_samples"] >= 10)]
    print(f"\n{len(strong)} with AS + control hints and n>=10:\n")
    with pd.option_context("display.max_colwidth", 95, "display.width", 250):
        print(strong[["accession", "n_samples", "gdstype", "title"]].to_string(index=False))


if __name__ == "__main__":
    main()
