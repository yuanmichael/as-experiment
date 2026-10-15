"""Step 1 of the bulk meta-analysis: verify sample labels directly.

For each candidate series, pull every GSM's title + characteristics from GEO
(text view of the accession page) and classify each sample as AS / other-disease
/ healthy-control by regex over its own metadata. Prints a per-study inventory
and writes one long table of every sample for manual inspection.

Usage:  python scripts/geo_inventory.py [GSE1 GSE2 ...]
"""

from __future__ import annotations

import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

# Everything from the GEO search that is human, expression/methylation, and
# could plausibly contain both AS cases and controls.
CANDIDATES = [
    # blood / immune cells -- the tissues we actually want
    "GSE73754", "GSE25101", "GSE18781", "GSE221786", "GSE141646",
    "GSE299639", "GSE243689", "GSE266295", "GSE232131", "GSE122624",
    "GSE11886", "GSE252708", "GSE223717", "GSE181364", "GSE205812",
    "GSE178408", "GSE220130",
    # tissue / other modality -- inventoried for completeness, not pooled
    "GSE41038", "GSE30023", "GSE58667", "GSE14442", "GSE179571",
    "GSE304063",
    # single cell, for project 2
    "GSE194315", "GSE288581", "GSE163314",
]

HEALTHY = re.compile(
    r"\b(healthy|normal|control|hc\b|nc\b|ctrl|unaffected|non[- ]?diseased)\b", re.I
)
AS_CASE = re.compile(
    r"(ankylosing\s+spondylitis|\bas\s+patient|\bas[-_ ]?\d|axial\s+spondyloarthrit"
    r"|\bax[- ]?spa\b|\bspondyloarthrit|\bspa\b|\bas\b)",
    re.I,
)
OTHER_DISEASE = re.compile(
    r"(rheumatoid|psoriatic|psoriasis|crohn|ulcerative|sarcoid|lupus|\bsle\b"
    r"|osteoarthrit|\bra\b|\bpsa\b|septic|gout|undifferentiated)",
    re.I,
)
TISSUE = re.compile(
    r"(whole\s+blood|peripheral\s+blood|\bpbmc\b|\bpbl\b|serum|plasma|monocyte"
    r"|macrophage|synovial\s+fluid|synovi\w*|gut|colon|ileum|bone\s+marrow"
    r"|entheseal|mesenchymal|neutrophil|\bt[- ]cell|cd8|cd4|exosome)",
    re.I,
)


def fetch_samples(gse: str) -> str:
    """All GSM records of a series, as GEO's brief text view."""
    url = f"{ACC_CGI}?{urllib.parse.urlencode({'acc': gse, 'targ': 'gsm', 'form': 'text', 'view': 'brief'})}"
    req = urllib.request.Request(url, headers={"User-Agent": "as-experiment/1.0"})
    with urllib.request.urlopen(req, timeout=120) as fh:
        return fh.read().decode("utf-8", "replace")


def parse_samples(text: str) -> list[dict]:
    """Split the text view into one dict per ^SAMPLE block."""
    samples: list[dict] = []
    cur: dict | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                samples.append(cur)
            cur = {"gsm": line.split("=", 1)[1].strip(), "fields": []}
        elif cur is not None and line.startswith("!Sample_"):
            key, _, val = line.partition("=")
            key = key.strip().removeprefix("!Sample_")
            val = val.strip()
            if key in {
                "title", "source_name_ch1", "characteristics_ch1",
                "description", "library_strategy", "platform_id",
            }:
                cur["fields"].append(f"{key}: {val}")
    if cur:
        samples.append(cur)
    return samples


def classify(blob: str) -> str:
    """Label one sample from its own metadata. Healthy wins over AS, because
    'healthy control' strings in an AS study also match the AS pattern."""
    if HEALTHY.search(blob):
        return "control"
    other = OTHER_DISEASE.search(blob)
    as_hit = AS_CASE.search(blob)
    if as_hit and not other:
        return "AS"
    if other and not as_hit:
        return f"other:{other.group(1).lower()}"
    if other and as_hit:
        return f"ambiguous:{other.group(1).lower()}"
    return "unknown"


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    accessions = sys.argv[1:] or CANDIDATES

    rows: list[dict] = []
    for gse in accessions:
        try:
            samples = parse_samples(fetch_samples(gse))
        except Exception as exc:  # network / malformed page
            print(f"{gse}: FETCH FAILED ({exc})")
            continue
        for s in samples:
            blob = " | ".join(s["fields"])
            tissues = sorted({m.group(0).lower() for m in TISSUE.finditer(blob)})
            rows.append(
                {
                    "gse": gse,
                    "gsm": s["gsm"],
                    "label": classify(blob),
                    "tissue": ";".join(tissues),
                    "meta": blob,
                }
            )
        print(f"{gse}: {len(samples)} samples fetched")
        time.sleep(0.5)

    df = pd.DataFrame(rows)
    out = Path(__file__).resolve().parents[1] / "data" / "geo_sample_inventory.tsv"
    df.to_csv(out, sep="\t", index=False)

    print(f"\nwrote {len(df)} samples -> {out}\n")
    print("=" * 100)
    for gse, g in df.groupby("gse", sort=False):
        counts = g["label"].value_counts().to_dict()
        n_as = counts.get("AS", 0)
        n_ctl = counts.get("control", 0)
        verdict = "USABLE" if (n_as >= 4 and n_ctl >= 4) else "drop/inspect"
        tis = ";".join(sorted({t for s in g["tissue"] for t in s.split(";") if t}))
        print(f"\n{gse}  n={len(g)}  AS={n_as}  control={n_ctl}  -> {verdict}")
        print(f"  labels: {counts}")
        print(f"  tissue: {tis[:120]}")
        print(f"  e.g.  : {g['meta'].iloc[0][:200]}")


if __name__ == "__main__":
    main()
