"""Download the per-study expression matrices for the AS-vs-control meta-analysis.

Microarray studies: GEO series matrix (submitter-normalised, log or linear).
RNA-seq studies:    supplementary count tables.

Nothing is merged here. Each study lands in its own directory, with a sample
sheet carrying the labels verified in geo_inventory.py.

Usage:  python scripts/bulk_download.py [GSE1 GSE2 ...]
"""

from __future__ import annotations

import gzip
import io
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "bulk"

# series matrix lives under a padded range directory, e.g. GSE73nnn/GSE73754
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{gse}/matrix/{gse}_series_matrix.txt.gz"
)
SUPPL_DIR = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{gse}/suppl/"

STUDIES = {
    # accession: (kind, tissue, platform note)
    "GSE73754": ("matrix", "whole_blood", "Illumina HT-12 v4 (GPL10558)"),
    "GSE25101": ("matrix", "whole_blood", "Illumina HT-12 v3 (GPL6947)"),
    "GSE18781": ("matrix", "peripheral_blood", "Affymetrix U133 Plus 2.0 (GPL570)"),
    "GSE11886": ("matrix", "mdm", "Affymetrix U133 Plus 2.0 (GPL570)"),
    "GSE299639": ("suppl", "peripheral_blood", "ONT long-read RNA-seq"),
    "GSE181364": ("suppl", "whole_blood", "Illumina RNA-seq"),
    "GSE221786": ("suppl", "il17_enriched", "Illumina RNA-seq, CytoSim-stimulated"),
    "GSE243689": ("suppl", "cd8_sorted", "Illumina RNA-seq"),
}


def stub(gse: str) -> str:
    """GSE73754 -> GSE73 (the FTP range directory prefix)."""
    num = gse.removeprefix("GSE")
    return "GSE" + (num[:-3] if len(num) > 3 else "")


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "as-experiment/1.0"})
    with urllib.request.urlopen(req, timeout=600) as fh:
        return fh.read()


def get_matrix(gse: str, outdir: Path) -> None:
    url = MATRIX_URL.format(stub=stub(gse), gse=gse)
    dest = outdir / f"{gse}_series_matrix.txt"
    if dest.exists():
        print(f"  {dest.name} already present")
        return
    blob = fetch(url)
    text = gzip.decompress(blob).decode("utf-8", "replace")
    dest.write_text(text, encoding="utf-8")
    n = sum(1 for _ in text.splitlines())
    print(f"  {dest.name}  ({len(blob) / 1e6:.1f} MB gz, {n} lines)")


def list_suppl(gse: str) -> list[str]:
    """Scrape the supplementary file listing off the FTP HTML index."""
    import re

    html = fetch(SUPPL_DIR.format(stub=stub(gse), gse=gse)).decode("utf-8", "replace")
    return sorted(set(re.findall(r'href="([^"?/][^"]*)"', html)))


def get_suppl(gse: str, outdir: Path) -> None:
    files = list_suppl(gse)
    print(f"  supplementary files: {files}")
    # Keep the tabular ones; skip BAMs, fastqs, per-sample archives.
    wanted = [
        f
        for f in files
        if f.lower().endswith((".txt.gz", ".csv.gz", ".tsv.gz", ".xlsx", ".txt", ".csv", ".tsv"))
    ]
    for f in wanted:
        dest = outdir / f
        if dest.exists():
            print(f"  {f} already present")
            continue
        blob = fetch(SUPPL_DIR.format(stub=stub(gse), gse=gse) + f)
        dest.write_bytes(blob)
        print(f"  {f}  ({len(blob) / 1e6:.2f} MB)")


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    targets = sys.argv[1:] or list(STUDIES)
    for gse in targets:
        kind, tissue, platform = STUDIES[gse]
        outdir = RAW / gse
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"\n{gse}  [{kind}]  tissue={tissue}  {platform}")
        try:
            (get_matrix if kind == "matrix" else get_suppl)(gse, outdir)
        except Exception as exc:
            print(f"  FAILED: {exc}")


if __name__ == "__main__":
    main()
