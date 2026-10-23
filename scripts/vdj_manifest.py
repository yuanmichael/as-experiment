"""Build a download manifest for paired GEX + V(D)J single-cell data.

Everything downstream depends on one operation: joining a cell's receptor to
that same cell's transcriptome by its 10x barcode. This script works out, for
each sample in each study, which files are needed and where they actually live
-- without downloading any of them.

Why this is the hard part. GEO has no schema for this. Across five studies we
have already hit:

  * files under the SAMPLE directory (GSM-prefixed, bundled in a series tar)
    versus the SERIES directory (GSE-prefixed, deposited loose)
  * series with filelist.txt and series without
  * filelist.txt columns ordered Archive/File | Name | Time | Size | Type,
    so the name is the SECOND field
  * sorted libraries: Bac17L_CD3pos and Bac17L_CD3neg share a sample name, and
    pairing the TCR file to the CD3-NEGATIVE barcodes gives ~1% overlap, which
    reads as poor recovery rather than as a mismatch
  * expression as a barcodes/features/matrix triple OR a single .h5

A sample is only usable if we can resolve a TCR contig file AND an expression
object for the SAME library. Anything else is reported as a gap rather than
silently dropped -- a silently dropped sample looks like a small cohort.

Usage:
    python -I scripts/vdj_manifest.py
    python -I scripts/vdj_manifest.py --series GSE167118 GSE266519
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "airr"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
UA = {"User-Agent": "as-experiment/1.0 (manifest build; metadata only)"}

# v1 training set. Chosen by what the manifest could actually RESOLVE, not by
# title or sample count -- of 19 candidates checked, 5 resolved nothing at all:
#   GSE161918  112 .h5 files, zero TCR files. Not a paired dataset; the earlier
#              text-based sweep false-positived on "CSP"/"feature_bc" names.
#   GSE270518  163 contig files, zero expression files. Its GEX sits in a
#              sibling series under a SuperSeries.
#   GSE270387, GSE296419, GSE350364, GSE311112, GSE234261  same shape.
# Deliberately spread across cancer, mucosal inflammation, infection challenge,
# alloimmunity and a developmental baseline, so the model does not learn one
# disease's expression signature as if it were T cell biology.
TRAIN = [
    "GSE272993",   # 121 libs, 1,321,532 cells, checkpoint blockade
    "GSE329527",   #  71 libs,   200,597 cells, Candida Th17, oral-gut
    "GSE271413",   #  27 libs,   106,819 cells, immune atlas, premature infants
    "GSE167118",   #  22 libs,   104,062 cells, COVID lung Th17
    "GSE301133",   #  34 libs,    28,055 cells, aerosolised BCG challenge
    "GSE266519",   #   1 lib,     73,674 cells, HSCT + healthy donors
    "GSE292177",   #  15 libs,    28,589 cells  -- DROPPED, see below
    "GSE256141",   #  30 libs,     3,029 cells  -- DROPPED, see below
]

# Dropped after measurement, not before. Both are real paired data; they are
# simply mostly not T cells, so the cost per usable cell is prohibitive and
# their pairing rates swamp the distribution used to spot genuine mismatches.
#   GSE256141  liver transplant biopsies: 3,029 cells from 30 libraries,
#              1.1% median pairing
#   GSE292177  CD45-sorted lymphoma: 28,589 cells from 15 libraries,
#              4.6% median pairing
# Leaves 276 training libraries and 1,862,114 cells.
DROP = ["GSE256141", "GSE292177"]

# Never enters training. GSE288581 is the AS out-of-distribution test.
TEST = ["GSE288581"]

BARCODES = re.compile(r"barcodes\.tsv(\.gz)?$", re.I)
FEATURES = re.compile(r"(features|genes)\.tsv(\.gz)?$", re.I)
MATRIX = re.compile(r"matrix\.mtx(\.gz)?$", re.I)
H5 = re.compile(r"(filtered|raw)_feature_bc_matrix\.h5$|\.h5ad(\.gz)?$", re.I)
CONTIG = re.compile(r"(filtered|all)_contig_annotations[^/]*\.csv(\.gz)?$", re.I)

# Tags that differ between the GEX and TCR file of the same library.
# CD3pos/CD3neg are NOT in here: they identify which sorted library the
# barcodes came from, and stripping them pairs TCR against non-T cells.
MODIFIER = re.compile(
    r"(?:^|[_-])(gex|rna|5gex|scrna|tcr|vdj[_-]?t|vdj|t[_-]?cell|"
    r"filtered|raw|all)(?=[_-]|$)", re.I)

# A TCR library belongs with the sample's T-cell GEX library: untagged, or
# CD3-positive. Tried in order.
T_SORTED = ["", "_cd3pos", "_cd3", "_cd3p", "_tcells", "_tcell", "_cd3_pos"]

# BCR files would pair against the same GEX barcodes and are not wanted here.
IS_BCR = re.compile(r"vdj[_-]?b|bcr|_ig[hkl][_-]|bcell", re.I)


def fetch_text(url: str, timeout: int = 90) -> str | None:
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=timeout) as fh:
            return fh.read().decode("utf-8", "replace")
    except Exception:
        return None


def discover(gse: str) -> list[tuple[str, int]]:
    """(filename, bytes) for a series, from filelist.txt or the directory."""
    stub = gse[:-3] + "nnn"
    out: list[tuple[str, int]] = []

    text = fetch_text(f"{FTP}/series/{stub}/{gse}/suppl/filelist.txt")
    if text:
        for line in text.splitlines()[1:]:
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 2 or not parts[1] or parts[1].endswith(".tar"):
                continue
            size = 0
            for p in parts[2:]:
                q = p.replace(",", "")
                if q.isdigit():
                    size = int(q)
                    break
            out.append((parts[1], size))

    if not out:
        html = fetch_text(f"{FTP}/series/{stub}/{gse}/suppl/")
        if html:
            for m in re.finditer(
                    r'<a href="([^"?/][^"]*)"[^>]*>[^<]*</a>\s*'
                    r'(?:\d{4}-\d\d-\d\d \d\d:\d\d)?\s*([\d.]+[KMGT]?)?', html):
                name, size = m.group(1), m.group(2)
                if name.startswith("http") or name.endswith(".tar") \
                        or name.lower() == "filelist.txt":
                    continue
                out.append((name, parse_size(size)))
    return out


def parse_size(s: str | None) -> int:
    if not s:
        return 0
    mult = {"K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
    try:
        return int(float(s[:-1]) * mult[s[-1]]) if s[-1] in mult else int(float(s))
    except Exception:
        return 0


def url_for(name: str, gse: str) -> str | None:
    """Where the file is actually servable from."""
    m = re.match(r"(GSM\d+)", name)
    if m:
        gsm = m.group(1)
        return f"{FTP}/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"
    if name.startswith(gse):
        return f"{FTP}/series/{gse[:-3]}nnn/{gse}/suppl/{name}"
    return None


def token(name: str, suffix: re.Pattern) -> str:
    """Library identity: strip the accession prefix and the file-type suffix,
    then strip modality tags until nothing changes (they stack)."""
    s = re.sub(r"^(GSM\d+|GSE\d+)[_-]", "", name)
    s = suffix.sub("", s)
    prev = None
    while prev != s:
        prev = s
        s = MODIFIER.sub("", s)
    return re.sub(r"[_-]+", "_", s).strip("_-.").lower()


BC_SFX = re.compile(r"[_-]?barcodes\.tsv(\.gz)?$", re.I)
FT_SFX = re.compile(r"[_-]?(features|genes)\.tsv(\.gz)?$", re.I)
MX_SFX = re.compile(r"[_-]?matrix\.mtx(\.gz)?$", re.I)
H5_SFX = re.compile(r"[_-]?((filtered|raw)_feature_bc_matrix\.h5|\.h5ad(\.gz)?)$",
                    re.I)
CT_SFX = re.compile(r"[_-]?(filtered|all)_contig_annotations[^.]*\.csv(\.gz)?$",
                    re.I)


def build(gse: str) -> tuple[list[dict], dict]:
    files = discover(gse)
    if not files:
        return [], dict(gse=gse, note="no file listing")

    size_of = {n: s for n, s in files}
    buckets: dict[str, dict[str, str]] = defaultdict(dict)
    for name, _ in files:
        if BARCODES.search(name):
            buckets[token(name, BC_SFX)].setdefault("barcodes", name)
        elif FEATURES.search(name):
            buckets[token(name, FT_SFX)].setdefault("features", name)
        elif MATRIX.search(name):
            buckets[token(name, MX_SFX)].setdefault("matrix", name)
        elif H5.search(name):
            buckets[token(name, H5_SFX)].setdefault("h5", name)
        elif CONTIG.search(name) and not IS_BCR.search(name):
            buckets[token(name, CT_SFX)].setdefault("contig", name)

    rows, unmatched_contig = [], []
    for tok, b in sorted(buckets.items()):
        if "contig" not in b:
            continue
        # Find the expression object for this library: same token, or the
        # T-cell-sorted variant of it.
        gex_tok = None
        for sfx in T_SORTED:
            cand = buckets.get(tok + sfx, {})
            if "h5" in cand or ("barcodes" in cand and "matrix" in cand):
                gex_tok = tok + sfx
                break
        if gex_tok is None:
            unmatched_contig.append(b["contig"])
            continue

        g = buckets[gex_tok]
        need = [b["contig"]] + ([g["h5"]] if "h5" in g else
                                [g.get("barcodes"), g.get("features"),
                                 g.get("matrix")])
        need = [n for n in need if n]
        urls = {n: url_for(n, gse) for n in need}
        rows.append(dict(
            gse=gse, library=tok or "(series)", gex_library=gex_tok or "(series)",
            contig=b["contig"],
            barcodes=g.get("barcodes", ""), features=g.get("features", ""),
            matrix=g.get("matrix", ""), h5=g.get("h5", ""),
            bytes=sum(size_of.get(n, 0) for n in need),
            unresolved_urls=sum(1 for u in urls.values() if not u),
        ))

    stats = dict(
        gse=gse, n_files=len(files), n_libraries=len(buckets),
        n_usable=len(rows), n_contig_without_gex=len(unmatched_contig),
        gb=round(sum(r["bytes"] for r in rows) / 1e9, 2),
        note="",
    )
    return rows, stats


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", nargs="*", default=None)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    series = args.series or [s for s in TRAIN + TEST if s not in DROP]
    all_rows, all_stats = [], []
    for gse in series:
        rows, stats = build(gse)
        split = "test" if gse in TEST else "train"
        for r in rows:
            r["split"] = split
        stats["split"] = split
        all_rows += rows
        all_stats.append(stats)
        print(f"{gse:12s} [{split:5s}] {stats.get('n_files', 0):4d} files, "
              f"{stats.get('n_libraries', 0):4d} libraries -> "
              f"{stats['n_usable']:4d} usable  {stats['gb']:7.2f} GB"
              + (f"   ({stats['n_contig_without_gex']} TCR files with no GEX match)"
                 if stats.get("n_contig_without_gex") else "")
              + (f"   {stats['note']}" if stats.get("note") else ""))

    man = pd.DataFrame(all_rows)
    st = pd.DataFrame(all_stats)
    print("\n" + "=" * 92)
    print("MANIFEST")
    print("=" * 92)
    if man.empty:
        print("  nothing resolved")
        return
    # 6,400 cells/library is MEASURED across the 291 libraries actually built
    # (1,862,114 cells), not the 12k/library guessed from 7 T-sorted samples.
    # The guess was 2x optimistic because tissue biopsies are mostly not T cells.
    for split, g in man.groupby("split"):
        print(f"  {split:5s}  {len(g):4d} libraries  "
              f"{g['bytes'].sum()/1e9:7.2f} GB  "
              f"~{6_400*len(g):,} cells at 6.4k/library (measured)")
    bad = int(man["unresolved_urls"].sum())
    if bad:
        print(f"\n  WARNING: {bad} files have no resolvable URL")

    man.to_csv(OUT / "vdj_manifest.tsv", sep="\t", index=False)
    st.to_csv(OUT / "vdj_manifest_stats.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT / 'vdj_manifest.tsv'}")
    print("\nCell counts extrapolate the measured 6.4k/library. Per-library truth")
    print("is in the build's stats.json on the Modal volume; a library far below")
    print("the median pairing rate needs its T-cell enrichment checked before use.")


if __name__ == "__main__":
    main()
