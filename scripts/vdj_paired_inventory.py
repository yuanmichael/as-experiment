"""Which public datasets can train a cell-level TCR-sequence -> cell-state model?

The requirement is narrow and non-negotiable: both the TCR sequence AND the
transcriptome for the SAME cell, linked by the 10x cell barcode. 3' chemistry
cannot do this (no V segment -- see LITERATURE_REVIEW.md Section 9). Only 5'
chemistry with V(D)J enrichment can.

Detection: GEO supplementary file names. 10x Cell Ranger emits
  filtered_contig_annotations.csv / clonotypes.csv   <- VDJ, carries the barcode
  matrix.mtx / features.tsv / barcodes.tsv           <- GEX, carries the barcode
Both present in one series => the barcode join is possible.

The point of looking at PROCESSED files rather than raw reads: a
filtered_contig_annotations.csv is a few MB. The FASTQs behind it are tens of
GB. If the processed outputs are deposited, the training data is cheap.

Downloads nothing. Reads directory listings and file sizes only.

Usage:  python scripts/vdj_paired_inventory.py [--limit 200]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "airr"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
UA = {"User-Agent": "as-experiment/1.0 (paired VDJ inventory; metadata only)"}

VDJ_FILE = re.compile(
    r"contig_annot|clonotype|_vdj|vdj_|airr|tcr.*\.csv|filtered_contig", re.I)
GEX_FILE = re.compile(
    r"matrix\.mtx|features\.tsv|genes\.tsv|barcodes\.tsv|_counts|expression"
    r"|\.h5ad|filtered_feature_bc", re.I)

# Series summary text that indicates 5' + VDJ rather than 3'
PAIRED_TEXT = re.compile(
    r"5'|5 prime|five.prime|v\(d\)j|vdj|immune profiling|tcr enrich", re.I)
THREE_PRIME = re.compile(r"\b3'\b|3 prime|three.prime", re.I)


def get(url: str, timeout: int = 60) -> bytes | None:
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                    timeout=timeout) as fh:
            return fh.read()
    except Exception:
        return None


def esearch(term: str, retmax: int) -> list[str]:
    raw = get(f"{EUTILS}/esearch.fcgi?db=gds&retmax={retmax}&retmode=json"
              f"&term={urllib.parse.quote(term)}")
    if not raw:
        return []
    try:
        return json.loads(raw)["esearchresult"].get("idlist", []) or []
    except Exception:
        return []


def esummary(uids: list[str]) -> list[dict]:
    rows = []
    for i in range(0, len(uids), 150):
        raw = get(f"{EUTILS}/esummary.fcgi?db=gds&retmode=json"
                  f"&id={','.join(uids[i:i+150])}")
        time.sleep(0.4)
        if not raw:
            continue
        try:
            res = json.loads(raw).get("result", {})
        except Exception:
            continue
        for uid in res.get("uids", []):
            r = res.get(uid, {})
            if not str(r.get("accession", "")).startswith("GSE"):
                continue
            rows.append(dict(acc=r["accession"],
                             title=(r.get("title") or "").strip(),
                             summary=(r.get("summary") or "").strip(),
                             n_samples=int(r.get("n_samples") or 0),
                             taxon=r.get("taxon", ""),
                             pdat=r.get("PDAT", "")))
    return rows


def suppl_listing(gse: str) -> list[tuple[str, int]]:
    """Filenames and byte sizes for a series' supplementary data.

    Most series deposit one `GSE*_RAW.tar` plus `filelist.txt`; the per-sample
    10x outputs are named only inside that filelist. So read the directory, and
    whenever a filelist is there, read it too -- otherwise every 10x series
    looks like a single opaque tarball and the VDJ/GEX detector sees nothing.
    """
    stub = gse[:-3] + "nnn" if len(gse) > 6 else "GSEnnn"
    base = f"{FTP}/{stub}/{gse}/suppl"
    raw = get(f"{base}/", timeout=45)
    if not raw:
        return []
    html = raw.decode("utf-8", "replace")
    out, has_filelist = [], False
    for m in re.finditer(r'<a href="([^"?/][^"]*)"[^>]*>[^<]*</a>\s*'
                         r'(\d{4}-\d\d-\d\d \d\d:\d\d)?\s*([\d.]+[KMGT]?)?',
                         html):
        name, _, size = m.groups()
        if name.startswith("http"):
            continue
        if name.lower() == "filelist.txt":
            has_filelist = True
            continue
        out.append((name, parse_size(size)))

    if has_filelist:
        fl = get(f"{base}/filelist.txt", timeout=45)
        if fl:
            text = fl.decode("utf-8", "replace")
            # Columns are: Archive/File | Name | Time | Size | Type
            # -- the NAME is the second field, not the first.
            for line in text.splitlines()[1:]:
                parts = [p.strip() for p in line.split("\t")]
                if len(parts) < 2 or not parts[1]:
                    continue
                sz = 0
                for p in parts[2:]:
                    q = p.replace(",", "")
                    if q.isdigit():
                        sz = int(q)
                        break
                out.append((parts[1], sz))
    return out


def parse_size(s: str | None) -> int:
    if not s:
        return 0
    mult = {"K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
    try:
        if s[-1] in mult:
            return int(float(s[:-1]) * mult[s[-1]])
        return int(float(s))
    except Exception:
        return 0


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=200,
                    help="how many candidate series to check file listings for")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    # Text queries have poor recall here: GSE288581 deposits
    # filtered_contig_annotations AND matrix.mtx but never says "VDJ" in its
    # summary, so a keyword sweep misses it. The filelist IS the ground truth
    # and costs ~2 requests per series, so cast the net wide and let the file
    # names decide.
    queries = {
        "10x 5' VDJ, T cells": (
            '("V(D)J"[All Fields] OR "VDJ"[All Fields]) AND '
            '("single cell"[All Fields]) AND ("T cell"[All Fields] OR "TCR"[All Fields]) '
            'AND "Homo sapiens"[Organism]'),
        "immune profiling 5'": (
            '("immune profiling"[All Fields] OR "5\' gene expression"[All Fields]) '
            'AND "Homo sapiens"[Organism]'),
        "paired GEX + TCR": (
            '("paired"[All Fields] AND "TCR"[All Fields] AND "transcriptome"[All Fields]) '
            'AND "Homo sapiens"[Organism]'),
        "sc T cell / PBMC (broad)": (
            '("single cell RNA"[All Fields] OR "scRNA-seq"[All Fields]) AND '
            '("T cell"[All Fields] OR "PBMC"[All Fields] OR "lymphocyte"[All Fields] '
            'OR "synovial"[All Fields]) AND "Homo sapiens"[Organism] '
            'AND "gse"[Entry Type]'),
        "clonal expansion / clonotype": (
            '("clonally expanded"[All Fields] OR "clonotype"[All Fields] OR '
            '"clonal expansion"[All Fields]) AND "Homo sapiens"[Organism] '
            'AND "gse"[Entry Type]'),
    }

    # Accessions we already care about, in case the queries miss them.
    SEEDS = ["GSE288581", "GSE272431", "GSE178991", "GSE207037", "GSE290921",
             "GSE216914", "GSE214974", "GSE232802", "GSE274471", "GSE134649",
             "GSE314223", "GSE342939", "GSE132832", "GSE308968", "GSE197456",
             "GSE291328", "GSE157277", "GSE320359", "GSE283500", "GSE256256"]

    seen, rows = set(), []
    for acc in SEEDS:
        uids = esearch(f"{acc}[ACCN]", 5)
        time.sleep(0.3)
        for r in esummary([u for u in uids if u not in seen]):
            seen.update(uids)
            r["found_by"] = "seed"
            rows.append(r)
    print(f"  {'seeds':26s} {len(rows):4d} series")
    print("=" * 92)
    print("GEO CANDIDATES")
    print("=" * 92)
    for label, term in queries.items():
        uids = [u for u in esearch(term, 400) if u not in seen]
        seen.update(uids)
        got = esummary(uids)
        for r in got:
            r["found_by"] = label
        rows += got
        print(f"  {label:26s} {len(uids):4d} new uids -> {len(got):4d} GSE")

    df = pd.DataFrame(rows).drop_duplicates(subset=["acc"])
    df["text"] = df["title"] + " || " + df["summary"]
    df["looks_5prime"] = df["text"].str.contains(PAIRED_TEXT)
    df["is_seed"] = df["found_by"] == "seed"
    df = df.sort_values(["is_seed", "looks_5prime", "n_samples"],
                        ascending=[False, False, False])
    print(f"\n{len(df)} distinct series, {int(df['looks_5prime'].sum())} "
          f"mentioning 5'/VDJ in the text")

    # The expensive-ish part: one directory listing per series
    cand = df.head(args.limit)
    print("\n" + "=" * 92)
    print(f"CHECKING SUPPLEMENTARY FILES for {len(cand)} series")
    print("=" * 92)
    recs = []
    for i, r in enumerate(cand.itertuples(), 1):
        files = suppl_listing(r.acc)
        time.sleep(0.25)
        names = " ".join(n for n, _ in files)
        has_vdj = bool(VDJ_FILE.search(names))
        has_gex = bool(GEX_FILE.search(names))
        total = sum(sz for _, sz in files)
        recs.append(dict(acc=r.acc, title=r.title, n_samples=r.n_samples,
                         n_files=len(files), suppl_gb=round(total / 1e9, 3),
                         has_vdj=has_vdj, has_gex=has_gex,
                         paired=has_vdj and has_gex,
                         looks_5prime=r.looks_5prime, pdat=r.pdat))
        if has_vdj and has_gex:
            print(f"  [{i:3d}] PAIRED  {r.acc:12s} n={r.n_samples:4d} "
                  f"{total/1e9:7.2f} GB  {r.title[:70]}")

    res = pd.DataFrame(recs)
    res.to_csv(OUT / "vdj_paired_inventory.tsv", sep="\t", index=False)

    paired = res[res["paired"]].sort_values("n_samples", ascending=False)
    print("\n" + "=" * 92)
    print("SUMMARY")
    print("=" * 92)
    print(f"  checked                       {len(res)}")
    print(f"  VDJ files present             {int(res['has_vdj'].sum())}")
    print(f"  GEX files present             {int(res['has_gex'].sum())}")
    print(f"  BOTH (barcode join possible)  {len(paired)}")
    print(f"  total samples in paired set   {int(paired['n_samples'].sum()):,}")
    print(f"  total processed bytes         {paired['suppl_gb'].sum():.1f} GB")

    if len(paired):
        print("\nTop 30 paired series by sample count:")
        with pd.option_context("display.width", 200, "display.max_colwidth", 68):
            print(paired[["acc", "n_samples", "suppl_gb", "title"]]
                  .head(30).to_string(index=False))

    print(f"\nwrote {OUT / 'vdj_paired_inventory.tsv'}")
    print("\nNOTE: file-name detection is a proxy. 'paired' means a VDJ-shaped and")
    print("a GEX-shaped file both exist in the series, NOT that their barcodes")
    print("actually overlap. That has to be checked by opening one of each.")


if __name__ == "__main__":
    main()
