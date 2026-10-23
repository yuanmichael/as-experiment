"""Do the GEX and V(D)J barcodes actually join?

The one load-bearing assumption behind a cell-level `TCR sequence -> cell state`
model: that 10x deposits let you recover WHICH receptor belongs to WHICH cell.
"Both file types present in the series" is not that. This opens one of each and
measures the overlap.

What it reports, and why each number matters:

  n_gex          cells in the expression matrix
  n_vdj_cells    distinct barcodes in the VDJ table
  overlap        cells with BOTH -- the actual training set size
  pct_gex_paired overlap / n_gex. If this is 10%, the usable data is 10x
                 smaller than the cell count suggests.
  productive_trb cells with a productive TRB CDR3 -- the real input
  suffix         barcode suffix convention on each side ("-1" etc). A mismatch
                 here silently zeroes the overlap and looks like biology.

Downloads only barcodes.tsv.gz and filtered_contig_annotations.csv.gz -- a few
hundred KB each. Nothing large. Files land in their own directory under
data/vdj_check/ and are parsed, never executed.

Usage:  python -I scripts/vdj_join_check.py
        python -I scripts/vdj_join_check.py --series GSE288581 GSE167118
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import sys
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DL = ROOT / "data" / "vdj_check"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
UA = {"User-Agent": "as-experiment/1.0 (barcode join check)"}

# GSE288581 first: it is the AS test set, so if its barcodes do not join the
# out-of-distribution test is dead and nothing else matters.
DEFAULT = ["GSE288581", "GSE167118", "GSE175450", "GSE222633", "GSE266519"]

MAX_BYTES = 30_000_000   # per file; these should be ~10 KB - 1 MB

BARCODE_FILE = re.compile(r"barcodes\.tsv\.gz$", re.I)
CONTIG_FILE = re.compile(r"(filtered_contig_annotations|all_contig_annotations)"
                         r"[^/]*\.csv\.gz$", re.I)


def fetch(url: str, dest: Path, timeout: int = 120) -> Path | None:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=timeout) as fh:
            size = int(fh.headers.get("Content-Length") or 0)
            if size > MAX_BYTES:
                print(f"    skip (too big: {size/1e6:.1f} MB) {url.split('/')[-1]}")
                return None
            data = fh.read(MAX_BYTES + 1)
        if len(data) > MAX_BYTES:
            print(f"    skip (over cap) {url.split('/')[-1]}")
            return None
        dest.write_bytes(data)
        return dest
    except Exception as exc:
        print(f"    fetch failed ({type(exc).__name__}) {url.split('/')[-1]}")
        return None


def filelist(gse: str) -> list[str]:
    """Series file names, from filelist.txt if there is one and otherwise from
    the directory itself -- series that deposit per-sample files directly have
    no tar and therefore no filelist."""
    stub = gse[:-3] + "nnn"
    names: list[str] = []

    try:
        req = urllib.request.Request(
            f"{FTP}/series/{stub}/{gse}/suppl/filelist.txt", headers=UA)
        with urllib.request.urlopen(req, timeout=90) as fh:
            text = fh.read().decode("utf-8", "replace")
        for line in text.splitlines()[1:]:
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) >= 2 and parts[1] and not parts[1].endswith(".tar"):
                names.append(parts[1])
    except Exception:
        pass

    if not names:
        try:
            req = urllib.request.Request(
                f"{FTP}/series/{stub}/{gse}/suppl/", headers=UA)
            with urllib.request.urlopen(req, timeout=90) as fh:
                html = fh.read().decode("utf-8", "replace")
            for m in re.finditer(r'href="([^"?/][^"]*)"', html):
                n = m.group(1)
                if not n.startswith("http") and not n.endswith(".tar") \
                        and n.lower() != "filelist.txt":
                    names.append(n)
        except Exception:
            pass
    return names


def sample_url(fname: str, gse: str) -> str | None:
    """Where one supplementary file is actually servable.

    GSM-prefixed names sit in the sample directory (files bundled inside a
    series RAW.tar are not individually servable from the series directory, but
    they are from the sample's). GSE-prefixed names sit in the series
    directory."""
    m = re.match(r"(GSM\d+)", fname)
    if m:
        gsm = m.group(1)
        return f"{FTP}/samples/{gsm[:-3]}nnn/{gsm}/suppl/{fname}"
    if fname.startswith(gse):
        return f"{FTP}/series/{gse[:-3]}nnn/{gse}/suppl/{fname}"
    return None


# Modality/sort tags that differ between the GEX and VDJ file of the SAME
# sample, so matching on the raw token fails:
#   ..._Bac17B_CD3pos_barcodes.tsv.gz   vs  ..._Bac17B_TCR_all_contig_...csv.gz
#
# CD3pos/CD3neg are NOT modality tags -- they identify which sorted library the
# barcodes came from, and a sample can have both. Stripping them pairs a TCR
# file against the CD3-NEGATIVE (non-T-cell) barcodes, which yields ~1% overlap
# and reads as poor TCR recovery rather than as a mismatch. Leave them in the
# token and resolve the sorted case explicitly in pair_up().
MODIFIER = re.compile(
    r"(?:^|[_-])(gex|rna|5gex|tcr|vdj[_-]?t|vdj[_-]?b|vdj|"
    r"bcr|t[_-]?cell|filtered|all)(?=[_-]|$)", re.I)

# GEX libraries that a TCR library belongs with, when the sample has sorted ones.
T_SORTED = ["", "_cd3pos", "_cd3", "_cd3p", "_tcells", "_tcell"]

# BCR files are not useful here and would pair against the same GEX barcodes.
IS_BCR = re.compile(r"vdj[_-]?b|bcr|_ig[hkl]_", re.I)


def pair_up(names: list[str], gse: str) -> list[tuple[str, str, str]]:
    """Match a barcodes file to a T-cell contig file for the same sample."""
    def token(n: str, suffix_re: re.Pattern) -> str:
        s = re.sub(r"^(GSM\d+|GSE\d+)_", "", n)
        s = suffix_re.sub("", s)
        prev = None
        while prev != s:                       # modifiers can stack
            prev = s
            s = MODIFIER.sub("", s)
        return re.sub(r"[_-]+", "_", s).strip("_-.").lower()

    bc_sfx = re.compile(r"_?barcodes\.tsv\.gz$", re.I)
    ct_sfx = re.compile(r"_?(filtered|all)_contig_annotations[^.]*\.csv\.gz$", re.I)

    bcs, ctg = {}, {}
    for n in names:
        if BARCODE_FILE.search(n):
            bcs.setdefault(token(n, bc_sfx), n)
        elif CONTIG_FILE.search(n) and not IS_BCR.search(n):
            ctg.setdefault(token(n, ct_sfx), n)

    out = []
    for tok, c in ctg.items():
        # A TCR library pairs with the sample's T-cell GEX library: the one with
        # no sort tag, or the CD3-POSITIVE one. Never the CD3-negative one.
        for sfx in T_SORTED:
            if tok + sfx in bcs:
                out.append((tok + sfx or "(series)", bcs[tok + sfx], c))
                break
    # Single-sample series: one GEX, one TCR, tokens need not agree.
    if not out and len(bcs) == 1 and len(ctg) == 1:
        out = [("(series)", next(iter(bcs.values())), next(iter(ctg.values())))]
    return out


def read_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt", errors="replace") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def read_contigs(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as fh:
        return pd.read_csv(io.StringIO(fh.read()), low_memory=False)


def suffix_of(barcodes: list[str]) -> str:
    sfx = Counter()
    for b in barcodes[:3000]:
        sfx["-" + b.split("-")[-1] if "-" in b else "(none)"] += 1
    return ", ".join(f"{k}:{v}" for k, v in sfx.most_common(3))


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", nargs="*", default=DEFAULT)
    ap.add_argument("--per-series", type=int, default=1)
    args = ap.parse_args()

    rows = []
    for gse in args.series:
        print(f"\n{'='*92}\n{gse}\n{'='*92}")
        names = filelist(gse)
        if not names:
            print("  no filelist")
            continue
        pairs = pair_up(names, gse)
        print(f"  {len(names)} files, {len(pairs)} barcode/contig pairs matched")
        if not pairs:
            bc = sum(1 for n in names if BARCODE_FILE.search(n))
            ct = sum(1 for n in names if CONTIG_FILE.search(n))
            print(f"  (barcodes files: {bc}, contig files: {ct}) -- "
                  f"naming does not pair up; needs manual look")
            continue

        for tok, bname, cname in pairs[:args.per_series]:
            print(f"\n  sample '{tok}'")
            burl, curl = sample_url(bname, gse), sample_url(cname, gse)
            if not burl or not curl:
                print("    cannot build GSM urls")
                continue
            bp = fetch(burl, DL / gse / bname)
            cp = fetch(curl, DL / gse / cname)
            if not bp or not cp:
                continue

            gex = read_barcodes(bp)
            ctg = read_contigs(cp)
            bc_col = next((c for c in ctg.columns
                           if c.lower() in ("barcode", "cell_id", "cellbarcode")), None)
            if bc_col is None:
                print(f"    no barcode column in contigs: {list(ctg.columns)[:8]}")
                continue

            vdj_all = set(ctg[bc_col].astype(str))
            gex_set = set(gex)
            overlap = gex_set & vdj_all

            # productive TRB only -- that is the actual model input
            prod_col = next((c for c in ctg.columns
                             if c.lower() == "productive"), None)
            chain_col = next((c for c in ctg.columns
                              if c.lower() == "chain"), None)
            cdr3_col = next((c for c in ctg.columns
                             if c.lower() in ("cdr3", "cdr3_aa", "cdr3aa")), None)
            sub = ctg
            if chain_col:
                sub = sub[sub[chain_col].astype(str).str.upper() == "TRB"]
            if prod_col:
                sub = sub[sub[prod_col].astype(str).str.lower().isin(
                    ("true", "t", "yes", "productive"))]
            if cdr3_col:
                sub = sub[sub[cdr3_col].notna()]
            trb = set(sub[bc_col].astype(str)) & gex_set

            print(f"    GEX barcodes            {len(gex_set):>8,}   suffix "
                  f"[{suffix_of(gex)}]")
            print(f"    VDJ barcodes            {len(vdj_all):>8,}   suffix "
                  f"[{suffix_of(sorted(vdj_all))}]")
            print(f"    OVERLAP                 {len(overlap):>8,}   "
                  f"({100*len(overlap)/max(len(gex_set),1):.1f}% of GEX cells)")
            print(f"    with productive TRB     {len(trb):>8,}   "
                  f"({100*len(trb)/max(len(gex_set),1):.1f}% of GEX cells)")
            if cdr3_col and len(sub):
                ex = sub[cdr3_col].dropna().astype(str).head(3).tolist()
                print(f"    example CDR3s           {', '.join(ex)}")
            print(f"    contig columns          {', '.join(list(ctg.columns)[:10])}")

            rows.append(dict(gse=gse, sample=tok, n_gex=len(gex_set),
                             n_vdj_cells=len(vdj_all), overlap=len(overlap),
                             pct_gex_paired=round(100*len(overlap)/max(len(gex_set),1), 1),
                             productive_trb=len(trb),
                             pct_trb=round(100*len(trb)/max(len(gex_set),1), 1)))

    if rows:
        out = pd.DataFrame(rows)
        print(f"\n{'='*92}\nSUMMARY\n{'='*92}")
        with pd.option_context("display.width", 200):
            print(out.to_string(index=False))
        tot = out["productive_trb"].sum()
        print(f"\n  cells usable for training, in these {len(out)} samples: {tot:,}")
        print(f"  median pct of GEX cells with a productive TRB: "
              f"{out['pct_trb'].median():.1f}%")
        out.to_csv(ROOT / "data" / "airr" / "vdj_join_check.tsv",
                   sep="\t", index=False)
    else:
        print("\nnothing checked successfully")


if __name__ == "__main__":
    main()
