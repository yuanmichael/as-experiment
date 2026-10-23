"""Are there CDR3s shared by AS patients and absent from healthy controls?

This is the question the project started with, and it is Faham et al.'s
statistic: not a model, but a 2x2 table per sequence -- present/absent against
case/control -- and the share of patients carrying it.

Data: GSE288581, which deposits 10x `filtered_contig_annotations.csv` files, so
no FASTQ processing is needed. All samples are the same sort gate (CD45RO+ CD8+
memory T cells), which matters: the comparison is like-for-like, but these are
not whole repertoires and frequencies here are not comparable to whole-blood
TCRb studies.

    4 AS patients, blood (PB) and synovial fluid (SF) each
    4 healthy controls, blood

Four against four cannot produce a percentage. With 4 patients the only
observable rates are 0, 25, 50, 75, 100%, and 1-of-4 against 0-of-4 is p = 1.0
by Fisher's exact test. So the within-study test is descriptive.

What makes it worth running anyway is the second control: 786 healthy blood
repertoires from Emerson 2017 are already on disk. Faham's headline rested on
"0 of 227 controls". Any candidate sequence can be checked against 786 healthy
people here for free. That control has its own caveat -- Emerson is whole-blood
genomic DNA TCRb, not sorted CD8 memory, and a different platform -- so a
sequence missing from it could be missing for technical reasons. It is read as
a *negative* control: a candidate that shows up frequently in 786 healthy
people is dead, while absence is suggestive rather than proof.

Usage:
    python scripts/as_clonotypes.py                 # download + analyse
    python scripts/as_clonotypes.py --skip-download
"""

from __future__ import annotations

import argparse
import gzip
import io
import sys
import urllib.request
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "as_vdj"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{stem}nnn/{gsm}/suppl/{fname}"

# GSM id -> (sample name, group, tissue). From the series' own sample titles;
# every one is CD45RO+ CD8+ sorted.
SAMPLES = {
    "GSM8770237": ("SF73-4", "AS", "synovial"),
    "GSM8770239": ("SF76", "AS", "synovial"),
    "GSM8770241": ("SF74", "AS", "synovial"),
    "GSM8770243": ("SF81", "AS", "synovial"),
    "GSM8770245": ("PB73-4", "AS", "blood"),
    "GSM8770247": ("PB76", "AS", "blood"),
    "GSM8770249": ("PB74", "AS", "blood"),
    "GSM8770251": ("PB81", "AS", "blood"),
    "GSM8770253": ("HC1564", "HC", "blood"),
    "GSM8770255": ("HC546", "HC", "blood"),
    "GSM8770257": ("HC1660", "HC", "blood"),
    "GSM8770259": ("HC1788", "HC", "blood"),
    # Same three patients again, sorted further on exhaustion markers. Kept
    # out of the case/control test: they are not independent people, and
    # counting them as extra patients would be pseudo-replication -- the error
    # that halved this project's cell-level result once already.
    "GSM8770261": ("73-4_exh", "AS_exhausted", "synovial"),
    "GSM8770263": ("76_exh", "AS_exhausted", "synovial"),
    "GSM8770265": ("81_exh", "AS_exhausted", "synovial"),
}

# Which patient each sample came from, so a clone is never counted twice.
SUBJECT = {
    "SF73-4": "73-4", "PB73-4": "73-4", "73-4_exh": "73-4",
    "SF74": "74", "PB74": "74",
    "SF76": "76", "PB76": "76", "76_exh": "76",
    "SF81": "81", "PB81": "81", "81_exh": "81",
    "HC1564": "HC1564", "HC546": "HC546", "HC1660": "HC1660", "HC1788": "HC1788",
}


def download(raw: Path) -> None:
    raw.mkdir(parents=True, exist_ok=True)
    listing = urllib.request.urlopen(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288581/suppl/filelist.txt",
        timeout=120).read().decode("utf-8", "replace")
    # GEO's filelist is Archive/File | Name | Time | Size | Type -- the name is
    # the SECOND field. Reading the first once gave every file the name "File".
    names = [ln.split("\t")[1] for ln in listing.splitlines()[1:] if "\t" in ln]
    want = {n.split("_")[0]: n for n in names if "contig_annotations" in n}
    for gsm, (name, _, _) in SAMPLES.items():
        fname = want.get(gsm)
        if not fname:
            print(f"  {gsm} ({name}): no contig file in listing")
            continue
        dest = raw / f"{name}.csv.gz"
        if dest.exists():
            continue
        url = GEO.format(stem=gsm[:-3], gsm=gsm, fname=fname)
        urllib.request.urlretrieve(url, dest)
        print(f"  {name}: {dest.stat().st_size / 1e6:.1f} MB")


def load_sample(path: Path) -> pd.DataFrame:
    """One row per cell-level TRB clonotype in one sample.

    A 10x contig file has one row per assembled contig, so a cell with both
    chains appears twice and a cell can carry two beta contigs. Only
    productive, high-confidence TRB rows with a called CDR3 are kept, and the
    result is deduplicated to one row per (barcode, CDR3) so a cell counts once.
    """
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh)
    cols = {c.lower(): c for c in df.columns}
    chain = cols.get("chain")
    df = df[df[chain].astype(str).str.upper() == "TRB"]
    for flag in ("productive", "high_confidence", "is_cell"):
        if flag in cols:
            df = df[df[cols[flag]].astype(str).str.lower().isin(("true", "yes", "1"))]
    cdr3 = cols.get("cdr3")
    df = df[df[cdr3].notna() & (df[cdr3].astype(str).str.len() >= 6)]
    df = df[~df[cdr3].astype(str).str.contains(r"[*_]", regex=True)]
    keep = df[[cols["barcode"], cdr3, cols.get("v_gene", cdr3), cols.get("j_gene", cdr3)]]
    keep.columns = ["barcode", "cdr3", "v_gene", "j_gene"]
    return keep.drop_duplicates(["barcode", "cdr3"])


def build_table(raw: Path) -> pd.DataFrame:
    rows = []
    for gsm, (name, group, tissue) in SAMPLES.items():
        p = raw / f"{name}.csv.gz"
        if not p.exists():
            continue
        d = load_sample(p)
        clones = (d.groupby("cdr3")
                    .agg(cells=("barcode", "nunique"),
                         v_gene=("v_gene", "first"), j_gene=("j_gene", "first"))
                    .reset_index())
        clones["sample"] = name
        clones["group"] = group
        clones["tissue"] = tissue
        clones["subject"] = SUBJECT[name]
        rows.append(clones)
        print(f"  {name:<10} {group:<13} {tissue:<9} "
              f"{d['barcode'].nunique():>7,} cells  {len(clones):>7,} clonotypes")
    return pd.concat(rows, ignore_index=True)


def case_control(tab: pd.DataFrame, case_tissue: str) -> pd.DataFrame:
    """Per-CDR3 2x2 table: carried / not carried, against AS / healthy.

    The unit is the *subject*, never the cell. A clone expanded to 5,000 cells
    in one patient is still one patient carrying it.
    """
    cases = tab[(tab["group"] == "AS") & (tab["tissue"] == case_tissue)]
    ctrls = tab[(tab["group"] == "HC")]
    n_case = cases["subject"].nunique()
    n_ctrl = ctrls["subject"].nunique()

    c_carriers = cases.groupby("cdr3")["subject"].nunique()
    h_carriers = ctrls.groupby("cdr3")["subject"].nunique()
    allc = c_carriers.index.union(h_carriers.index)
    df = pd.DataFrame({
        "cdr3": allc,
        "as_carriers": c_carriers.reindex(allc).fillna(0).astype(int),
        "hc_carriers": h_carriers.reindex(allc).fillna(0).astype(int),
    }).reset_index(drop=True)
    df["n_as"] = n_case
    df["n_hc"] = n_ctrl
    df["as_pct"] = 100 * df["as_carriers"] / n_case
    df["hc_pct"] = 100 * df["hc_carriers"] / n_ctrl
    # Only sequences seen in more than one patient can distinguish anything,
    # and testing all of them would be millions of tests for no gain.
    cand = df[df["as_carriers"] >= 2].copy()
    cand["fisher_p"] = [
        fisher_exact([[a, n_case - a], [h, n_ctrl - h]], alternative="greater")[1]
        for a, h in zip(cand["as_carriers"], cand["hc_carriers"])
    ]
    return cand.sort_values(["as_carriers", "hc_carriers"], ascending=[False, True])


def check_against_emerson(cdr3s: list[str]) -> pd.DataFrame:
    """How many of 786 healthy blood repertoires carry each candidate?

    This is the control Faham's claim turns on, with 786 people instead of 227.
    Caveat, and it is a real one: Emerson is whole-blood genomic DNA TCRb while
    these are sorted CD8 memory cells on a different platform. A candidate that
    is common in 786 healthy people is refuted; one that is absent is
    encouraging but could be absent for technical reasons.
    """
    import duckdb

    reps = ROOT / "data" / "airr" / "emerson-2017-natgen" / "reps"
    if not reps.exists():
        print("  Emerson repertoires not on disk; skipping the 786-control check")
        return pd.DataFrame()
    glob = str(reps / "*.parquet").replace("\\", "/")
    con = duckdb.connect()
    con.execute("SET memory_limit='4GB'")
    con.register("cand", pd.DataFrame({"cdr3": cdr3s}))
    return con.execute(f"""
        SELECT c.cdr3,
               count(DISTINCT parse_filename(filename, true)) AS healthy_carriers
        FROM cand c
        LEFT JOIN read_parquet('{glob}', filename=true) r
               ON r.junction_aa = c.cdr3
        GROUP BY c.cdr3
    """).df()


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-download", action="store_true")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    raw = OUT / "raw"
    if not args.skip_download:
        print("downloading GSE288581 clonotype tables")
        download(raw)

    print("\nper-sample recovery")
    tab = build_table(raw)
    OUT.mkdir(parents=True, exist_ok=True)
    tab.to_parquet(OUT / "clonotypes.parquet", index=False)

    for tissue in ("blood", "synovial"):
        print(f"\n{'=' * 78}\nAS {tissue} vs healthy blood\n{'=' * 78}")
        if tissue == "synovial":
            print("NOTE: different tissue from the controls, so any difference\n"
                  "      confounds disease with compartment. Read as exploratory.")
        cand = case_control(tab, tissue)
        if cand.empty:
            print("  no CDR3 shared by 2 or more patients")
            continue
        n_as = int(cand["n_as"].iloc[0])
        print(f"  {len(cand):,} CDR3s carried by >=2 of {n_as} AS patients")
        for k in range(n_as, 1, -1):
            hits = cand[(cand["as_carriers"] == k) & (cand["hc_carriers"] == 0)]
            print(f"    in {k}/{n_as} patients and 0/4 controls: {len(hits):,}")

        shared = cand[(cand["as_carriers"] >= max(2, n_as - 1))
                      & (cand["hc_carriers"] == 0)]
        if shared.empty:
            continue
        print(f"\n  checking {len(shared):,} candidates against 786 healthy repertoires")
        em = check_against_emerson(shared["cdr3"].tolist())
        if not em.empty:
            shared = shared.merge(em, on="cdr3", how="left")
            shared["healthy_carriers"] = shared["healthy_carriers"].fillna(0).astype(int)
            shared["healthy_pct"] = 100 * shared["healthy_carriers"] / 786
            shared = shared.sort_values(["as_carriers", "healthy_carriers"],
                                        ascending=[False, True])
        shared.to_csv(OUT / f"candidates_{tissue}.tsv", sep="\t", index=False)
        cols = ["cdr3", "as_carriers", "as_pct", "hc_carriers", "fisher_p"]
        if "healthy_carriers" in shared:
            cols += ["healthy_carriers", "healthy_pct"]
        print(shared[cols].head(args.top).to_string(index=False))
        print(f"  wrote {OUT / f'candidates_{tissue}.tsv'}")

    # Komech reported a TRBV9-associated clonal cluster in AS synovium, and
    # GSE277117 could not test it because 3' chemistry discards V segments.
    # This data has V calls, so it can.
    print(f"\n{'=' * 78}\nTRBV9 usage (Komech's AS-associated segment)\n{'=' * 78}")
    v = (tab.assign(v=tab["v_gene"].astype(str).str.split("*").str[0])
            .groupby(["group", "tissue", "sample"])
            .apply(lambda d: 100 * (d["v"] == "TRBV9").mean(), include_groups=False)
            .rename("trbv9_pct_of_clonotypes").reset_index())
    print(v.to_string(index=False))


if __name__ == "__main__":
    main()
