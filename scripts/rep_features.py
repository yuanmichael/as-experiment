"""Patient-level features from repertoires. One row per subject.

Two blocks, kept apart on purpose.

**Summary features** (this script, pass A) are unsupervised: they never see a
label, so computing them once over everything is safe. V and J gene usage, the
CDR3 length distribution, amino-acid composition, and clonality.

**Public-clone incidence** (pass B) is the raw material for the feature family
that actually carries the signal in Emerson's paper: which specific CDR3s a
subject has. Turning that into a feature *requires* choosing which CDR3s to
look at, and any selection that sees the labels is leakage. So this script
stops at an unsupervised, label-blind object -- a subjects x CDR3 presence
matrix restricted to CDR3s seen in at least `--min-subjects` subjects -- and
the label-dependent selection happens inside each cross-validation fold in
`rep_classify.py`. The `--min-subjects` filter is itself label-blind, which is
what makes doing it once legitimate.

Depth matters. Clone counts run 36,000 to 428,000 across these subjects, and
anything that scales with depth (richness, number of public clones hit) is
partly a sequencing-depth feature. Depth-sensitive features are flagged by name
so a run can exclude them.

Usage:
    python scripts/rep_features.py                      # both passes
    python scripts/rep_features.py --pass summary
    python scripts/rep_features.py --pass public --min-subjects 7
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
SCRATCH = Path(r"C:/Users/myuan/AppData/Local/Temp/claude/C--Coding-as-experiment"
               r"/2f391c41-410f-452c-9331-4e63aba51115/scratchpad")

AA = list("ACDEFGHIKLMNPQRSTVWY")
LEN_BINS = list(range(8, 25))  # CDR3 amino-acid lengths; outside this is lumped

# Features whose value scales with how deeply the sample was sequenced.
DEPTH_SENSITIVE = ("n_clones", "n_templates", "n_unique_cdr3", "log_n_clones",
                   "log_n_templates", "shannon")


# --------------------------------------------------------------------------
# pass A: summary features


def clonality(counts: np.ndarray) -> dict:
    """Shape of the clone-size distribution.

    Raw Shannon entropy rises with the number of clones observed, so it is
    reported both raw (depth-sensitive) and divided by log(n) (normalised, and
    the one to trust across samples of different depth). Clonality is
    1 - normalised entropy: 0 is a perfectly even repertoire, 1 is one clone.
    """
    total = counts.sum()
    p = counts / total
    n = len(counts)
    shannon = float(-(p * np.log(p)).sum())
    norm = shannon / np.log(n) if n > 1 else 0.0
    order = np.sort(counts)[::-1]
    cum = np.cumsum(order) / total
    # Gini over clone sizes: 0 even, 1 maximally unequal.
    asc = np.sort(counts)
    idx = np.arange(1, n + 1)
    gini = float((2 * idx - n - 1).dot(asc) / (n * asc.sum()))
    return dict(
        shannon=shannon,
        shannon_norm=norm,
        clonality=1.0 - norm,
        gini=gini,
        top1_share=float(cum[0]),
        top10_share=float(cum[min(9, n - 1)]),
        top100_share=float(cum[min(99, n - 1)]),
        n_expanded_gt10=int((counts > 10).sum()) / n,
        singleton_frac=float((counts == 1).mean()),
    )


def summarise(path: Path) -> dict:
    df = pq.read_table(path, columns=["v_gene", "j_gene", "junction_aa", "count"]).to_pandas()
    counts = df["count"].to_numpy(dtype=np.float64)
    cdr3 = df["junction_aa"]
    n = len(df)

    row = dict(subject_id=path.stem,
               n_clones=n,
               n_templates=float(counts.sum()),
               n_unique_cdr3=int(cdr3.nunique()))
    row["log_n_clones"] = float(np.log10(n))
    row["log_n_templates"] = float(np.log10(row["n_templates"]))
    row.update(clonality(counts))

    # Gene usage as a fraction of unique clones, not of templates: template
    # weighting lets a single huge clone dominate a subject's V profile.
    #
    # Only TRB genes get their own column. A few clones per subject carry a
    # stray TRAJ or TRDJ call from the aligner; left alone those become columns
    # that exist for one subject and not the next, so the feature matrix would
    # depend on which subjects happen to be in it. They are lumped into
    # `:other` instead, which is itself a mild alignment-quality feature.
    for prefix, col, locus in (("v", "v_gene", "TRBV"), ("j", "j_gene", "TRBJ")):
        frac = df[col].value_counts(normalize=True)
        is_trb = frac.index.str.startswith(locus)
        for gene_name, val in frac[is_trb].items():
            row[f"{prefix}:{gene_name}"] = float(val)
        row[f"{prefix}:other"] = float(frac[~is_trb].sum())

    lengths = cdr3.str.len().to_numpy()
    for L in LEN_BINS:
        row[f"len:{L}"] = float((lengths == L).mean())
    row["len:other"] = float(((lengths < LEN_BINS[0]) | (lengths > LEN_BINS[-1])).mean())
    row["len_mean"] = float(lengths.mean())
    row["len_sd"] = float(lengths.std())

    # Amino-acid composition over the whole junction, each subject's own
    # residue pool as the denominator. The conserved leading C and trailing
    # F/W are left in; they are near-constant and so carry no variance.
    joined = "".join(cdr3.to_numpy(dtype=str))
    arr = np.frombuffer(joined.encode("ascii", "ignore"), dtype=np.uint8)
    tally = np.bincount(arr, minlength=128).astype(np.float64)
    denom = tally[[ord(a) for a in AA]].sum()
    for a in AA:
        row[f"aa:{a}"] = float(tally[ord(a)] / denom) if denom else 0.0
    return row


def pass_summary(reps_dir: Path, out: Path) -> None:
    files = sorted(reps_dir.glob("*.parquet"))
    print(f"pass A: summary features over {len(files)} repertoires")
    rows, t0 = [], time.time()
    for i, f in enumerate(files, 1):
        rows.append(summarise(f))
        if i % 100 == 0 or i == len(files):
            print(f"  {i}/{len(files)}  {time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows).fillna(0.0)
    # Keep a stable column order: identity, then scalars, then each block.
    front = ["subject_id"]
    blocks = (sorted(c for c in df.columns if ":" not in c and c != "subject_id")
              + sorted(c for c in df.columns if c.startswith("v:"))
              + sorted(c for c in df.columns if c.startswith("j:"))
              + sorted(c for c in df.columns if c.startswith("len:"))
              + sorted(c for c in df.columns if c.startswith("aa:")))
    df = df[front + blocks]
    df.to_parquet(out, index=False)
    print(f"  wrote {out}: {df.shape[0]} subjects x {df.shape[1] - 1} features")
    print(f"  depth spread: {df['n_clones'].min():,.0f} to {df['n_clones'].max():,.0f} clones "
          f"(median {df['n_clones'].median():,.0f})")


# --------------------------------------------------------------------------
# pass B: public-clone incidence


def pass_public(reps_dir: Path, outdir: Path, min_subjects: int, key: str) -> None:
    """Subjects x CDR3 presence matrix, restricted to shared CDR3s.

    Done in duckdb because the full table is ~150M rows and a Python dict over
    ~100M distinct CDR3s would not fit in 15 GB. Two passes: count how many
    subjects carry each CDR3, keep those at or above the threshold, then emit
    the (subject, CDR3) pairs for the survivors.
    """
    import duckdb

    glob = str(reps_dir / "*.parquet").replace("\\", "/")
    tmp = SCRATCH / "duckdb_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("SET memory_limit='6GB'")
    con.execute(f"SET temp_directory='{str(tmp).replace(chr(92), '/')}'")
    con.execute("SET preserve_insertion_order=false")

    # The subject is the file stem; duckdb's filename column gives the path.
    keyexpr = {"cdr3": "junction_aa",
               "vcdr3j": "v_gene || '|' || junction_aa || '|' || j_gene"}[key]
    # parse_filename(path, true) strips directory and extension, and handles
    # Windows backslashes -- a hand-rolled regex on '/' silently matched the
    # whole path and produced subject ids nothing could be joined on.
    con.execute(f"""
        CREATE VIEW reps AS
        SELECT parse_filename(filename, true) AS subject_id,
               {keyexpr} AS clone
        FROM read_parquet('{glob}', filename=true)
    """)

    t0 = time.time()
    n_rows = con.execute("SELECT count(*) FROM reps").fetchone()[0]
    print(f"pass B: {n_rows:,} clone rows, {time.time() - t0:.0f}s to count")

    t0 = time.time()
    con.execute(f"""
        CREATE TABLE shared AS
        SELECT clone, count(DISTINCT subject_id) AS n_subj
        FROM reps GROUP BY clone
        HAVING count(DISTINCT subject_id) >= {min_subjects}
    """)
    n_shared = con.execute("SELECT count(*) FROM shared").fetchone()[0]
    print(f"  {n_shared:,} clones in >= {min_subjects} subjects  ({time.time() - t0:.0f}s)")

    t0 = time.time()
    pairs = con.execute("""
        SELECT r.subject_id, r.clone
        FROM reps r JOIN shared s USING (clone)
        GROUP BY 1, 2
    """).fetch_arrow_table()
    print(f"  {pairs.num_rows:,} (subject, clone) pairs  ({time.time() - t0:.0f}s)")

    subj = pd.Index(sorted(p.stem for p in reps_dir.glob("*.parquet")))
    clones = pd.Index(con.execute("SELECT clone FROM shared ORDER BY clone").df()["clone"])
    si = subj.get_indexer(pairs.column("subject_id").to_pandas())
    ci = clones.get_indexer(pairs.column("clone").to_pandas())
    assert (si >= 0).all() and (ci >= 0).all(), "unmapped subject or clone"

    M = sp.csr_matrix((np.ones(len(si), dtype=np.int8), (si, ci)),
                      shape=(len(subj), len(clones)))
    outdir.mkdir(parents=True, exist_ok=True)
    sp.save_npz(outdir / f"public_{key}_min{min_subjects}.npz", M)
    pd.Series(subj, name="subject_id").to_csv(
        outdir / f"public_{key}_min{min_subjects}.subjects.tsv", sep="\t", index=False)
    pd.Series(clones, name="clone").to_frame().to_parquet(
        outdir / f"public_{key}_min{min_subjects}.clones.parquet", index=False)
    dens = M.nnz / (M.shape[0] * M.shape[1])
    print(f"  wrote {M.shape[0]} x {M.shape[1]} incidence matrix, density {dens:.4f}")
    print(f"  clones per subject: median {np.median(M.sum(axis=1).A1):,.0f}")


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="emerson-2017-natgen")
    ap.add_argument("--pass", dest="which", default="both",
                    choices=["summary", "public", "both"])
    ap.add_argument("--min-subjects", type=int, default=7,
                    help="a CDR3 must appear in this many subjects to be kept")
    ap.add_argument("--key", default="cdr3", choices=["cdr3", "vcdr3j"],
                    help="what counts as the same clone across people")
    args = ap.parse_args()

    base = ROOT / "data" / "airr" / args.study
    reps = base / "reps"
    if not reps.exists():
        print(f"no repertoires at {reps} -- run adc_fetch.py first")
        return

    if args.which in ("summary", "both"):
        pass_summary(reps, base / "features_summary.parquet")
    if args.which in ("public", "both"):
        pass_public(reps, base, args.min_subjects, args.key)


if __name__ == "__main__":
    main()
