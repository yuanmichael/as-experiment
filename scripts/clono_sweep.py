"""Open-ended clonotype discovery: which sequences mark AS, not just the one we knew.

The target clonotype `CASSVGLYSTDTQYF` came from prior studies, so finding it
again tests a hypothesis someone else formed. This asks the harder question:
given millions of clonotypes and no prior, which ones separate AS patients from
healthy controls?

Four design decisions, each earned by something that went wrong earlier.

**Screen on expansion, not only presence.** Our own clone would be missed by a
presence sweep: in PRJEB51315 blood it sits in 5 AS and 4 PsA patients, which is
unremarkable. Its entire signal is in frequency. So both sweeps run, and where
they disagree is informative rather than a nuisance.

**Normalise depth per subject.** Library depth spans 17-fold within a group, and
a deeper library shows more distinct clonotypes for purely sampling reasons. Any
presence-based comparison between groups of unequal depth is measuring depth.
Counts are therefore thinned to a common number of junction reads before the
presence sweep. (At the group level there is no confound here -- AS vs healthy
depth gives p = 0.39 -- but per-subject variation is real.)

**Permutation null, always.** On the PRJEB51315 blood arm, clonotypes present in
>=4 AS and 0 PsA patients numbered 6 with true labels and a median of 11 under
shuffled labels: fewer than chance. Without the null that list of 6 would have
looked like a finding. The null is cheap and non-negotiable.

**A healthy-versus-healthy arm, to measure batch rather than argue about it.**
Cohort membership is predictable from repertoires at AUC 0.994, so a sweep that
crosses deposits can return assay artefacts. Running the identical test on
Komech's healthy controls against Emerson's healthy subjects answers how much:
every hit there is batch by construction, because both groups are healthy. If
that list is short while the AS list is long, cross-study discovery is sound. If
they are comparable, only the within-study comparison can be trusted.

What this cannot do: V-gene usage. Read-level extraction recovers the junction
but not the germline V call, so TRBV9 enrichment is not testable from these
tables -- it needs alignment, which is a separate job.

    python scripts/clono_sweep.py --help
"""
import argparse
import collections
import glob
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

EXPANSION_THRESHOLD = 1e-5      # see the note on denominators below
MIN_READS_PRESENT = 2           # a clonotype needs 2 reads to count as present
N_PERMUTATIONS = 200
MIN_CDR3_LEN = 10               # shorter junctions are extraction fragments

# Why the expansion threshold moved from 1e-4 to 1e-5.
#
# Frequencies were originally computed as motif reads / CASS-containing reads,
# taken from the read-level scan. That denominator was wrong: the scan searched
# both strands for the motif but only the forward strand for CASS, and in these
# libraries read 2 carries the CDR3 reverse-complemented (47.2% of R2 reads
# contain CASS on the reverse complement, 0.06% forward). The scan therefore
# counted CASS in 0.72% of reads where the true figure is 24.3%, understating
# the denominator ~34x and inflating every frequency by the same factor -- and
# not by a constant, since the inflation tracks each library's R1/R2 mix.
#
# These tables have an honest denominator (all junction reads; 58.6% of them
# carry a CDR3 beginning CASS, as expected biologically). On that scale no AS
# patient in either cohort reaches 1e-4, so a 1e-4 screen selects nothing.
# 1e-5 is the operating point where the cohorts actually separate -- chosen
# post hoc, so the threshold-free Mann-Whitney is the number to trust.


def load_tables(root: str, group_of=None, subsets=None) -> pd.DataFrame:
    """One row per (subject, clonotype). Replicate runs are summed per subject.

    `subsets` restricts to sorted fractions by label ('F' unsorted, '4' CD4,
    '8' CD8). Summing across fractions mixes denominators -- a CD8 library is
    intrinsically more clonal than whole blood -- so any diversity comparison
    between groups with different fraction mixes needs this.
    """
    meta = {}
    for f in glob.glob(f"{root}/**/*_cdr3.json", recursive=True):
        with open(f) as fh:
            d = json.load(fh)
        if d.get("status") == "ok":
            meta[d["run"]] = d

    frames = []
    for f in sorted(glob.glob(f"{root}/**/*_cdr3.tsv.gz", recursive=True)):
        run = os.path.basename(f).split("_")[0]
        m = meta.get(run)
        if m is None:
            continue
        if subsets and m.get("subset", "") not in subsets:
            continue
        t = pd.read_csv(f, sep="\t")
        # Junctions under 10 aa are extraction fragments, not CDR3s. They are
        # only 0.73% of reads but they are degenerate, so near-identical stubs
        # (CASSLW, CASSLF) turn up in almost every subject and monopolise the
        # top of an expansion sweep. The length histogram is otherwise textbook
        # -- modal 14-15 aa -- so this is the whole of the cleanup needed.
        t = t[t["cdr3"].str.len() >= MIN_CDR3_LEN]
        # `subject` in the 2018 deposit, `patient` in PRJEB51315.
        t["subject"] = m.get("subject") or m.get("patient")
        t["group"] = group_of(m) if group_of else m.get("disease")
        frames.append(t)
    if not frames:
        sys.exit(f"no clonotype tables under {root}")
    df = pd.concat(frames, ignore_index=True)
    return (df.groupby(["group", "subject", "cdr3"], as_index=False)["reads"].sum())


def thin(df: pd.DataFrame, depth: int, rng) -> pd.DataFrame:
    """Subsample every subject to `depth` junction reads, so that distinct-clonotype
    counts are comparable. Subjects below `depth` are dropped, not padded."""
    out = []
    dropped = []
    for (grp, sub), t in df.groupby(["group", "subject"]):
        total = int(t["reads"].sum())
        if total < depth:
            dropped.append((grp, sub, total))
            continue
        counts = rng.multivariate_hypergeometric(t["reads"].to_numpy(), depth)
        keep = t.assign(reads=counts)
        out.append(keep[keep["reads"] > 0])
    if dropped:
        print(f"  dropped {len(dropped)} subjects below depth {depth:,}: "
              + ", ".join(f"{g}/{s} ({n:,})" for g, s, n in dropped[:8]))
    return pd.concat(out, ignore_index=True)


def carrier_matrix(df: pd.DataFrame, mode: str):
    """A boolean subjects-by-clonotypes carrier matrix, plus its labels.

    The sweep is a Fisher test per clonotype, repeated under hundreds of label
    shuffles. Done literally that is hundreds of millions of `fisher_exact`
    calls and it does not finish. Two facts make it cheap instead:

      * With the group sizes fixed, the one-sided p-value for [[a, nc-a],
        [b, nk-b]] depends only on `a` (carriers among cases) and `a+b` (total
        carriers). So every p-value any permutation can ever need comes from one
        small precomputed table -- see `p_table`.
      * `a` for every clonotype at once is a single sparse matrix-vector
        product: an indicator row vector for the case subjects times the
        carrier matrix.

    That turns each shuffle into one matvec over ~7M non-zeros plus an integer
    lookup: milliseconds instead of hours.
    """
    from scipy import sparse

    if mode == "presence":
        sel = df[df["reads"] >= MIN_READS_PRESENT]
    else:
        tot = df.groupby("subject")["reads"].transform("sum")
        sel = df[df["reads"] / tot >= EXPANSION_THRESHOLD]

    subj_codes, subjects = pd.factorize(sel["subject"], sort=True)
    cdr3_codes, cdr3s = pd.factorize(sel["cdr3"], sort=True)
    m = sparse.csr_matrix(
        (np.ones(len(sel), dtype=np.int32), (subj_codes, cdr3_codes)),
        shape=(len(subjects), len(cdr3s)))
    m.data[:] = 1          # a subject carries a clonotype once, however many runs
    return m, np.asarray(subjects), np.asarray(cdr3s)


def p_table(n_case: int, n_ctrl: int) -> np.ndarray:
    """p[m, a] = one-sided Fisher p for `a` case carriers out of `m` total."""
    from scipy.stats import hypergeom

    n = n_case + n_ctrl
    tab = np.ones((n + 1, n_case + 1))
    for m in range(n + 1):
        a = np.arange(0, min(m, n_case) + 1)
        # P(X >= a) for X hypergeometric: m draws from n, n_case of them cases.
        tab[m, a] = hypergeom.sf(a - 1, n, n_case, m)
    return tab


def sweep_matrix(mat, is_case: np.ndarray, min_case: int, ptab: np.ndarray):
    """Carriers-in-cases, total carriers and p for every clonotype at once."""
    a = np.asarray(is_case.astype(np.int32) @ mat).ravel()
    tot = np.asarray(mat.sum(axis=0)).ravel()
    p = ptab[tot, a]
    tested = a >= min_case
    return a, tot, p, tested


def permutation_null(mat, n_case: int, min_case: int, alpha: float,
                     n_perm: int, ptab: np.ndarray, rng) -> np.ndarray:
    """How many clonotypes clear `alpha` when the labels mean nothing."""
    n = mat.shape[0]
    counts = []
    for _ in range(n_perm):
        idx = rng.permutation(n)[:n_case]
        is_case = np.zeros(n, dtype=bool)
        is_case[idx] = True
        a, tot, p, tested = sweep_matrix(mat, is_case, min_case, ptab)
        counts.append(int(((p <= alpha) & tested).sum()))
    return np.array(counts)


def aggregate_features(df: pd.DataFrame) -> pd.DataFrame:
    """Per-subject repertoire summaries. One test each, so these have power at
    n=25 where a per-clonotype sweep does not.

    These MUST be computed on depth-thinned counts. Shannon entropy, clonality
    and distinct-clonotype counts all move with sequencing depth: a shallower
    library looks less diverse and more clonal because fewer rare clones were
    sampled, not because the repertoire differs. On the PRJEB51315 blood arm the
    unthinned version reported AS as markedly more clonal than PsA (top-100
    share 29% vs 13%, p = 0.015) while AS libraries carried 2.8x fewer junction
    reads -- exactly the artefact this docstring exists to prevent.
    """
    rows = []
    for (grp, sub), t in df.groupby(["group", "subject"]):
        r = t["reads"].to_numpy().astype(float)
        f = r / r.sum()
        f_sorted = np.sort(f)[::-1]
        rows.append({
            "group": grp, "subject": sub,
            "junction_reads": int(r.sum()),
            "n_clonotypes": len(r),
            "top1_share": f_sorted[0],
            "top100_share": f_sorted[:100].sum(),
            "shannon": float(-(f * np.log(f)).sum()),
            "simpson_clonality": float((f ** 2).sum()),
            "mean_cdr3_len": float(np.average([len(c) for c in t["cdr3"]], weights=f)),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True, nargs="+",
                    help="directories of *_cdr3.tsv.gz plus *_cdr3.json")
    ap.add_argument("--case", required=True, help="group label treated as cases")
    ap.add_argument("--ctrl", required=True, help="group label treated as controls")
    ap.add_argument("--depth", type=int, default=800_000,
                    help="junction reads to thin each subject to (presence sweep). "
                         "Set this near the SHALLOWEST library, not far below it: "
                         "thinning to 30,000 when the libraries hold 5,000,000 "
                         "hides every clonotype below ~1e-4 and silently "
                         "destroys the sweep's sensitivity.")
    ap.add_argument("--control-cdr3", default=None,
                    help="a clonotype known to separate the groups, reported as a "
                         "positive control. A sweep that cannot recover it is "
                         "underpowered, and its null result means nothing.")
    ap.add_argument("--min-case", type=int, default=4,
                    help="minimum case carriers for a clonotype to be tested")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    ap.add_argument("--out", default=None, help="write ranked hits here (tsv)")
    ap.add_argument("--subsets", nargs="*", default=None,
                    help="restrict to sorted fractions, e.g. F (unsorted) or 8 (CD8)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    df = pd.concat([load_tables(t, subsets=args.subsets) for t in args.tables],
                   ignore_index=True)
    df = df[df["group"].isin([args.case, args.ctrl])]
    counts = df.groupby("group")["subject"].nunique()
    print(f"loaded {len(df):,} (subject, clonotype) rows")
    print(f"  {dict(counts)}  |  {df['cdr3'].nunique():,} distinct clonotypes")

    from scipy.stats import mannwhitneyu

    print(f"\nraw depth per subject (the confounder, reported before anything else)")
    raw = df.groupby(["group", "subject"])["reads"].sum().reset_index()
    for g in (args.case, args.ctrl):
        v = raw.loc[raw["group"] == g, "reads"]
        print(f"  {g:<8} median {int(v.median()):>12,}  min {int(v.min()):>12,}"
              f"  max {int(v.max()):>12,}")
    print("  case vs ctrl depth: p = %.3f"
          % mannwhitneyu(raw.loc[raw["group"] == args.case, "reads"],
                         raw.loc[raw["group"] == args.ctrl, "reads"]).pvalue)

    print(f"\naggregate repertoire features, thinned to {args.depth:,} junction reads")
    thinned = thin(df, args.depth, rng)
    feats = aggregate_features(thinned)
    a = feats[feats["group"] == args.case]
    b = feats[feats["group"] == args.ctrl]
    print(f"  n = {len(a)} cases, {len(b)} controls after thinning")
    print(f"  {'feature':<20}{'case median':>14}{'ctrl median':>14}{'p':>9}")
    for col in ("n_clonotypes", "top1_share", "top100_share",
                "shannon", "simpson_clonality", "mean_cdr3_len"):
        p = mannwhitneyu(a[col], b[col]).pvalue
        print(f"  {col:<20}{a[col].median():>14.4g}{b[col].median():>14.4g}{p:>9.3f}"
              + ("  *" if p < 0.05 else ""))

    for mode in ("presence", "expansion"):
        print(f"\n{'=' * 72}\n{mode.upper()} SWEEP\n{'=' * 72}")
        work = thin(df, args.depth, rng) if mode == "presence" else df
        mat, subjects, cdr3s = carrier_matrix(work, mode)
        grp = (work.drop_duplicates("subject").set_index("subject")["group"]
               .reindex(subjects).to_numpy())
        is_case = grp == args.case
        n_case, n_ctrl = int(is_case.sum()), int((~is_case).sum())
        ptab = p_table(n_case, n_ctrl)
        a, tot, p, tested = sweep_matrix(mat, is_case, args.min_case, ptab)
        print(f"  {n_case} cases, {n_ctrl} controls, {mat.shape[1]:,} clonotypes "
              f"observed, {int(tested.sum()):,} testable at >= {args.min_case} "
              f"case carriers")
        res = pd.DataFrame({
            "cdr3": cdr3s[tested], "case": a[tested], "n_case": n_case,
            "ctrl": tot[tested] - a[tested], "n_ctrl": n_ctrl, "p": p[tested],
        }).sort_values("p")
        hits = res[res["p"] <= args.alpha]
        null = permutation_null(mat, n_case, args.min_case, args.alpha,
                                args.permutations, ptab, rng)
        exp = float(null.mean())
        fdr = exp / max(len(hits), 1)
        print(f"  clonotypes with p <= {args.alpha}: {len(hits):,}")
        print(f"  expected under {args.permutations} label shuffles: "
              f"{exp:.1f} (median {np.median(null):.0f}, 95th pct {np.percentile(null, 95):.0f})")
        print(f"  --> empirical FDR {min(fdr, 1.0):.0%}"
              + ("   NOTHING HERE: observed is at or below chance" if len(hits) <= exp
                 else f"   {len(hits) / max(exp, 1e-9):.1f}x enrichment over chance"))
        if args.control_cdr3:
            c = args.control_cdr3
            row = res[res["cdr3"] == c]
            if row.empty:
                where = np.flatnonzero(cdr3s == c)
                n_car = int(a[where[0]]) if len(where) else 0
                print(f"\n  POSITIVE CONTROL {c}: NOT TESTABLE "
                      f"({n_car} case carriers, below --min-case "
                      f"{args.min_case}). This sweep could not have found it, so "
                      f"a null result above is uninformative.")
            else:
                r = row.iloc[0]
                rank = int((res["p"] < r["p"]).sum()) + 1
                print(f"\n  POSITIVE CONTROL {c}: {int(r['case'])}/{int(r['n_case'])} "
                      f"case, {int(r['ctrl'])}/{int(r['n_ctrl'])} ctrl, p = {r['p']:.4g}, "
                      f"rank {rank:,} of {len(res):,}"
                      + ("  <- recovered" if r["p"] <= args.alpha
                         else "  <- present but does not clear alpha"))
        if len(hits):
            print(f"\n  top hits ({mode}):")
            print(hits.head(20).to_string(index=False))
        if args.out:
            path = args.out.replace(".tsv", f"_{mode}.tsv")
            res.to_csv(path, sep="\t", index=False)
            print(f"\n  wrote {path}")


if __name__ == "__main__":
    main()
