"""Open-ended clonotype sweep in Hanson 2020, with a permutation null.

Every clonotype in the sorted-CD8 compartment is tested for enrichment in AS
against the cohort's own controls. Nothing is pre-filtered by a candidate list,
and the carrier threshold is a reporting parameter rather than a significance
floor -- the permutation null decides what is real at any threshold.

Three controls, as the project's rules require:
  * a permutation null over patient labels, so the discovery count at every
    threshold is read against what shuffling produces;
  * the CD4 fraction, swept identically, where a CD8-restricted signal must
    not appear;
  * the motif family, which must be recovered or the sweep is uninformative.

The p-value trick from clono_sweep applies here too: with group sizes fixed,
the one-sided Fisher p for a clonotype depends only on its carrier count among
cases and its total carrier count, so one precomputed table serves every
permutation and each shuffle is a sparse matrix-vector product.

    python scripts/hanson_sweep.py --export <dir> --min-carriers 3
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clono_sweep import p_table  # noqa: E402
from hanson_scan import parse_name  # noqa: E402
from motif_family_test import PRESPECIFIED  # noqa: E402

SEED = 0


def load_fraction(export: str, fraction: str, min_len: int) -> pd.DataFrame:
    rows = []
    files = sorted(f for f in os.listdir(export) if f.endswith(".tsv"))
    for k, fn in enumerate(files, 1):
        meta = parse_name(os.path.splitext(fn)[0])
        if meta is None or meta["prefix"] not in ("A", "K"):
            continue
        if meta["fraction"] != fraction:
            continue
        path = os.path.join(export, fn)
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            i_aa, i_n = hdr.index("aminoAcid"), hdr.index("count (templates/reads)")
            i_st = hdr.index("sequenceStatus")
            counts: dict[str, int] = {}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) <= i_st or p[i_st] != "In":
                    continue
                aa = p[i_aa]
                if len(aa) < min_len:
                    continue
                try:
                    counts[aa] = counts.get(aa, 0) + int(p[i_n])
                except ValueError:
                    pass
        rows.append(pd.DataFrame({
            "subject": meta["subject"],
            "group": "as" if meta["prefix"] == "A" else "healthy",
            "cdr3": list(counts), "reads": list(counts.values())}))
        if k % 50 == 0:
            print(f"    {k}/{len(files)}")
    return pd.concat(rows, ignore_index=True)


def thin_to(df: pd.DataFrame, depth: int, rng) -> pd.DataFrame:
    """Subsample every subject to `depth` templates. `depth <= 0` skips it.

    Thinning is mandatory where library depth differs by group, because
    presence is depth-sensitive and label shuffling would then manufacture
    enrichment. Here CD8 depth does not differ (median 66K AS vs 62K control,
    Mann-Whitney p = 0.55), so thinning costs 14 of 48 AS patients to remove a
    confound that is not present. Both settings are reported.
    """
    if depth <= 0:
        print(f"  no thinning; depths {df.groupby('subject')['reads'].sum().min():,}"
              f"-{df.groupby('subject')['reads'].sum().max():,} templates")
        return df
    out, dropped = [], []
    for sub, t in df.groupby("subject"):
        total = int(t["reads"].sum())
        if total < depth:
            dropped.append((sub, total))
            continue
        c = rng.multivariate_hypergeometric(t["reads"].to_numpy(), depth)
        keep = t.assign(reads=c)
        out.append(keep[keep["reads"] > 0])
    if dropped:
        print(f"  dropped {len(dropped)} below {depth:,} templates: "
              + ", ".join(f"{s}({n:,})" for s, n in dropped[:6]))
    return pd.concat(out, ignore_index=True)


def matrix(df: pd.DataFrame):
    sc, subjects = pd.factorize(df["subject"], sort=True)
    cc, cdr3s = pd.factorize(df["cdr3"], sort=True)
    m = sparse.csr_matrix((np.ones(len(df), np.int32), (sc, cc)),
                          shape=(len(subjects), len(cdr3s)))
    m.data[:] = 1
    return m, np.asarray(subjects), np.asarray(cdr3s)


def sweep(df: pd.DataFrame, label: str, min_car: int, n_perm: int, rng):
    m, subjects, cdr3s = matrix(df)
    grp = (df.drop_duplicates("subject").set_index("subject")["group"]
           .reindex(subjects).to_numpy())
    is_case = grp == "as"
    n_case, n_ctrl = int(is_case.sum()), int((~is_case).sum())
    tot = np.asarray(m.sum(0)).ravel()
    a = np.asarray(is_case.astype(np.int32) @ m).ravel()
    tab = p_table(n_case, n_ctrl)
    p = tab[tot, a]
    print(f"\n{label}: {n_case} AS, {n_ctrl} healthy, {m.shape[1]:,} clonotypes")

    null = np.zeros((n_perm, 3), int)
    idx = np.arange(len(is_case))
    for i in range(n_perm):
        rng.shuffle(idx)
        sh = is_case[idx]
        aa = np.asarray(sh.astype(np.int32) @ m).ravel()
        pp = tab[tot, aa]
        null[i] = [(pp <= 0.05).sum(), (pp <= 0.01).sum(), (pp <= 0.001).sum()]
    obs = [(p <= 0.05).sum(), (p <= 0.01).sum(), (p <= 0.001).sum()]
    print(f"  {'alpha':>8}{'observed':>11}{'null mean':>12}{'null 95th':>12}"
          f"{'empirical p':>14}")
    for j, al in enumerate((0.05, 0.01, 0.001)):
        emp = (1 + (null[:, j] >= obs[j]).sum()) / (n_perm + 1)
        print(f"  {al:>8}{obs[j]:>11,}{null[:, j].mean():>12,.1f}"
              f"{np.percentile(null[:, j], 95):>12,.1f}{emp:>14.3f}")

    keep = (a >= min_car) & (tot == a)
    res = pd.DataFrame({"cdr3": cdr3s[keep], "as_carriers": a[keep],
                        "ctrl_carriers": (tot - a)[keep], "p": p[keep]})
    print(f"  clonotypes in >= {min_car} AS and 0 controls: {len(res):,}")
    fam = np.isin(cdr3s, list(PRESPECIFIED))
    if fam.any():
        print("  positive control (motif family):")
        for c, ac, pc in zip(cdr3s[fam], a[fam], p[fam]):
            print(f"    {c:<20} {ac:>2}/{n_case} AS  p = {pc:.4g}")
    return res.sort_values("p"), n_case, n_ctrl


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", required=True)
    ap.add_argument("--min-carriers", type=int, default=3)
    ap.add_argument("--depth", type=int, default=40_000)
    ap.add_argument("--min-len", type=int, default=10)
    ap.add_argument("--permutations", type=int, default=200)
    ap.add_argument("--out", default="data/panel/hanson_sweep.tsv")
    args = ap.parse_args()

    rng = np.random.default_rng(SEED)
    print("loading CD8")
    cd8 = thin_to(load_fraction(args.export, "CD8", args.min_len), args.depth, rng)
    res8, n_case, n_ctrl = sweep(cd8, "CD8", args.min_carriers,
                                 args.permutations, rng)
    print("\nloading CD4 (negative control)")
    cd4 = thin_to(load_fraction(args.export, "CD4", args.min_len), args.depth, rng)
    sweep(cd4, "CD4", args.min_carriers, args.permutations, rng)

    print(f"\nTOP CD8 HITS (>= {args.min_carriers} AS carriers, 0 controls)")
    print(f"  {'cdr3':<22}{'AS':>5}/{n_case:<4}{'p':>10}")
    for _, r in res8.head(20).iterrows():
        fam = " *family*" if r["cdr3"] in PRESPECIFIED else ""
        print(f"  {r['cdr3']:<22}{r['as_carriers']:>5}/{n_case:<4}{r['p']:>10.4g}{fam}")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    res8.to_csv(args.out, sep="\t", index=False, float_format="%.6g")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
