"""Is the AS-vs-control signal real, or a global shift?

Trap 5 from the single-cell project: FDR correction filters noise, it cannot
detect a normalisation artefact that moves every gene in one group. A genome-wide
one-directional result is exactly what such an artefact looks like, and also
exactly what a real cell-composition shift looks like. These checks separate them.

  1. Global shift     -- is the whole distribution of g displaced from zero?
  2. Sample loadings  -- does total signal / library size track disease status?
  3. Label permutation -- how many 'hits' does the identical pipeline produce
                          when the disease labels are shuffled within study?
  4. Composition      -- are the hits a cell-type signature, or diffuse?

Usage:  python scripts/bulk_diagnostics.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bulk_meta import (  # noqa: E402
    PREP,
    PRIMARY,
    build_design,
    dersimonian_laird,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "bulk" / "meta"

RNG = np.random.default_rng(0)

# Marker sets for check 4. Deliberately small and unambiguous.
LINEAGE = {
    "T cell": ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK", "TRAC", "IL7R", "CD28"],
    "NK": ["KLRD1", "KLRF1", "NKG7", "GNLY", "PRF1", "EOMES", "KLRB1", "NCR3", "KLRC1"],
    "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "BLNK", "TCL1A", "BANK1"],
    "monocyte": ["CD14", "LYZ", "VCAN", "FCN1", "CSF1R", "CD68", "MNDA", "CEBPB"],
    "neutrophil": ["ELANE", "DEFA3", "DEFA4", "MPO", "CSF3R", "FCGR3B", "CEACAM8", "LCN2"],
    "platelet": ["PF4", "PPBP", "ITGA2B", "GP9", "TUBB1"],
    "erythroid": ["HBB", "HBD", "ALAS2", "AHSP", "SLC4A1"],
}


def load(gse: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(PREP / f"{gse}_expr.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(PREP / f"{gse}_samples.tsv", sep="\t", index_col=0)
    expr = expr[[c for c in expr.columns if c in meta.index]]
    return expr, meta.loc[expr.columns]


def effects_for_labels(expr: pd.DataFrame, meta: pd.DataFrame, labels: np.ndarray):
    """Same estimator as bulk_meta, but with the group vector supplied."""
    m = meta.copy()
    m["group"] = np.where(labels, "AS", "control")
    X, names, _ = build_design(m)
    n, p = X.shape
    dof = n - p
    j = names.index("_AS")

    Y = expr.to_numpy(float).T
    keep = Y.std(axis=0) > 0
    Y, genes = Y[:, keep], expr.index[keep]

    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    mse = ((Y - X @ beta) ** 2).sum(axis=0) / dof
    g = (1.0 - 3.0 / (4.0 * dof - 1.0)) * beta[j] / np.sqrt(mse)
    var_g = np.linalg.inv(X.T @ X)[j, j] + g**2 / (2.0 * dof)
    return pd.Series(g, index=genes), pd.Series(var_g, index=genes)


def pooled_hits(gs: dict, vs: dict, min_k: int = 3) -> pd.DataFrame:
    gw = pd.DataFrame(gs)
    vw = pd.DataFrame(vs)
    gw, vw = gw.align(vw, join="outer")
    keep = gw.notna().sum(axis=1) >= min_k
    gw, vw = gw[keep], vw[keep]
    res = dersimonian_laird(gw.to_numpy().T, vw.to_numpy().T)
    out = pd.DataFrame(res, index=gw.index)
    out["q"] = multipletests(out["p"].fillna(1.0), method="fdr_bh")[1]
    return out


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    data = {g: load(g) for g in PRIMARY}

    # ------------------------------------------------------------------
    print("=" * 100)
    print("CHECK 1 -- global shift: where is the whole distribution of g?")
    print("=" * 100)
    print("  A real signal moves a minority of genes. An artefact moves all of them.")
    print(f"\n  {'study':11s} {'median g':>9s} {'mean g':>8s} {'sd g':>6s} "
          f"{'frac g<0':>9s} {'|g|>0.5':>8s}")
    obs_g, obs_v = {}, {}
    for gse, (expr, meta) in data.items():
        g, v = effects_for_labels(expr, meta, (meta["group"] == "AS").to_numpy())
        obs_g[gse], obs_v[gse] = g, v
        print(f"  {gse:11s} {g.median():9.3f} {g.mean():8.3f} {g.std():6.3f} "
              f"{(g < 0).mean():9.1%} {(g.abs() > 0.5).mean():8.1%}")

    print("\n  Interpretation: frac(g<0) near 50% = no global shift. Far from 50% =")
    print("  every gene is being moved, and FDR will not notice.")

    # ------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("CHECK 2 -- does overall signal level track disease status?")
    print("=" * 100)
    print("  If AS samples are systematically dimmer/brighter overall, every gene")
    print("  inherits that, and the contrast measures the normalisation, not biology.")
    for gse, (expr, meta) in data.items():
        tot = expr.mean(axis=0)
        a = tot[(meta["group"] == "AS").to_numpy()]
        c = tot[(meta["group"] == "control").to_numpy()]
        t, p = stats.ttest_ind(a, c, equal_var=False)
        d = (a.mean() - c.mean()) / np.sqrt((a.var(ddof=1) + c.var(ddof=1)) / 2)
        print(f"  {gse:11s} mean per-sample signal: AS {a.mean():7.3f}  HC {c.mean():7.3f}  "
              f"d = {d:+.2f}  p = {p:.3g}")

    # ------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("CHECK 3 -- permutation null: hits when the labels are shuffled")
    print("=" * 100)
    obs = pooled_hits(obs_g, obs_v)
    n_obs = int((obs["q"] < 0.05).sum())
    print(f"  observed: {n_obs} genes at q<0.05")

    n_perm = 20
    counts = []
    for i in range(n_perm):
        pg, pv = {}, {}
        for gse, (expr, meta) in data.items():
            # Shuffle labels WITHIN study, preserving each study's group sizes.
            lab = (meta["group"] == "AS").to_numpy()
            pg[gse], pv[gse] = effects_for_labels(expr, meta, RNG.permutation(lab))
        counts.append(int((pooled_hits(pg, pv)["q"] < 0.05).sum()))
        print(f"    permutation {i + 1:2d}: {counts[-1]} hits")
    counts = np.array(counts)
    print(f"\n  permuted hits: median {np.median(counts):.0f}, "
          f"max {counts.max()}, mean {counts.mean():.1f}")
    print(f"  observed {n_obs} vs permuted max {counts.max()}  -> "
          f"{'SIGNAL EXCEEDS NULL' if n_obs > counts.max() else 'INDISTINGUISHABLE FROM NULL'}")

    # ------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("CHECK 4 -- are the hits a cell-composition signature?")
    print("=" * 100)
    print("  A shift in blood cell proportions moves whole lineage programmes")
    print("  coherently. That is real biology, but it is composition, not regulation.")
    hits = obs[obs["q"] < 0.05]
    up = set(hits.index[hits["pooled_g"] > 0])
    dn = set(hits.index[hits["pooled_g"] < 0])
    print(f"\n  {len(hits)} hits: {len(up)} up in AS, {len(dn)} down in AS")
    print(f"\n  {'lineage':12s} {'measured':>8s} {'mean g':>7s} {'n up':>5s} {'n down':>7s}  genes down in AS")
    for name, genes in LINEAGE.items():
        have = [g for g in genes if g in obs.index]
        if not have:
            continue
        sub = obs.loc[have]
        d = sorted(g for g in have if g in dn)
        print(f"  {name:12s} {len(have):8d} {sub['pooled_g'].mean():+7.3f} "
              f"{len([g for g in have if g in up]):5d} {len(d):7d}  {','.join(d)}")

    obs.to_csv(OUT / "meta_primary_rerun.tsv", sep="\t")
    print(f"\nwrote {OUT / 'meta_primary_rerun.tsv'}")


if __name__ == "__main__":
    main()
