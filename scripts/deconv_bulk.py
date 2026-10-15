"""Composition or regulation? Deconvolve the bulk cohorts and test both.

The NK transcriptional programme is down in AS blood. Two explanations:

  composition  -- there are fewer NK cells in the sample
  regulation   -- there are the same NK cells, each making less cytotoxic cargo

A published flow-cytometry meta-analysis (n ~ 7,000) reports NK PERCENTAGES are
unchanged, which favours regulation. This script tests it from the expression
side.

The trap: deconvolution reads proportions out of the same transcripts as the
expression signature, so "adjust the signature for the estimated proportion"
is circular by construction. The fix used here is DISJOINT GENE SETS --

  proportion  is estimated from NK IDENTITY genes (receptors, lineage TFs),
              with every effector gene removed from the signature matrix
  the outcome is the NK EFFECTOR programme (granzymes, perforin, granulysin)

so the question becomes: given how many NK cells the identity markers say are
there, is the cytotoxic cargo lower than it should be?

Usage:  python scripts/deconv_bulk.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bulk_meta import PREP, PRIMARY, build_design, dersimonian_laird  # noqa: E402
from deconv_validate import COARSE, GRANULOCYTE, deconvolve  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEC = ROOT / "data" / "deconv"
OUT = ROOT / "data" / "bulk" / "meta"

# The cytotoxic cargo itself -- what an NK cell makes, not what it is.
EFFECTOR = [
    "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "NKG7",
    "CTSW", "SERPINB9", "FASLG",
]
# What an NK cell is, regardless of how armed it is.
IDENTITY = [
    "KLRD1", "KLRF1", "NCR1", "NCR3", "EOMES", "TBX21", "IL2RB", "CD247",
    "S1PR5", "KLRC1", "PLEKHF1", "SH2D1B", "PRSS23", "CD160", "KIR2DL3",
    "FCGR3A", "TXK", "CLIC3",
]

LINEAGES = ["NK", "monocyte", "B", "CD4 T", "CD8 T", "DC"]
RELIABLE = {"NK", "monocyte", "B", "CD4 T"}  # r > 0.8 in validation


def to_linear(expr: pd.DataFrame, gse: str) -> pd.DataFrame:
    """Prepared matrices are log2; deconvolution needs linear scale."""
    lin = np.power(2.0, expr) - 1.0
    return lin.clip(lower=0)


def load_bulk(gse: str):
    expr = pd.read_csv(PREP / f"{gse}_expr.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(PREP / f"{gse}_samples.tsv", sep="\t", index_col=0)
    expr = expr[[c for c in expr.columns if c in meta.index]]
    return expr, meta.loc[expr.columns]


def effect_of(values: pd.Series, meta: pd.DataFrame) -> tuple[float, float]:
    """Standardised AS-vs-control effect for one per-sample measurement,
    using the same covariate-adjusted estimator as the gene-level analysis."""
    X, names, idx = build_design(meta)
    y = values.loc[idx].to_numpy(float)
    n, p = X.shape
    dof = n - p
    j = names.index("_AS")
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    s = np.sqrt((resid**2).sum() / dof)
    if s <= 0:
        return np.nan, np.nan
    g = (1.0 - 3.0 / (4.0 * dof - 1.0)) * beta[j] / s
    var = np.linalg.inv(X.T @ X)[j, j] + g**2 / (2.0 * dof)
    return float(g), float(var)


def pool(gv: dict[str, tuple[float, float]]) -> dict:
    keys = [k for k, v in gv.items() if np.isfinite(v[0])]
    g = np.array([[gv[k][0] for k in keys]]).T
    v = np.array([[gv[k][1] for k in keys]]).T
    r = dersimonian_laird(g, v)
    return {k: float(np.ravel(val)[0]) for k, val in r.items()}


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    sig_full = pd.read_csv(DEC / "signature_matrix.tsv", sep="\t", index_col=0)
    sig_full = sig_full[~sig_full.index.isin(GRANULOCYTE)]

    # The signature used for proportions has every effector gene removed, so
    # the proportion estimate cannot be driven by the thing we then test.
    sig = sig_full[~sig_full.index.isin(EFFECTOR)]
    dropped = len(sig_full) - len(sig)
    print(f"signature: {len(sig_full)} markers, {dropped} effector genes removed "
          f"-> {len(sig)} used for proportion estimation")

    grp = pd.Series(COARSE)
    props_all, rows = {}, []

    print("\n" + "=" * 96)
    print("DECONVOLUTION OF THE BULK COHORTS")
    print("=" * 96)

    for gse in PRIMARY:
        expr, meta = load_bulk(gse)
        lin = to_linear(expr, gse)
        shared = sig.index.intersection(lin.index)
        est = pd.DataFrame(
            {s: deconvolve(lin[s], sig.loc[shared], "svr") for s in lin.columns}
        )
        coarse = est.groupby(grp.reindex(est.index)).sum()
        props_all[gse] = coarse

        n_as = int(meta["group"].eq("AS").sum())
        print(f"\n  {gse}  ({len(shared)} of {len(sig)} signature genes present, "
              f"AS {n_as} / HC {len(meta) - n_as})")
        print(f"    {'lineage':10s} {'AS mean':>8s} {'HC mean':>8s} {'g':>7s} {'p':>8s}")
        for lin_name in LINEAGES:
            if lin_name not in coarse.index:
                continue
            v = coarse.loc[lin_name]
            g, var = effect_of(v, meta)
            a = v[(meta["group"] == "AS").to_numpy()]
            c = v[(meta["group"] == "control").to_numpy()]
            p = stats.mannwhitneyu(a, c).pvalue
            flag = "" if lin_name in RELIABLE else "  (unreliable)"
            print(f"    {lin_name:10s} {a.mean():8.3f} {c.mean():8.3f} "
                  f"{g:+7.2f} {p:8.3f}{flag}")
            rows.append(dict(gse=gse, lineage=lin_name, g=g, var_g=var,
                             as_mean=a.mean(), hc_mean=c.mean(), p=p))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "deconv_per_study.tsv", sep="\t", index=False)

    # ----------------------------------------------------------------- pooled
    print("\n" + "=" * 96)
    print("POOLED -- is the cell-type PROPORTION different in AS?")
    print("=" * 96)
    print(f"  {'lineage':10s} {'k':>2s} {'pooled g':>9s} {'95% CI':>18s} "
          f"{'p':>8s} {'I2':>6s}   reliability")
    prop_pooled = {}
    for lin_name in LINEAGES:
        sub = df[df["lineage"] == lin_name]
        if sub.empty:
            continue
        r = pool({row.gse: (row.g, row.var_g) for row in sub.itertuples()})
        prop_pooled[lin_name] = r
        rel = "validated r>0.8" if lin_name in RELIABLE else "NOT validated"
        print(f"  {lin_name:10s} {int(r['k']):2d} {r['pooled_g']:+9.3f} "
              f"[{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}]  {r['p']:8.4f} "
              f"{r['i2']:6.1f}   {rel}")

    # ------------------------------------------------- composition vs regulation
    print("\n" + "=" * 96)
    print("COMPOSITION OR REGULATION?")
    print("=" * 96)
    print("  Effector cargo (granzymes/perforin/granulysin), before and after")
    print("  adjusting for the NK proportion estimated from DISJOINT identity genes.\n")
    print(f"  {'study':10s} {'eff g raw':>10s} {'eff g | NKprop':>15s} "
          f"{'NK prop g':>10s} {'r(eff,prop)':>12s}")

    raw_gv, adj_gv = {}, {}
    for gse in PRIMARY:
        expr, meta = load_bulk(gse)
        z = zscore_rows(expr)
        eff = [g for g in EFFECTOR if g in z.index]
        if len(eff) < 4:
            print(f"  {gse:10s} only {len(eff)} effector genes measured; skipped")
            continue
        eff_score = z.loc[eff].mean(axis=0)
        nk_prop = props_all[gse].loc["NK"]

        g_raw, v_raw = effect_of(eff_score, meta)
        gp, vp = effect_of(nk_prop, meta)

        # Adjust: regress the effector score on the NK proportion, then test the
        # residual. What is left is cargo per NK cell, not total cargo.
        m2 = meta.copy()
        m2["_nkprop"] = (nk_prop - nk_prop.mean()) / nk_prop.std()
        X, names, idx = build_design(m2)
        extra = m2.loc[idx, "_nkprop"].to_numpy(float)[:, None]
        Xa = np.hstack([X, extra])
        y = eff_score.loc[idx].to_numpy(float)
        dof = len(y) - Xa.shape[1]
        beta, *_ = np.linalg.lstsq(Xa, y, rcond=None)
        s = np.sqrt(((y - Xa @ beta) ** 2).sum() / dof)
        j = names.index("_AS")
        g_adj = (1.0 - 3.0 / (4.0 * dof - 1.0)) * beta[j] / s
        v_adj = np.linalg.inv(Xa.T @ Xa)[j, j] + g_adj**2 / (2.0 * dof)

        r_ep = float(stats.pearsonr(eff_score.loc[idx], nk_prop.loc[idx])[0])
        print(f"  {gse:10s} {g_raw:+10.2f} {g_adj:+15.2f} {gp:+10.2f} {r_ep:+12.2f}")
        raw_gv[gse] = (g_raw, v_raw)
        adj_gv[gse] = (g_adj, float(v_adj))

    pr, pa = pool(raw_gv), pool(adj_gv)
    print(f"\n  pooled effector score, unadjusted      g = {pr['pooled_g']:+.3f} "
          f"[{pr['ci_lo']:+.2f}, {pr['ci_hi']:+.2f}]  p = {pr['p']:.4g}  I2 = {pr['i2']:.0f}%")
    print(f"  pooled effector score, NK-prop-adjusted g = {pa['pooled_g']:+.3f} "
          f"[{pa['ci_lo']:+.2f}, {pa['ci_hi']:+.2f}]  p = {pa['p']:.4g}  I2 = {pa['i2']:.0f}%")
    nk = prop_pooled.get("NK", {})
    if nk:
        print(f"  pooled NK PROPORTION                    g = {nk['pooled_g']:+.3f} "
              f"[{nk['ci_lo']:+.2f}, {nk['ci_hi']:+.2f}]  p = {nk['p']:.4g}")

    retained = (pa["pooled_g"] / pr["pooled_g"] * 100) if pr["pooled_g"] else np.nan
    print(f"\n  {retained:.0f}% of the effector effect survives adjustment for NK proportion.")
    print("  Near 0% would mean composition explains it. Near 100% would mean the")
    print("  cargo is reduced per cell, independent of how many cells there are.")

    pd.DataFrame(prop_pooled).T.to_csv(OUT / "deconv_pooled.tsv", sep="\t")
    for gse, c in props_all.items():
        c.to_csv(DEC / f"props_{gse}.tsv", sep="\t")
    print(f"\nwrote {OUT / 'deconv_pooled.tsv'} and per-study proportions to {DEC}")


if __name__ == "__main__":
    main()
