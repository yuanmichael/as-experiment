"""AS vs healthy control: per-study effect sizes, then random-effects meta-analysis.

Design (the point of the whole exercise): every AS-vs-control comparison happens
*inside* one study, so the batch effect cancels. Only standardised effect sizes
cross study boundaries -- never expression values.

Per study, per gene, fit  log2(expr) ~ AS + covariates  by least squares, then
convert the AS coefficient to a standardised mean difference:

    g      = J * beta_AS / sqrt(MSE)              (SD units, covariate-adjusted)
    var(g) = (X'X)^-1[AS,AS] + g^2 / (2 * dof)    (design term + SD-estimate term)
    J      = 1 - 3 / (4*dof - 1)                  (Hedges' small-sample correction)

All genes in a study share one design matrix, so one lstsq solves every gene at
once -- no per-gene regression loop.

Then pool across studies with DerSimonian-Laird random effects, and report
pooled g, CI, p, Benjamini-Hochberg q, and between-study heterogeneity I^2.

Usage:  python scripts/bulk_meta.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "data" / "bulk" / "prepared"
OUT = ROOT / "data" / "bulk" / "meta"

# Primary: circulating blood, unstimulated. One biological question.
PRIMARY = ["GSE73754", "GSE25101", "GSE18781", "GSE299639", "GSE181364"]

# Secondary: internally valid AS-vs-HC contrasts, but of a different quantity
# (sorted or perturbed cells). Run separately; never pooled with the primary.
SECONDARY = ["GSE221786", "GSE11886"]

# Covariates to use where the study reports them. Disease-defining variables
# (HLA-B27, disease activity) are deliberately excluded -- adjusting for them
# would adjust away part of the contrast.
COVARS = ["sex", "age", "set"]

# Genes the previous single-cell project implicated, tested as a named panel
# rather than discovered post hoc.
PANEL = {
    "AIF1 (the paper's marker)": ["AIF1"],
    "monocyte identity": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN", "FCN1", "CSF1R"],
    "non-classical monocyte": ["FCGR3A", "MS4A7", "CDKN1C"],
    "type-1 IFN response": [
        "IFI6", "IFI44", "IFI44L", "IFI27", "ISG15", "MX1", "MX2", "OAS1",
        "OAS2", "OAS3", "IFIT1", "IFIT3", "RSAD2", "STAT1", "STAT2", "IRF7",
        "XAF1", "EPSTI1", "SIGLEC1", "USP18", "LY6E", "SPATS2L", "HERC5",
    ],
    "IL-17 / IL-23 axis": ["IL17A", "IL17F", "IL23A", "IL23R", "RORC", "IL6R", "TNF"],
}


# --------------------------------------------------------------------------
# per-study effect sizes
# --------------------------------------------------------------------------

def build_design(meta: pd.DataFrame) -> tuple[np.ndarray, list[str], pd.Index]:
    """Intercept + AS indicator + usable covariates. Returns (X, names, samples)."""
    df = pd.DataFrame(index=meta.index)
    df["_AS"] = (meta["group"] == "AS").astype(float)

    for cv in COVARS:
        if cv not in meta.columns:
            continue
        col = meta[cv]
        if cv == "age":
            num = pd.to_numeric(col, errors="coerce")
            if num.notna().sum() < len(num) or num.nunique() < 3:
                continue
            df["age"] = (num - num.mean()) / num.std()
        else:
            lvl = col.astype(str).str.strip().str.lower().replace({"nan": np.nan})
            if lvl.isna().any() or lvl.nunique() < 2:
                continue
            # one level is absorbed into the intercept
            for v in sorted(lvl.unique())[1:]:
                df[f"{cv}_{v}"] = (lvl == v).astype(float)

    # A covariate that is collinear with the AS indicator cannot be adjusted for
    # -- it would be absorbing the contrast itself. Drop it and say so.
    for c in [c for c in df.columns if c != "_AS"]:
        if abs(np.corrcoef(df["_AS"], df[c])[0, 1]) > 0.95:
            print(f"      dropping covariate {c!r}: collinear with disease status")
            df = df.drop(columns=c)

    X = np.column_stack([np.ones(len(df)), df.to_numpy(float)])
    names = ["intercept"] + list(df.columns)
    return X, names, df.index


def study_effects(gse: str) -> pd.DataFrame:
    expr = pd.read_csv(PREP / f"{gse}_expr.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(PREP / f"{gse}_samples.tsv", sep="\t", index_col=0)
    expr = expr[[c for c in expr.columns if c in meta.index]]
    meta = meta.loc[expr.columns]

    X, names, _ = build_design(meta)
    n, p = X.shape
    dof = n - p
    j_as = names.index("_AS")

    Y = expr.to_numpy(float).T  # samples x genes

    # Drop genes with no variance -- the model is undefined and they are not
    # measuring anything.
    keep = Y.std(axis=0) > 0
    Y, genes = Y[:, keep], expr.index[keep]

    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    mse = (resid**2).sum(axis=0) / dof
    s = np.sqrt(mse)

    # g = beta / s, so a gene whose residual SD is near zero reports an enormous
    # effect from a trivial difference. At these sample sizes that is an
    # artefact of the fit, not a finding. Drop the degenerate tail.
    floor = 0.05 * np.median(s)
    ok = s > floor
    if not ok.all():
        print(f"      dropping {int((~ok).sum())} genes with residual SD < {floor:.4g}")
        Y, genes, beta, s = Y[:, ok], genes[ok], beta[:, ok], s[ok]

    xtx_inv = np.linalg.inv(X.T @ X)
    design_var = xtx_inv[j_as, j_as]

    g = beta[j_as] / s
    j_corr = 1.0 - 3.0 / (4.0 * dof - 1.0)
    g = j_corr * g
    var_g = design_var + g**2 / (2.0 * dof)

    n_as = int(meta["group"].eq("AS").sum())
    n_ctl = int(meta["group"].eq("control").sum())
    print(
        f"      n={n} (AS {n_as}, HC {n_ctl})  dof={dof}  "
        f"covars={[x for x in names if x not in ('intercept', '_AS')]}  "
        f"genes={len(genes)}"
    )

    return pd.DataFrame(
        {
            "gse": gse,
            "gene": genes,
            "g": g,
            "var_g": var_g,
            "se_g": np.sqrt(var_g),
            "log2fc": beta[j_as],
            "n_as": n_as,
            "n_ctl": n_ctl,
        }
    ).set_index("gene")


# --------------------------------------------------------------------------
# DerSimonian-Laird random effects
# --------------------------------------------------------------------------

def dersimonian_laird(g: np.ndarray, v: np.ndarray) -> dict:
    """g, v: (k_studies x n_genes), NaN where a gene is absent from a study."""
    mask = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g = np.where(mask, g, 0.0)
    w = np.where(mask, 1.0 / np.where(mask, v, 1.0), 0.0)

    k = mask.sum(axis=0)
    sw = w.sum(axis=0)
    fixed = (w * g).sum(axis=0) / np.where(sw > 0, sw, np.nan)

    q = (w * (g - fixed) ** 2).sum(axis=0)
    c = sw - (w**2).sum(axis=0) / np.where(sw > 0, sw, np.nan)
    tau2 = np.maximum(0.0, (q - (k - 1)) / np.where(c > 0, c, np.nan))

    wr = np.where(mask, 1.0 / (np.where(mask, v, 1.0) + tau2), 0.0)
    swr = wr.sum(axis=0)
    pooled = (wr * g).sum(axis=0) / np.where(swr > 0, swr, np.nan)
    se = np.sqrt(1.0 / np.where(swr > 0, swr, np.nan))

    z = pooled / se
    pval = 2.0 * stats.norm.sf(np.abs(z))
    i2 = np.where(q > 0, np.maximum(0.0, (q - (k - 1)) / np.where(q > 0, q, np.nan)), 0.0)

    # direction agreement across the studies that measured the gene
    sign = np.where(mask, np.sign(g), np.nan)
    with np.errstate(invalid="ignore"):
        same_dir = np.nanmax(
            np.stack([(sign > 0).sum(axis=0), (sign < 0).sum(axis=0)]), axis=0
        )

    return dict(
        k=k, pooled_g=pooled, se=se, ci_lo=pooled - 1.96 * se, ci_hi=pooled + 1.96 * se,
        z=z, p=pval, tau2=tau2, q=q, i2=i2 * 100.0, fixed_g=fixed,
        n_same_direction=same_dir,
    )


def meta_analyse(studies: list[str], label: str, min_k: int) -> pd.DataFrame:
    print(f"\n{label}")
    per = {}
    for gse in studies:
        print(f"    {gse}")
        per[gse] = study_effects(gse)

    g_wide = pd.DataFrame({k: v["g"] for k, v in per.items()})
    v_wide = pd.DataFrame({k: v["var_g"] for k, v in per.items()})
    g_wide, v_wide = g_wide.align(v_wide, join="outer")

    present = g_wide.notna().sum(axis=1)
    keep = present >= min_k
    print(
        f"    genes measured in >={min_k} of {len(studies)} studies: "
        f"{int(keep.sum())} of {len(keep)}"
    )
    g_wide, v_wide = g_wide[keep], v_wide[keep]

    res = dersimonian_laird(g_wide.to_numpy().T, v_wide.to_numpy().T)
    out = pd.DataFrame(res, index=g_wide.index)
    out["q"] = multipletests(out["p"].fillna(1.0), method="fdr_bh")[1]
    for gse in studies:
        out[f"g_{gse}"] = g_wide.get(gse)
    return out.sort_values("p")


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def detectable_effect(n1: int, n2: int, alpha: float = 0.05, power: float = 0.80) -> float:
    """Two-sided two-sample t-test: smallest d detectable at the given power."""
    z_a = stats.norm.ppf(1 - alpha / 2)
    z_b = stats.norm.ppf(power)
    return (z_a + z_b) * np.sqrt(1.0 / n1 + 1.0 / n2)


def report(out: pd.DataFrame, label: str, top: int = 25) -> None:
    sig = out[out["q"] < 0.05]
    print(f"\n  {label}: {len(sig)} genes at q<0.05, "
          f"{int((out['q'] < 0.10).sum())} at q<0.10")
    cols = ["k", "pooled_g", "ci_lo", "ci_hi", "p", "q", "i2", "n_same_direction"]
    show = (sig if len(sig) else out).head(top)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(show[cols].round(4).to_string())


def report_panel(out: pd.DataFrame, studies: list[str]) -> None:
    print("\n" + "=" * 100)
    print("NAMED PANEL -- hypotheses carried over from the single-cell project")
    print("=" * 100)
    cols = ["k", "pooled_g", "ci_lo", "ci_hi", "p", "q", "i2"] + [
        f"g_{s}" for s in studies
    ]
    for name, genes in PANEL.items():
        have = [g for g in genes if g in out.index]
        missing = [g for g in genes if g not in out.index]
        print(f"\n{name}  ({len(have)}/{len(genes)} measured)")
        if missing:
            print(f"  not measured: {missing}")
        if not have:
            continue
        sub = out.loc[have, cols].sort_values("p")
        with pd.option_context("display.width", 220, "display.max_columns", 30):
            print(sub.round(3).to_string())
        # Is the panel as a whole shifted? Sign test on pooled direction.
        pos = int((sub["pooled_g"] > 0).sum())
        n = len(sub)
        print(
            f"  direction: {pos}/{n} up in AS   "
            f"mean pooled g = {sub['pooled_g'].mean():+.3f}   "
            f"sign-test p = {stats.binomtest(pos, n, 0.5).pvalue:.3g}"
        )


LINEAGE = {
    "T cell": ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK", "TRAC", "IL7R", "CD28"],
    "NK": ["KLRD1", "KLRF1", "NKG7", "GNLY", "PRF1", "EOMES", "KLRB1", "NCR3", "KLRC1"],
    "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "BLNK", "TCL1A", "BANK1"],
    "monocyte": ["CD14", "LYZ", "VCAN", "FCN1", "CSF1R", "CD68", "MNDA", "CEBPB"],
    "neutrophil": ["ELANE", "DEFA3", "DEFA4", "MPO", "CSF3R", "FCGR3B", "CEACAM8", "LCN2"],
    "platelet": ["PF4", "PPBP", "ITGA2B", "GP9", "TUBB1"],
    "erythroid": ["HBB", "HBD", "ALAS2", "AHSP", "SLC4A1"],
}


def sensitivity() -> None:
    """Does the conclusion depend on any single study, or on the two small
    RNA-seq cohorts that carry a residual global shift?"""
    print("\n" + "=" * 100)
    print("SENSITIVITY -- leave-one-study-out, and arrays only")
    print("=" * 100)

    subsets = {"all 5 studies": PRIMARY, "arrays only (3)": PRIMARY[:3]}
    for drop in PRIMARY:
        subsets[f"drop {drop}"] = [s for s in PRIMARY if s != drop]

    rows = []
    for label, studies in subsets.items():
        min_k = min(3, len(studies))
        out = meta_analyse(studies, f"  [{label}]", min_k=min_k)
        hits = out[out["q"] < 0.05]
        row = {
            "subset": label,
            "k_studies": len(studies),
            "hits_q05": len(hits),
            "frac_down": round(float((hits["pooled_g"] < 0).mean()), 3) if len(hits) else np.nan,
        }
        for lin, genes in LINEAGE.items():
            have = [g for g in genes if g in out.index]
            row[lin] = round(float(out.loc[have, "pooled_g"].mean()), 3) if have else np.nan
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "sensitivity.tsv", sep="\t", index=False)
    print("\n  Mean pooled g per lineage marker set, by subset:")
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(df.to_string(index=False))
    print("\n  A finding that survives every leave-one-out column is not one")
    print("  study's artefact. A lineage column that moves while others stay")
    print("  flat is composition, not a global shift.")


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("POWER -- the smallest standardised effect each study could detect")
    print("=" * 100)
    for gse in PRIMARY + SECONDARY:
        meta = pd.read_csv(PREP / f"{gse}_samples.tsv", sep="\t", index_col=0)
        n1 = int(meta["group"].eq("AS").sum())
        n2 = int(meta["group"].eq("control").sum())
        print(f"  {gse:10s} AS={n1:3d} HC={n2:3d}   d(80% power) = {detectable_effect(n1, n2):.2f}")
    tot_as = sum(
        int(pd.read_csv(PREP / f"{g}_samples.tsv", sep="\t")["group"].eq("AS").sum())
        for g in PRIMARY
    )
    tot_ctl = sum(
        int(pd.read_csv(PREP / f"{g}_samples.tsv", sep="\t")["group"].eq("control").sum())
        for g in PRIMARY
    )
    print(
        f"  {'POOLED':10s} AS={tot_as:3d} HC={tot_ctl:3d}   "
        f"d(80% power) = {detectable_effect(tot_as, tot_ctl):.2f}  "
        f"(upper bound: random effects with tau^2>0 has less power than this)"
    )

    print("\n" + "=" * 100)
    primary = meta_analyse(PRIMARY, "PRIMARY -- circulating blood, unstimulated", min_k=3)
    primary.to_csv(OUT / "meta_primary.tsv", sep="\t")
    report(primary, "PRIMARY")
    report_panel(primary, PRIMARY)

    print("\n" + "=" * 100)
    secondary = meta_analyse(
        SECONDARY, "SECONDARY -- sorted / perturbed cells (different quantity)", min_k=2
    )
    secondary.to_csv(OUT / "meta_secondary.tsv", sep="\t")
    report(secondary, "SECONDARY")

    # Does the primary signal reappear in an independent measurement system?
    print("\n" + "=" * 100)
    print("CONVERGENCE -- primary hits vs the secondary (sorted/perturbed) studies")
    print("=" * 100)
    hits = primary[primary["q"] < 0.05]
    shared = hits.index.intersection(secondary.index)
    if len(shared) == 0:
        print("  no primary hits to check")
    else:
        agree = np.sign(hits.loc[shared, "pooled_g"]) == np.sign(
            secondary.loc[shared, "pooled_g"]
        )
        print(
            f"  {len(shared)} of {len(hits)} primary hits also measured in secondary; "
            f"{int(agree.sum())} ({100 * agree.mean():.0f}%) agree in direction "
            f"(chance = 50%)"
        )
        print(f"  sign-test p = {stats.binomtest(int(agree.sum()), len(shared), 0.5).pvalue:.3g}")

    sensitivity()

    print(f"\nwrote {OUT / 'meta_primary.tsv'} and {OUT / 'meta_secondary.tsv'}")


if __name__ == "__main__":
    main()
