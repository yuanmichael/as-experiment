"""How well does the NK signature actually separate AS from healthy?

A pooled effect size says how far apart two group means are. It does not say
whether you could classify a person. That question needs an AUC, computed
per study -- the probability that a random AS patient scores below a random
healthy control.

Also diagnoses GSE25101, whose EOMES effect is the weakest of the five.

Usage:  python scripts/bulk_score.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "data" / "bulk" / "prepared"

PRIMARY = ["GSE73754", "GSE25101", "GSE18781", "GSE299639", "GSE181364"]

SIGNATURES = {
    "NK / cytotoxic": ["KLRD1", "KLRF1", "NKG7", "GNLY", "PRF1", "EOMES",
                       "KLRB1", "NCR3", "KLRC1", "GZMA", "GZMK", "FGFBP5",
                       "S1PR5", "FGFBP2", "IL2RB", "CD247", "PLEKHF1"],
    "T cell": ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK", "TRAC", "IL7R", "CD28"],
    "type-1 IFN": ["IFI6", "IFI44", "IFI44L", "IFI27", "ISG15", "MX1", "MX2",
                   "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT3", "RSAD2", "IRF7",
                   "LY6E", "HERC5", "EPSTI1"],
}


def load(gse: str):
    expr = pd.read_csv(PREP / f"{gse}_expr.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(PREP / f"{gse}_samples.tsv", sep="\t", index_col=0)
    expr = expr[[c for c in expr.columns if c in meta.index]]
    return expr, meta.loc[expr.columns]


def score(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, int]:
    """Mean z-score across the signature's genes, one value per sample.

    Each gene is z-scored within the study first, so a highly-expressed gene
    does not dominate the average purely because of its scale.
    """
    have = [g for g in genes if g in expr.index]
    sub = expr.loc[have]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0), len(have)


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    print("=" * 96)
    print("SIGNATURE AUC -- could you classify a person with this?")
    print("=" * 96)
    print("  AUC = P(a random AS patient scores more extreme than a random control).")
    print("  0.50 = coin flip. Computed inside each study; nothing is trained.\n")

    data = {g: load(g) for g in PRIMARY}
    rows = []
    for sig, genes in SIGNATURES.items():
        print(f"  {sig}")
        for gse, (expr, meta) in data.items():
            s, n_have = score(expr, genes)
            y = (meta["group"] == "AS").to_numpy().astype(int)
            auc = roc_auc_score(y, -s)  # signature is LOW in AS, so negate
            a, c = s[y == 1], s[y == 0]
            sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(c) - 1) * c.var(ddof=1))
                         / (len(a) + len(c) - 2))
            d = (a.mean() - c.mean()) / sp
            print(f"    {gse:10s} genes {n_have:2d}/{len(genes)}  "
                  f"AUC {auc:.3f}   score d = {d:+.2f}   "
                  f"(AS {int(y.sum())}, HC {int((1 - y).sum())})")
            rows.append(dict(signature=sig, gse=gse, auc=auc, d=d,
                             n_genes=n_have, n_as=int(y.sum()),
                             n_hc=int((1 - y).sum())))
        print()

    df = pd.DataFrame(rows)
    print("  sample-size-weighted mean AUC per signature:")
    for sig, g in df.groupby("signature", sort=False):
        w = g["n_as"] + g["n_hc"]
        print(f"    {sig:16s} {np.average(g['auc'], weights=w):.3f}   "
              f"(range {g['auc'].min():.3f}-{g['auc'].max():.3f})")
    df.to_csv(ROOT / "data" / "bulk" / "meta" / "signature_auc.tsv", sep="\t", index=False)

    # ----------------------------------------------------------------------
    print("\n" + "=" * 96)
    print("GSE25101 -- why is its NK effect the weakest of the five?")
    print("=" * 96)

    nk = SIGNATURES["NK / cytotoxic"]
    print("\n  per-study g for each NK gene (from the meta table):")
    meta_tbl = pd.read_csv(
        ROOT / "data" / "bulk" / "meta" / "meta_primary.tsv", sep="\t", index_col=0
    )
    cols = [f"g_{s}" for s in PRIMARY]
    have = [g for g in nk if g in meta_tbl.index]
    sub = meta_tbl.loc[have, cols]
    with pd.option_context("display.width", 200):
        print(sub.round(2).to_string())
    print("\n  mean g per study over the NK panel:")
    print(sub.mean().round(3).to_string())

    print("\n  Hypothesis 1: GSE25101 attenuates everything (a scale problem).")
    print("    If so its whole distribution of |g| would be compressed.")
    for gse in PRIMARY:
        col = meta_tbl[f"g_{gse}"].dropna()
        print(f"    {gse:10s} sd(g) over all genes = {col.std():.3f}   "
              f"IQR = {col.quantile(.75) - col.quantile(.25):.3f}")

    print("\n  Hypothesis 2: the NK genes sit near background on that platform,")
    print("    so the measurement is noise-dominated.")
    for gse in PRIMARY:
        expr, _ = data[gse]
        allm = expr.mean(axis=1)
        h = [g for g in nk if g in expr.index]
        pct = [float((allm < allm[g]).mean()) for g in h]
        print(f"    {gse:10s} NK genes at mean expression percentile "
              f"{100 * np.mean(pct):.0f} (median {100 * np.median(pct):.0f}), "
              f"{len(h)}/{len(nk)} measured")

    print("\n  Hypothesis 3: the tissue differs. GSE25101's own records contradict")
    print("    themselves (whole blood vs PBMC). Neutrophil content is the tell:")
    neut = ["FCGR3B", "CSF3R", "CXCR2", "S100A12", "MMP9", "LCN2", "DEFA3", "CEACAM8"]
    for gse in PRIMARY:
        expr, _ = data[gse]
        allm = expr.mean(axis=1)
        h = [g for g in neut if g in expr.index]
        if not h:
            print(f"    {gse:10s} no neutrophil markers measured")
            continue
        pct = [float((allm < allm[g]).mean()) for g in h]
        print(f"    {gse:10s} neutrophil markers at expression percentile "
              f"{100 * np.mean(pct):.0f}  ({len(h)}/{len(neut)} measured)")
    print("\n    High percentile = neutrophil transcripts abundant = whole blood.")
    print("    Low = neutrophils absent = PBMC (they are removed by Ficoll).")


if __name__ == "__main__":
    main()
