"""Can this reference actually recover cell-type proportions? Validate first.

Deconvolution is reliable for abundant, distinct cell types and poor for rare
or similar ones, and nothing about the output tells you which case you are in.
So before any bulk cohort is touched, deconvolve held-out pseudobulk built from
patients the signature never saw, where the true proportions ARE known.

A cell type only earns the right to be interpreted in the bulk data if it is
recovered here.

Usage:  python scripts/deconv_validate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy import stats
from sklearn.svm import NuSVR

ROOT = Path(__file__).resolve().parents[1]
DEC = ROOT / "data" / "deconv"

# Granulocyte / neutrophil genes. The reference is PBMC, so it has no
# neutrophil component at all -- but GSE73754 and GSE25101 are whole blood,
# where neutrophils are the majority leukocyte. Their transcripts would be
# forced onto whichever PBMC type looks least unlike them. Excluding these
# genes from the marker set means the fit is driven by mononuclear signal, and
# proportions are then reported WITHIN the PBMC compartment.
GRANULOCYTE = {
    "FCGR3B", "CSF3R", "CXCR2", "CXCR1", "S100A12", "MMP8", "MMP9", "LCN2",
    "DEFA1", "DEFA1B", "DEFA3", "DEFA4", "ELANE", "MPO", "PRTN3", "CEACAM8",
    "CAMP", "LTF", "CTSG", "BPI", "ARG1", "OLFM4", "PGLYRP1", "TCN1",
    "ALPL", "SIGLEC5", "TREM1", "VNN2", "IL1R2", "CSF2RB", "SELL",
    "HBB", "HBA1", "HBA2", "HBD", "ALAS2", "AHSP", "SLC4A1",  # globin: PAXgene
}

COARSE = {
    "CD16+ NK cells": "NK", "CD16- NK cells": "NK", "NK cells": "NK",
    "Tcm/Naive helper T cells": "CD4 T", "Tem/Effector helper T cells": "CD4 T",
    "Regulatory T cells": "CD4 T",
    "Tcm/Naive cytotoxic T cells": "CD8 T", "Tem/Temra cytotoxic T cells": "CD8 T",
    "Tem/Trm cytotoxic T cells": "CD8 T", "MAIT cells": "CD8 T",
    "Naive B cells": "B", "Memory B cells": "B", "Age-associated B cells": "B",
    "Plasma cells": "B",
    "Classical monocytes": "monocyte", "Non-classical monocytes": "monocyte",
    "DC1": "DC", "DC2": "DC", "pDC": "DC",
    "Megakaryocytes/platelets": "platelet", "Megakaryocyte precursor": "platelet",
    "ILC3": "other", "Late erythroid": "other",
}


def load_signature() -> pd.DataFrame:
    sig = pd.read_csv(DEC / "signature_matrix.tsv", sep="\t", index_col=0)
    before = len(sig)
    sig = sig[~sig.index.isin(GRANULOCYTE)]
    print(f"signature: {before} markers -> {len(sig)} after dropping "
          f"granulocyte/globin genes, {sig.shape[1]} cell types")
    return sig


def deconvolve(mixture: pd.Series, sig: pd.DataFrame, method: str = "svr") -> pd.Series:
    """Estimate cell-type fractions for one sample.

    Both matrices are standardised before fitting (CIBERSORT's trick): it makes
    the fit invariant to the fact that a microarray intensity and a single-cell
    CPM have completely different scales and dynamic ranges.
    """
    shared = sig.index.intersection(mixture.index)
    S = sig.loc[shared].to_numpy(float)
    m = mixture.loc[shared].to_numpy(float)

    keep = np.isfinite(m) & np.isfinite(S).all(axis=1)
    S, m = S[keep], m[keep]

    if method == "nnls":
        coef, _ = nnls(S, m)
    else:
        Sz = (S - S.mean()) / S.std()
        mz = (m - m.mean()) / m.std()
        best, best_rmse = None, np.inf
        for nu in (0.25, 0.5, 0.75):
            try:
                mdl = NuSVR(nu=nu, C=1.0, kernel="linear").fit(Sz, mz)
            except Exception:
                continue
            w = mdl.coef_.ravel()
            rmse = float(np.sqrt(np.mean((Sz @ w - mz) ** 2)))
            if rmse < best_rmse:
                best, best_rmse = w, rmse
        coef = np.clip(best, 0, None)

    total = coef.sum()
    return pd.Series(coef / total if total > 0 else coef, index=sig.columns)


def score(est: pd.DataFrame, true: pd.DataFrame, label: str) -> pd.DataFrame:
    """Per cell type: correlation with truth, error, and bias."""
    types = [t for t in true.index if t in est.index]
    rows = []
    for t in types:
        e, v = est.loc[t], true.loc[t]
        ok = np.isfinite(e) & np.isfinite(v)
        if ok.sum() < 3 or v[ok].std() == 0:
            r = np.nan
        else:
            r = float(stats.pearsonr(e[ok], v[ok])[0])
        rows.append(dict(
            cell_type=t,
            true_mean=float(v.mean()),
            est_mean=float(e.mean()),
            bias=float(e.mean() - v.mean()),
            r=r,
            rmse=float(np.sqrt(np.mean((e[ok] - v[ok]) ** 2))),
        ))
    df = pd.DataFrame(rows).sort_values("true_mean", ascending=False)
    print(f"\n--- {label} ---")
    with pd.option_context("display.width", 200):
        print(df.round(4).to_string(index=False))
    return df


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    sig = load_signature()
    pb = pd.read_csv(DEC / "heldout_pseudobulk.tsv.gz", sep="\t", index_col=0)
    true = pd.read_csv(DEC / "heldout_true_props.tsv", sep="\t", index_col=0)
    print(f"held-out pseudobulk: {pb.shape[1]} patients x {pb.shape[0]:,} genes")

    # pseudobulk is summed raw counts -> CPM, matching the signature's scale
    pb = pb / pb.sum(axis=0) * 1e6

    for method in ("svr", "nnls"):
        est = pd.DataFrame(
            {s: deconvolve(pb[s], sig, method) for s in pb.columns}
        )
        fine = score(est, true, f"FINE, 23 types, method={method}")

        # Collapse to coarse lineages and re-score: summing related types
        # cancels the confusion between them, which is where most of the error
        # in fine-grained deconvolution lives.
        grp = pd.Series(COARSE)
        est_c = est.groupby(grp.reindex(est.index)).sum()
        true_c = true.groupby(grp.reindex(true.index)).sum()
        coarse = score(est_c, true_c, f"COARSE, 7 lineages, method={method}")

        fine.to_csv(DEC / f"validation_fine_{method}.tsv", sep="\t", index=False)
        coarse.to_csv(DEC / f"validation_coarse_{method}.tsv", sep="\t", index=False)

        good = coarse[(coarse["r"] > 0.7) & (coarse["true_mean"] > 0.01)]
        print(f"\n  method={method}: lineages recovered with r>0.7 and >1% abundance: "
              f"{list(good['cell_type'])}")
        nk = coarse[coarse["cell_type"] == "NK"]
        if len(nk):
            print(f"  NK specifically: r={nk['r'].iloc[0]:.3f}, "
                  f"true mean {nk['true_mean'].iloc[0]:.3f}, "
                  f"estimated {nk['est_mean'].iloc[0]:.3f}, "
                  f"bias {nk['bias'].iloc[0]:+.3f}")


if __name__ == "__main__":
    main()
