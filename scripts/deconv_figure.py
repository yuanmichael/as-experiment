"""Figure for the deconvolution: what is trustworthy, and what it says.

Same colour language as figures/bulk_meta.png -- blue = lower in AS, red =
higher in AS, grey = neutral/diagnostic.

Usage:  python scripts/deconv_figure.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEC = ROOT / "data" / "deconv"
META = ROOT / "data" / "bulk" / "meta"
FIG = ROOT / "figures"

DOWN = "#2a78d6"
UP = "#e34948"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"

mpl.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": BASE,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
})

ORDER = ["NK", "monocyte", "B", "CD4 T", "CD8 T", "DC", "platelet"]
RELIABLE = {"NK", "monocyte", "B", "CD4 T"}

# from deconv_bulk.py / the sink check
EFFECTOR_RAW = (-0.978, -1.33, -0.63)
EFFECTOR_ADJ = (-0.804, -1.17, -0.44)
NK_PROP = (-0.540, -0.91, -0.17)


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    val = pd.read_csv(DEC / "validation_coarse_svr.tsv", sep="\t").set_index("cell_type")
    pooled = pd.read_csv(META / "deconv_pooled.tsv", sep="\t", index_col=0)

    fig = plt.figure(figsize=(14.0, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.05], wspace=0.42,
                          left=0.075, right=0.975, top=0.76, bottom=0.14)

    # ------------------------------------------------------------- panel A
    ax = fig.add_subplot(gs[0, 0])
    lin = [l for l in ORDER if l in val.index]
    r = val.loc[lin, "r"].to_numpy()
    y = np.arange(len(lin))[::-1]
    ax.barh(y, r, height=0.55, color=[INK2 if v >= 0.8 else BASE for v in r],
            edgecolor=SURFACE, linewidth=2)
    ax.axvline(0.8, color=UP, lw=1.3, ls=(0, (4, 3)))
    ax.text(0.815, y[-1] - 0.30, "r = 0.8 cutoff", fontsize=7.5, color=UP,
            va="center", ha="left")
    for yi, v, l in zip(y, r, lin):
        ax.text(max(v, 0) + 0.02, yi, f"{v:.2f}", va="center", fontsize=8.5,
                color=INK2)
    ax.set_yticks(y, [f"{l}  ({val.loc[l, 'true_mean']:.0%})" for l in lin], fontsize=8.5)
    ax.set_xlim(0, 1.22)
    ax.set_xlabel("correlation with the true proportion")
    ax.set_title("A   Validation on held-out pseudobulk",
                 loc="left", fontsize=11, fontweight="bold", pad=34)
    ax.text(0, 1.012, "patients the signature never saw; true proportions known\n"
                      "dark = trustworthy, pale = do not interpret",
            transform=ax.transAxes, fontsize=7.5, color=MUTED, va="bottom")
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    # ------------------------------------------------------------- panel B
    ax = fig.add_subplot(gs[0, 1])
    lin2 = [l for l in ORDER if l in pooled.index]
    y = np.arange(len(lin2))[::-1]
    g = pooled.loc[lin2, "pooled_g"].to_numpy()
    lo = pooled.loc[lin2, "ci_lo"].to_numpy()
    hi = pooled.loc[lin2, "ci_hi"].to_numpy()
    cols = [DOWN if v < 0 else UP for v in g]
    alphas = [1.0 if l in RELIABLE else 0.33 for l in lin2]
    for yi, gi, li, hi_, c, a in zip(y, g, lo, hi, cols, alphas):
        ax.hlines(yi, li, hi_, color=c, lw=2, alpha=0.5 * a)
        ax.scatter(gi, yi, s=52, color=c, alpha=a, edgecolor=SURFACE,
                   linewidth=1.6, zorder=3)
    ax.axvline(0, color=BASE, lw=1)
    labels = [l if l in RELIABLE else f"{l}  (unvalidated)" for l in lin2]
    ax.set_yticks(y, labels, fontsize=8.5)
    ax.set_xlabel("pooled $g$ for the estimated proportion")
    ax.set_xlim(-1.45, 1.05)
    ax.set_title("B   Is the cell-type proportion different?",
                 loc="left", fontsize=11, fontweight="bold", pad=34)
    ax.text(0, 1.012, "5 studies, 97 AS vs 70 HC, random effects\n"
                      "lower in AS  ←                  →  higher in AS",
            transform=ax.transAxes, fontsize=7.5, color=MUTED, va="bottom")
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    # ------------------------------------------------------------- panel C
    ax = fig.add_subplot(gs[0, 2])
    items = [
        ("NK effector cargo\n(raw)", EFFECTOR_RAW),
        ("NK effector cargo\nadjusted for NK proportion", EFFECTOR_ADJ),
        ("NK proportion\n(disjoint identity genes)", NK_PROP),
    ]
    y = np.arange(len(items))[::-1]
    for yi, (lab, (gi, li, hi_)) in zip(y, items):
        ax.hlines(yi, li, hi_, color=DOWN, lw=2.4, alpha=0.5)
        ax.scatter(gi, yi, s=68, color=DOWN, edgecolor=SURFACE, linewidth=1.8,
                   zorder=3)
        ax.text(gi, yi + 0.21, f"{gi:+.2f}", ha="center", fontsize=8.5, color=INK2)
    ax.axvline(0, color=BASE, lw=1)
    ax.set_yticks(y, [i[0] for i in items], fontsize=8.5)
    ax.set_xlim(-1.62, 0.30)
    ax.set_ylim(-0.85, len(items) - 0.35)
    ax.set_xlabel("pooled $g$")
    ax.set_title("C   Composition or regulation?",
                 loc="left", fontsize=11, fontweight="bold", pad=34)
    ax.text(0, 1.012, "if composition explained it, the adjusted effect would\n"
                      "collapse toward zero. It retains 82%.",
            transform=ax.transAxes, fontsize=7.5, color=MUTED, va="bottom")
    ax.annotate("", xy=(EFFECTOR_ADJ[0], y[1] + 0.40), xytext=(EFFECTOR_RAW[0], y[0] - 0.40),
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.1,
                                connectionstyle="arc3,rad=-0.25"))
    ax.text(-1.33, (y[0] + y[1]) / 2, "only 18%\nlost", fontsize=7.5,
            color=MUTED, ha="center", va="center")
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    fig.suptitle(
        "Deconvolving the AS blood cohorts with a 312,000-cell single-cell reference",
        x=0.075, ha="left", fontsize=13.5, fontweight="bold", y=0.955,
    )
    fig.text(0.075, 0.885,
             "NK cells are modestly fewer in AS, but most of the cytotoxic deficit is not explained by their number — "
             "it is reduced output per cell.",
             ha="left", fontsize=8.5, color=INK2)

    out = FIG / "deconvolution.png"
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
