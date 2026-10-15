"""Figure for the AS-vs-control bulk meta-analysis.

One colour language across all three panels: blue = lower in AS, red = higher
in AS, grey = the genome-wide background. Position carries the number; colour
only restates its sign, so nothing depends on colour alone.

Usage:  python scripts/bulk_figure.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "bulk" / "meta"
FIG = ROOT / "figures"

# validated diverging pair + chrome (see dataviz reference palette)
DOWN = "#2a78d6"
UP = "#e34948"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"

DIVERGING = LinearSegmentedColormap.from_list("as_div", [DOWN, "#f0efec", UP])

mpl.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
        "text.color": INK,
        "axes.labelcolor": INK2,
        "axes.edgecolor": BASE,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
    }
)

LINEAGE_ORDER = [
    "NK", "T cell", "B cell", "monocyte", "erythroid", "neutrophil", "platelet",
]

# Genes shown in the forest plot: the lineage markers that drive the result,
# plus the panel carried over from the single-cell project.
FOREST = {
    "NK / cytotoxic": ["EOMES", "GNLY", "NKG7", "PRF1", "KLRD1", "NCR3"],
    "T cell": ["CD247", "CD2", "CD3G", "IL7R", "LCK"],
    "monocyte": ["CSF1R", "CDKN1C", "MS4A7", "CD14", "LYZ", "VCAN"],
    "type-1 IFN": ["MX2", "IFI27", "IRF7", "ISG15", "RSAD2", "IFIT3"],
    "the paper's marker": ["AIF1"],
}

STUDIES = ["GSE73754", "GSE25101", "GSE18781", "GSE299639", "GSE181364"]
STUDY_N = {
    "GSE73754": "52/20", "GSE25101": "16/16", "GSE18781": "18/25",
    "GSE299639": "6/6", "GSE181364": "5/3",
}


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    m = pd.read_csv(META / "meta_primary.tsv", sep="\t", index_col=0)
    sens = pd.read_csv(META / "sensitivity.tsv", sep="\t")
    bg_median = float(m["pooled_g"].median())

    fig = plt.figure(figsize=(13.5, 9.2))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1.0, 1.25], height_ratios=[1.0, 1.35],
        hspace=0.34, wspace=0.30, left=0.085, right=0.965, top=0.885, bottom=0.075,
    )

    # ---------------------------------------------------------------- panel A
    ax = fig.add_subplot(gs[0, 0])
    vals = [float(sens.loc[sens["subset"] == "all 5 studies", l].iloc[0]) for l in LINEAGE_ORDER]
    y = np.arange(len(LINEAGE_ORDER))[::-1]
    ax.barh(
        y, vals, height=0.56,
        color=[DOWN if v < 0 else UP for v in vals],
        edgecolor=SURFACE, linewidth=2,
    )
    ax.axvline(0, color=BASE, lw=1)
    ax.axvline(bg_median, color=MUTED, lw=1.2, ls=(0, (4, 3)))
    ax.annotate(
        f"genome-wide median {bg_median:+.2f}",
        xy=(bg_median - 0.03, y[0] + 0.45), ha="right", va="center",
        fontsize=7.5, color=MUTED,
    )
    for yi, v in zip(y, vals):
        # Labels for near-zero bars go on the far side of the baseline, where
        # there is room, instead of colliding with it.
        side = -1 if v < 0 else 1
        ax.text(
            v + side * 0.035, yi, f"{v:+.2f}",
            va="center", ha="left" if side > 0 else "right",
            fontsize=8.5, color=INK2,
        )
    ax.set_yticks(y, LINEAGE_ORDER, fontsize=9, color=INK)
    ax.set_xlim(-1.02, 0.52)
    ax.set_xlabel("mean pooled effect size $g$  (SD units)")
    ax.set_title(
        "A   Which blood lineage moves",
        loc="left", fontsize=11, fontweight="bold", pad=10,
    )
    ax.text(
        0, 1.015, "lower in AS  ←                     →  higher in AS",
        transform=ax.transAxes, fontsize=7.5, color=MUTED,
    )
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    # ---------------------------------------------------------------- panel B
    ax = fig.add_subplot(gs[0, 1])
    # A blank row between groups holds the group label, so it never lands on a
    # confidence interval.
    rows, labels, yy, heads = [], [], [], []
    cursor = 0.0
    for grp, genes in FOREST.items():
        present = [g for g in genes if g in m.index]
        if not present:
            continue
        cursor -= 1.0
        heads.append((grp, cursor))
        for g in present:
            cursor -= 1.0
            rows.append(m.loc[g])
            labels.append(g)
            yy.append(cursor)
    fo = pd.DataFrame(rows)
    yy = np.array(yy)
    colors = [DOWN if v < 0 else UP for v in fo["pooled_g"]]
    ax.hlines(yy, fo["ci_lo"], fo["ci_hi"], color=colors, lw=2, alpha=0.55)
    ax.scatter(
        fo["pooled_g"], yy, s=46, color=colors,
        edgecolor=SURFACE, linewidth=1.6, zorder=3,
    )
    ax.axvline(0, color=BASE, lw=1)
    ax.set_yticks(yy, labels, fontsize=8)
    ax.set_xlim(-1.5, 1.15)
    ax.set_ylim(cursor - 0.6, 0.2)
    for yi, q in zip(yy, fo["q"]):
        if q < 0.05:
            ax.text(1.12, yi, "q<.05", va="center", ha="right",
                    fontsize=6.8, color=INK2)
    for grp, yh in heads:
        ax.text(-1.47, yh, grp, ha="left", va="center",
                fontsize=7.5, color=MUTED, style="italic")
        if yh < -1.5:
            ax.axhline(yh + 0.5, color=GRID, lw=0.8)
    ax.set_xlabel("pooled $g$ with 95% CI")
    ax.set_title(
        "B   Per-gene pooled effect, 5 studies, 97 AS vs 70 HC",
        loc="left", fontsize=11, fontweight="bold", pad=10,
    )
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    # ---------------------------------------------------------------- panel C
    ax = fig.add_subplot(gs[1, 0])
    cols = [c for c in STUDIES if f"g_{c}" in m.columns]
    hm = m.loc[labels, [f"g_{c}" for c in cols]].astype(float)
    norm = TwoSlopeNorm(vmin=-1.8, vcenter=0.0, vmax=1.8)
    im = ax.imshow(hm.to_numpy(), cmap=DIVERGING, norm=norm, aspect="auto")
    ax.set_xticks(
        range(len(cols)),
        [f"{c}\n{STUDY_N[c]}" for c in cols],
        fontsize=7.5, color=INK2,
    )
    ax.set_yticks(range(len(labels)), labels, fontsize=7.5)
    for i in range(len(labels)):
        for j in range(len(cols)):
            v = hm.iat[i, j]
            ax.text(
                j, i, "n/a" if not np.isfinite(v) else f"{v:+.1f}",
                ha="center", va="center", fontsize=6.2,
                color=INK if abs(v if np.isfinite(v) else 0) < 1.0 else SURFACE,
            )
    ax.set_title(
        "C   Same genes, each study separately  (AS/HC n below)",
        loc="left", fontsize=11, fontweight="bold", pad=10,
    )
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("$g$ within study", fontsize=8, color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, color=MUTED)

    # ---------------------------------------------------------------- panel D
    ax = fig.add_subplot(gs[1, 1])
    subsets = sens["subset"].tolist()
    x = np.arange(len(subsets))
    shown = ["NK", "T cell", "B cell", "monocyte", "neutrophil"]
    # Position carries the value; every line is directly labelled, so the
    # lineage identity never rests on colour.
    for lin in shown:
        v = sens[lin].to_numpy(float)
        col = DOWN if np.mean(v) < 0 else UP
        ax.plot(x, v, marker="o", ms=5.5, lw=2, color=col,
                markeredgecolor=SURFACE, markeredgewidth=1.4, alpha=0.9)
        ax.annotate(
            f" {lin}", xy=(x[-1], v[-1]), fontsize=8, color=INK2,
            va="center", ha="left",
        )
    ax.axhline(0, color=BASE, lw=1)
    ax.set_xticks(x, [s.replace("drop ", "−") for s in subsets],
                  rotation=32, ha="right", fontsize=7.5)
    ax.set_ylabel("mean pooled $g$")
    ax.set_xlim(-0.4, len(subsets) - 0.1)
    ax.margins(x=0.22)
    ax.set_title(
        "D   Leave-one-study-out: does any single study carry it?",
        loc="left", fontsize=11, fontweight="bold", pad=10,
    )
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    fig.suptitle(
        "Ankylosing spondylitis vs healthy controls in blood: random-effects meta-analysis of 5 studies",
        x=0.085, ha="left", fontsize=13.5, fontweight="bold", y=0.965,
    )
    fig.text(
        0.085, 0.925,
        "Effect size within each study, pooled by DerSimonian–Laird. "
        "Blue = lower in AS, red = higher in AS. "
        "The NK/cytotoxic decrease is the only lineage-scale signal, and it survives dropping any study.",
        ha="left", fontsize=8.5, color=INK2,
    )

    out = FIG / "bulk_meta.png"
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
