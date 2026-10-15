"""Composition at baseline, responders vs non-responders -- laid out like the
paper's Fig 1c so the two can be compared directly.

Maps our 23 CellTypist labels onto the paper's six coarse groups, then plots
every patient as a point with the Mann-Whitney p above each pair.

    modal run scripts/composition_fig.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-compfig")
image = modal_app.image.add_local_python_source("modal_app")

PB = "/data/pseudobulk_clean"
OUT = "/data/figures"

# our CellTypist labels -> the paper's Fig 1c groups
GROUPS = {
    "T cell": ["Tcm/Naive helper T cells", "Tem/Effector helper T cells",
               "Tem/Temra cytotoxic T cells", "Tem/Trm cytotoxic T cells",
               "Tcm/Naive cytotoxic T cells", "Regulatory T cells", "MAIT cells"],
    "NK cell": ["CD16+ NK cells", "CD16- NK cells", "NK cells", "ILC3"],
    "B cell": ["Naive B cells", "Memory B cells", "Age-associated B cells"],
    "CD14 M": ["Classical monocytes"],
    "CD16 M": ["Non-classical monocytes"],
    "DC": ["DC1", "DC2", "pDC"],
    "Plasma": ["Plasma cells"],
}
ORDER = ["T cell", "NK cell", "B cell", "CD14 M", "CD16 M", "DC", "Plasma"]


@app.function(image=image, volumes={"/data": vol}, memory=16384, timeout=1800)
def run():
    import os

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import stats

    os.makedirs(OUT, exist_ok=True)
    pd.set_option("display.width", 240)

    comp = pd.read_parquet(f"{PB}/composition.parquet")
    scores = pd.read_parquet(f"{PB}/scores.parquet")
    meta = scores.drop_duplicates("sample_id")[
        ["sample_id", "Patient", "Drug", "Stage", "TNFi_Respond"]]
    df = comp.reset_index().merge(meta, on="sample_id")

    for g, members in GROUPS.items():
        present = [m for m in members if m in comp.columns]
        df[g] = df[present].sum(axis=1)

    def one_per_patient(d, cohort):
        d = d[d["Stage"] == "base"].copy()
        if cohort == "A":
            d = d[d["Drug"] == "TNF"]
        d["_p"] = (d["Drug"] == "TNF").astype(int)
        return (d.sort_values("_p", ascending=False)
                 .drop_duplicates("Patient", keep="first"))

    for cohort in ["A", "B"]:
        d = one_per_patient(df, cohort)
        print(f"\n=== cohort {cohort}: "
              f"{(d['TNFi_Respond'] == 'Bad').sum()} NR vs "
              f"{(d['TNFi_Respond'] == 'Good').sum()} R ===")
        rows = []
        for g in ORDER:
            bad = d[d["TNFi_Respond"] == "Bad"][g]
            good = d[d["TNFi_Respond"] == "Good"][g]
            u, p = stats.mannwhitneyu(bad, good, alternative="two-sided")
            rows.append({"group": g, "pct_NR": 100 * bad.mean(),
                         "pct_R": 100 * good.mean(),
                         "delta": 100 * (bad.mean() - good.mean()), "p": p})
        print(pd.DataFrame(rows).round(4).to_string(index=False))

    # ---- figure: cohort B (and A) side by side ----
    fig, axes = plt.subplots(2, 1, figsize=(13, 11))
    for ax, cohort, title in [
            (axes[0], "B", "Cohort B (expanded): 13 non-responders vs 12 responders"),
            (axes[1], "A", "Cohort A (strict, TNF arm): 4 non-responders vs 12 responders")]:
        d = one_per_patient(df, cohort)
        R = d[d["TNFi_Respond"] == "Good"]
        NR = d[d["TNFi_Respond"] == "Bad"]
        width = 0.34
        for i, g in enumerate(ORDER):
            for k, (grp, col, lab) in enumerate(
                    [(R, "#4C72B0", "Responder"), (NR, "#C44E52", "Non-responder")]):
                x = i + (k - 0.5) * width
                v = 100 * grp[g].values
                ax.bar(x, v.mean(), width=width * 0.88, color=col, alpha=0.32,
                       edgecolor=col, linewidth=1.4,
                       label=lab if i == 0 else None)
                ax.errorbar(x, v.mean(), yerr=v.std() / np.sqrt(len(v)),
                            color=col, capsize=3, lw=1.4)
                ax.scatter(np.random.normal(x, 0.035, len(v)), v, s=26,
                           color=col, zorder=3, edgecolor="white", linewidth=0.5)
            bad = NR[g]; good = R[g]
            _, p = stats.mannwhitneyu(bad, good, alternative="two-sided")
            top = max(100 * bad.max(), 100 * good.max())
            ax.text(i, top * 1.08, f"p = {p:.3f}", ha="center",
                    fontsize=9.5,
                    fontweight="bold" if p < 0.05 else "normal",
                    color="black" if p < 0.05 else "0.45")
        ax.set_xticks(range(len(ORDER)))
        ax.set_xticklabels(ORDER)
        ax.set_ylabel("Proportion of PBMCs (%)")
        ax.set_title(title)
        ax.legend(frameon=False)
        ax.grid(axis="y", alpha=0.25)
        ax.set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(f"{OUT}/composition_baseline.png", dpi=125)
    vol.commit()

    with open(f"{OUT}/composition_baseline.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = run.remote()
    with open("figures/composition_baseline.png", "wb") as f:
        f.write(png)
    print("saved figures/composition_baseline.png")
