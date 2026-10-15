"""Test the paper's Fig 7b claim: activated Th17 and MAF-low Treg elevated in
non-responders, as a proportion of MEMORY CD4 T cells.

Differences from the earlier attempt:
  - NO Harmony on Sample_num. Between-patient variation is the signal here, and
    correcting it away inside the compartment deletes what we are testing.
  - Proportions computed within memory CD4 only (their denominator), after
    identifying and excluding naive clusters by CCR7/SELL/TCF7/LEF1.
  - Both cohort A and cohort B.
  - Plus a clustering-free threshold test, as we did for the ISG signal.

    modal run --detach scripts/cd4_nohrm.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-cd4")
image = modal_app.image.add_local_python_source("modal_app")

IN = "/data/annotated_clean/as_annotated.h5ad"
OUT = "/data/cd4"

CD4 = ["Tcm/Naive helper T cells", "Tem/Effector helper T cells",
       "Regulatory T cells"]

PANELS = {
    "TH17": ["RORA", "RORC", "CCR6", "KLRB1", "IL23R", "CCR4", "IL17A",
             "IL17F", "ADAM12", "CTSH", "IL22"],
    "ACTIV": ["CD69", "JUN", "FOS", "FOSB", "DUSP1", "NR4A1", "JUNB"],
    "TREG": ["FOXP3", "IL2RA", "IKZF2", "CTLA4", "TIGIT"],
    "MAF": ["MAF"],
    "NAIVE": ["CCR7", "SELL", "TCF7", "LEF1", "NOSIP"],
    "TH1": ["TBX21", "IFNG", "CXCR3"],
}


def _one_per_patient(df, cohort):
    d = df[df["Stage"] == "base"].copy()
    if cohort == "A":
        d = d[d["Drug"] == "TNF"]
    d["_p"] = (d["Drug"] == "TNF").astype(int)
    return (d.sort_values("_p", ascending=False)
             .drop_duplicates("Patient", keep="first").drop(columns="_p"))


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=8.0, timeout=10800)
def run():
    import os

    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import stats

    os.makedirs(OUT, exist_ok=True)
    sc.settings.n_jobs = 8
    pd.set_option("display.width", 270)
    pd.set_option("display.max_rows", 300)

    full = ad.read_h5ad(IN)
    sub = full[full.obs["ct_fine"].isin(CD4)].copy()
    del full
    sub.obs["sample_id"] = (sub.obs["Patient"].astype(str) + "|" +
                            sub.obs["Drug"].astype(str) + "|" +
                            sub.obs["Stage"].astype(str))
    print(f"CD4 compartment: {sub.n_obs:,} cells")

    sub.X = sub.layers["counts"].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    sc.pp.highly_variable_genes(sub, n_top_genes=2000, batch_key="Sample_num")
    hv = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=30, svd_solver="arpack")
    sub.obsm["X_pca"] = hv.obsm["X_pca"]
    del hv

    # No Harmony: cluster on the uncorrected PCs.
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(sub, resolution=1.0, key_added="sub",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(sub)
    print(f"{sub.obs['sub'].nunique()} sub-clusters (no batch correction)")

    for name, genes in PANELS.items():
        present = [g for g in genes if g in sub.var_names]
        sc.tl.score_genes(sub, present, score_name=name, ctrl_size=50)

    prof = sub.obs.groupby("sub", observed=True)[list(PANELS)].mean().round(3)
    prof["cells"] = sub.obs["sub"].value_counts()
    prof["pct"] = (100 * prof["cells"] / sub.n_obs).round(2)
    # how patient-specific is each cluster? high => batch artefact, not a state
    prof["max_patient"] = (sub.obs.groupby("sub", observed=True)["Patient"]
                              .apply(lambda s: s.value_counts(normalize=True).iloc[0])
                              .round(3))
    print("\n=== sub-cluster profile ===")
    print(prof.sort_values("TH17", ascending=False).to_string())

    # naive clusters: high NAIVE score and low activation
    naive_cl = prof.index[(prof["NAIVE"] > prof["NAIVE"].median()) &
                          (prof["TH17"] < prof["TH17"].median())]
    print(f"\nclusters called naive-like: {list(naive_cl)}")
    sub.obs["memory"] = ~sub.obs["sub"].isin(naive_cl)
    mem = sub[sub.obs["memory"]].copy()
    print(f"memory CD4: {mem.n_obs:,} cells "
          f"(paper's Fig 7 used 17,777)")

    meta = (sub.obs.groupby("sample_id", observed=True)
                   .agg(Patient=("Patient", "first"), Drug=("Drug", "first"),
                        Stage=("Stage", "first"),
                        TNFi_Respond=("TNFi_Respond", "first")))

    def test_props(obs, label, denom_name):
        prop = pd.crosstab(obs["sample_id"], obs["sub"], normalize="index")
        n = obs.groupby("sample_id", observed=True).size()
        tab = prop.join(meta)
        tab["n"] = n
        tab = tab[tab["n"] >= 50]
        rows = []
        for cohort in ["A", "B"]:
            d = _one_per_patient(tab.reset_index(), cohort)
            bad = d[d["TNFi_Respond"] == "Bad"]
            good = d[d["TNFi_Respond"] == "Good"]
            if len(bad) < 3 or len(good) < 3:
                continue
            for cl in prop.columns:
                if cl not in d.columns:
                    continue
                _, p = stats.mannwhitneyu(bad[cl], good[cl],
                                          alternative="two-sided")
                rows.append({"cohort": cohort, "sub": cl,
                             "n_bad": len(bad), "n_good": len(good),
                             "pct_bad": 100 * bad[cl].mean(),
                             "pct_good": 100 * good[cl].mean(),
                             "p": p,
                             **{k: prof.loc[cl, k] for k in
                                ["TH17", "ACTIV", "TREG", "MAF", "NAIVE"]}})
        res = pd.DataFrame(rows)
        print(f"\n=== {label} (denominator: {denom_name}) ===")
        for cohort in ["A", "B"]:
            c = res[res["cohort"] == cohort].sort_values("p")
            if len(c):
                print(f"\n  cohort {cohort}:")
                print(c.round(4).to_string(index=False))
        return res, tab

    res_mem, tab_mem = test_props(mem.obs, "sub-cluster proportions",
                                  "memory CD4")
    res_mem.to_csv(f"{OUT}/cd4_memory_tests.csv", index=False)

    # ---- clustering-free threshold test, mirroring the ISG analysis ----
    print("\n" + "=" * 74)
    print("THRESHOLD TEST: fraction of memory CD4 that is Th17-high / activated")
    print("=" * 74)
    rows = []
    for score in ["TH17", "TREG"]:
        for pct in [75, 90, 95, 99]:
            thr = mem.obs[score].quantile(pct / 100)
            hi = (mem.obs[score] > thr)
            frac = (hi.groupby(mem.obs["sample_id"], observed=True).mean()
                      .rename("hi").to_frame().join(meta))
            frac["n"] = mem.obs.groupby("sample_id", observed=True).size()
            frac = frac[frac["n"] >= 50]
            for cohort in ["A", "B"]:
                d = _one_per_patient(frac.reset_index(), cohort)
                bad = d[d["TNFi_Respond"] == "Bad"]["hi"]
                good = d[d["TNFi_Respond"] == "Good"]["hi"]
                if len(bad) < 3 or len(good) < 3:
                    continue
                _, p = stats.mannwhitneyu(bad, good, alternative="two-sided")
                rows.append({"score": score, "threshold": f"top {100 - pct}%",
                             "cohort": cohort, "n_bad": len(bad),
                             "n_good": len(good), "pct_bad": 100 * bad.mean(),
                             "pct_good": 100 * good.mean(), "p": p})
    thr_res = pd.DataFrame(rows)
    print(thr_res.round(4).to_string(index=False))
    thr_res.to_csv(f"{OUT}/cd4_threshold_tests.csv", index=False)

    # ---- figure ----
    xy = mem.obsm["X_umap"]
    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    ax = axes[0, 0]
    for j, cl in enumerate(sorted(mem.obs["sub"].unique(), key=int)):
        m = (mem.obs["sub"] == cl).values
        ax.scatter(xy[m, 0], xy[m, 1], s=0.4,
                   color=plt.get_cmap("tab20")(j % 20), rasterized=True)
        ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), cl, fontsize=10,
                ha="center", va="center",
                bbox=dict(fc="white", alpha=0.7, ec="none", pad=0.5))
    ax.set_title(f"memory CD4, no batch correction ({mem.n_obs:,} cells)")
    ax.set_xticks([]); ax.set_yticks([])

    for ax, name in zip(axes[0, 1:].tolist() + [axes[1, 0]],
                        ["TH17", "TREG", "ACTIV"]):
        v = mem.obs[name].values
        o = np.argsort(v)
        s = ax.scatter(xy[o, 0], xy[o, 1], s=0.4, c=v[o], cmap="magma",
                       rasterized=True)
        fig.colorbar(s, ax=ax); ax.set_title(f"{name} score")
        ax.set_xticks([]); ax.set_yticks([])

    top = res_mem[res_mem["cohort"] == "B"].nsmallest(2, "p")
    for i, (_, r) in enumerate(top.iterrows()):
        ax = axes[1, i + 1]
        d = _one_per_patient(tab_mem.reset_index(), "B")
        groups = [100 * d[d["TNFi_Respond"] == g][r["sub"]].values
                  for g in ["Good", "Bad"]]
        bp = ax.boxplot(groups, tick_labels=["Responder", "Non-responder"],
                        widths=0.55, patch_artist=True, showfliers=False)
        for patch, col in zip(bp["boxes"], ["#4C72B0", "#C44E52"]):
            patch.set_facecolor(col); patch.set_alpha(0.4)
        for k, (v, col) in enumerate(zip(groups, ["#4C72B0", "#C44E52"])):
            ax.scatter(np.random.normal(k + 1, 0.06, len(v)), v, color=col,
                       s=38, zorder=3, edgecolor="white", lw=0.6)
        ax.set_title(f"sub-cluster {r['sub']} (TH17={r['TH17']:.2f})\n"
                     f"p = {r['p']:.4f}")
        ax.set_ylabel("% of memory CD4"); ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(f"{OUT}/cd4.png", dpi=115)
    vol.commit()
    with open(f"{OUT}/cd4.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = run.remote()
    with open("figures/cd4.png", "wb") as f:
        f.write(png)
    print("saved figures/cd4.png")
