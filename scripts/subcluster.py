"""Reproduce the paper's actual design: sub-cluster compartments, then compare
the PROPORTIONS of sub-populations between responders and non-responders.

Song et al. Fig 2 sub-clusters 52,218 monocytes into 13 states and finds the
ISG-high non-classical monocyte fraction elevated in non-responders (p=0.0322).
Fig 7 sub-clusters 17,777 memory CD4 T cells and finds activated Th17 elevated
(p=0.049) and MAF-low Treg elevated (p=0.007).

Averaging a signature over a whole compartment -- what we did before -- cannot
see a 1-3% subpopulation. This re-clusters within compartments instead.

    modal run --detach scripts/subcluster.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-subcluster")
image = modal_app.image.add_local_python_source("modal_app")

IN = "/data/annotated_clean/as_annotated.h5ad"
OUT = "/data/subclusters"

MONO = ["Classical monocytes", "Non-classical monocytes"]
CD4 = ["Tcm/Naive helper T cells", "Tem/Effector helper T cells",
       "Regulatory T cells"]

# Marker panels used to name sub-clusters after the fact.
ISG = ["ISG15", "IFI6", "IFI44L", "MX1", "MX2", "OAS1", "OAS3", "IFIT1",
       "IFIT3", "RSAD2", "IFI44", "STAT1", "IRF7", "LY6E", "EPSTI1", "HERC5"]
TH17 = ["RORA", "RORC", "CCR6", "KLRB1", "IL23R", "CCR4", "IL17A", "IL17F",
        "ADAM12", "CTSH"]
ACTIV = ["CD69", "JUN", "FOS", "FOSB", "DUSP1", "NR4A1"]
TREG = ["FOXP3", "IL2RA", "IKZF2", "CTLA4", "TIGIT"]
NCMO = ["FCGR3A", "CDKN1C", "MS4A7", "LST1", "TCF7L2"]
CMO = ["CD14", "S100A8", "S100A9", "S100A12", "VCAN", "LYZ"]


def _one_per_patient(df):
    """One baseline row per patient; prefer the TNF draw where both exist."""
    d = df[df["Stage"] == "base"].copy()
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
    pd.set_option("display.width", 260)
    pd.set_option("display.max_rows", 300)

    full = ad.read_h5ad(IN)
    full.obs["sample_id"] = (full.obs["Patient"].astype(str) + "|" +
                             full.obs["Drug"].astype(str) + "|" +
                             full.obs["Stage"].astype(str))
    print(f"loaded {full.n_obs:,} cells")

    results = {}

    def subcluster(name, types, resolution, panels):
        print("\n" + "=" * 76)
        print(f"SUB-CLUSTERING: {name}")
        print("=" * 76)
        sub = full[full.obs["ct_fine"].isin(types)].copy()
        print(f"  {sub.n_obs:,} cells from {types}")

        # Re-do HVG/PCA WITHIN the compartment: the genes that distinguish
        # monocyte states are not the genes that distinguish monocytes from
        # T cells, so the global HVGs are the wrong basis here.
        sub.X = sub.layers["counts"].copy()
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sc.pp.highly_variable_genes(sub, n_top_genes=2000, batch_key="Sample_num")
        hv = sub[:, sub.var["highly_variable"]].copy()
        sc.pp.scale(hv, max_value=10)
        sc.tl.pca(hv, n_comps=30, svd_solver="arpack")
        sub.obsm["X_pca"] = hv.obsm["X_pca"]
        del hv

        import harmonypy
        ho = harmonypy.run_harmony(sub.obsm["X_pca"], sub.obs, ["Sample_num"])
        Z = ho.Z_corr
        if hasattr(Z, "detach"):
            Z = Z.detach().cpu().numpy()
        Z = np.asarray(Z, dtype=np.float32)
        if Z.shape[0] != sub.n_obs:
            Z = Z.T
        sub.obsm["X_harmony"] = np.ascontiguousarray(Z)

        sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_harmony")
        sc.tl.leiden(sub, resolution=resolution, key_added="sub",
                     flavor="igraph", n_iterations=2, directed=False)
        sc.tl.umap(sub)
        print(f"  {sub.obs['sub'].nunique()} sub-clusters")

        for pname, genes in panels.items():
            present = [g for g in genes if g in sub.var_names]
            sc.tl.score_genes(sub, present, score_name=pname, ctrl_size=50)

        prof = (sub.obs.groupby("sub", observed=True)[list(panels)]
                       .mean().round(3))
        prof["cells"] = sub.obs["sub"].value_counts()
        prof["pct"] = (100 * prof["cells"] / sub.n_obs).round(2)
        print("\n  sub-cluster marker profile:")
        print(prof.to_string())

        # proportion of each sub-cluster within the compartment, per sample
        prop = pd.crosstab(sub.obs["sample_id"], sub.obs["sub"], normalize="index")
        meta = (sub.obs.groupby("sample_id", observed=True)
                       .agg(Patient=("Patient", "first"), Drug=("Drug", "first"),
                            Stage=("Stage", "first"),
                            TNFi_Respond=("TNFi_Respond", "first"),
                            n=("sub", "size")))
        tab = prop.join(meta)
        results[name] = (sub, prof, tab)

        print("\n  --- proportions at baseline, non-responder vs responder ---")
        rows = []
        for cohort in ["A", "B"]:
            d = _one_per_patient(tab)
            if cohort == "A":
                d = d[d["Drug"] == "TNF"]
            bad = d[d["TNFi_Respond"] == "Bad"]
            good = d[d["TNFi_Respond"] == "Good"]
            for cl in prop.columns:
                if len(bad) < 3 or len(good) < 3:
                    continue
                u, p = stats.mannwhitneyu(bad[cl], good[cl],
                                          alternative="two-sided")
                rows.append({"cohort": cohort, "sub": cl,
                             "n_bad": len(bad), "n_good": len(good),
                             "pct_bad": 100 * bad[cl].mean(),
                             "pct_good": 100 * good[cl].mean(),
                             "p": p,
                             **{k: prof.loc[cl, k] for k in panels}})
        res = pd.DataFrame(rows)
        for cohort in ["A", "B"]:
            c = res[res["cohort"] == cohort].sort_values("p")
            print(f"\n  cohort {cohort}:")
            print(c.round(4).to_string(index=False))
        res.to_csv(f"{OUT}/{name}_tests.csv", index=False)
        return sub, prof, tab, res

    mono = subcluster("monocytes", MONO, 0.8,
                      {"ISG": ISG, "ncMo": NCMO, "cMo": CMO})
    cd4 = subcluster("cd4", CD4, 0.8,
                     {"TH17": TH17, "ACTIV": ACTIV, "TREG": TREG, "ISG": ISG})

    # ---------------- figures ----------------
    for name, (sub, prof, tab, res) in [("monocytes", mono), ("cd4", cd4)]:
        panels = [c for c in prof.columns if c not in ("cells", "pct")]
        n_panel = len(panels)
        fig, axes = plt.subplots(2, max(n_panel, 3), figsize=(5.2 * max(n_panel, 3), 10))

        xy = sub.obsm["X_umap"]
        ax = axes[0, 0]
        for j, cl in enumerate(sorted(sub.obs["sub"].unique(), key=int)):
            m = (sub.obs["sub"] == cl).values
            ax.scatter(xy[m, 0], xy[m, 1], s=0.3,
                       color=plt.get_cmap("tab20")(j % 20), rasterized=True)
            ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), cl, fontsize=11,
                    ha="center", va="center",
                    bbox=dict(fc="white", alpha=0.7, ec="none", pad=0.6))
        ax.set_title(f"{name}: {sub.obs['sub'].nunique()} sub-clusters")
        ax.set_xticks([]); ax.set_yticks([])

        for i, pname in enumerate(panels):
            ax = axes[0, i + 1] if i + 1 < axes.shape[1] else None
            if ax is None:
                break
            v = sub.obs[pname].values
            o = np.argsort(v)
            s = ax.scatter(xy[o, 0], xy[o, 1], s=0.3, c=v[o], cmap="magma",
                           rasterized=True)
            fig.colorbar(s, ax=ax)
            ax.set_title(f"{pname} score")
            ax.set_xticks([]); ax.set_yticks([])

        # boxplots for the three most significant sub-clusters (cohort B)
        top = res[res["cohort"] == "B"].nsmallest(3, "p")
        d = _one_per_patient(tab)
        for i, (_, r) in enumerate(top.iterrows()):
            if i >= axes.shape[1]:
                break
            ax = axes[1, i]
            cl = r["sub"]
            groups = [100 * d[d["TNFi_Respond"] == g][cl].values
                      for g in ["Good", "Bad"]]
            bp = ax.boxplot(groups, tick_labels=["Responder", "Non-responder"],
                            widths=0.55, patch_artist=True, showfliers=False)
            for patch, col in zip(bp["boxes"], ["#4C72B0", "#C44E52"]):
                patch.set_facecolor(col); patch.set_alpha(0.45)
            for k, (vals, col) in enumerate(zip(groups, ["#4C72B0", "#C44E52"])):
                ax.scatter(np.random.normal(k + 1, 0.06, len(vals)), vals,
                           color=col, s=34, zorder=3, edgecolor="white", lw=0.6)
            ax.set_ylabel("% of compartment")
            ax.set_title(f"sub-cluster {cl}   p = {r['p']:.4f}")
            ax.grid(axis="y", alpha=0.3)
        for j in range(len(top), axes.shape[1]):
            axes[1, j].axis("off")

        fig.tight_layout()
        fig.savefig(f"{OUT}/{name}.png", dpi=110)
        print(f"  wrote {OUT}/{name}.png")

    vol.commit()


@app.local_entrypoint()
def main():
    run.remote()
