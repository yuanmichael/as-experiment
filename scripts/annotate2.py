"""Re-annotate with CellTypist's own fine-grained over-clustering.

The first pass forced majority voting to use the 27-cluster res=1.0 Leiden
partition, which is far too coarse: a 536-cell population cannot be the
majority of an 11,600-cell cluster, so Th17/MAIT/gamma-delta/ILC3 were all
absorbed into their larger neighbours. Letting CellTypist pick its own
over-clustering (very high resolution, hundreds of small groups) denoises the
per-cell calls while leaving small populations intact.

Keeps both leiden_r0.6 and leiden_r1.0 for comparison.

    modal run --detach scripts/annotate2.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-annotate2")
image = modal_app.image.add_local_python_source("modal_app")

IN = "/data/annotated/as_annotated.h5ad"
OUT = "/data/annotated"
FINAL = f"{OUT}/as_annotated_v2.h5ad"

# Populations that matter mechanistically for spondyloarthritis (IL-17 axis).
WATCH = ["Type 17 helper T cells", "MAIT cells", "CRTAM+ gamma-delta T cells",
         "ILC3", "Regulatory T cells", "Tem/Effector helper T cells"]


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=8.0, timeout=10800)
def annotate2():
    import anndata as ad
    import celltypist
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc

    sc.settings.n_jobs = 8
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 250)

    adata = ad.read_h5ad(IN)
    print(f"loaded {adata.n_obs:,} x {adata.n_vars:,}")

    # No over_clustering argument -> CellTypist builds its own, at a resolution
    # scaled to dataset size (very fine for 314k cells).
    pred = celltypist.annotate(adata, model="Immune_All_Low.pkl",
                               majority_voting=True)
    adata.obs["ct_fine_v2"] = pred.predicted_labels["majority_voting"].values
    adata.obs["ct_oc"] = pred.predicted_labels["over_clustering"].values
    print(f"\nover-clusters used: {adata.obs['ct_oc'].nunique()}")
    print(f"median over-cluster size: "
          f"{adata.obs['ct_oc'].value_counts().median():,.0f} cells")
    print(f"labels after voting: {adata.obs['ct_fine_v2'].nunique()} "
          f"(was 14 with the coarse over-clustering)")

    print("\n=== label counts, v2 ===")
    vc = adata.obs["ct_fine_v2"].value_counts()
    print(pd.DataFrame({"cells": vc, "pct": (100 * vc / adata.n_obs).round(3)})
          .to_string())

    print("\n=== did the IL-17-axis populations survive? ===")
    raw = adata.obs["ct_fine_raw"].value_counts()
    v1 = adata.obs["ct_fine"].value_counts()
    v2 = adata.obs["ct_fine_v2"].value_counts()
    for pop in WATCH:
        print(f"  {pop:<32} raw={raw.get(pop, 0):>6,}  "
              f"v1={v1.get(pop, 0):>6,}  v2={v2.get(pop, 0):>6,}")

    print("\n=== where does low confidence live? ===")
    conf = (adata.obs.groupby("ct_fine_raw", observed=True)["ct_fine_conf"]
                     .agg(cells="size", median_conf="median",
                          frac_below_0_5=lambda s: (s < 0.5).mean())
                     .sort_values("cells", ascending=False).head(25).round(3))
    print(conf.to_string())

    print("\n=== does voting rescue the low-confidence cells? ===")
    low = adata.obs["ct_fine_conf"] < 0.5
    # .astype(str): both are Categorical with different category sets, which
    # pandas refuses to compare directly.
    agree_raw_v2 = (adata.obs["ct_fine_raw"].astype(str)
                    == adata.obs["ct_fine_v2"].astype(str))
    print(f"  cells below 0.5 confidence: {low.mean() * 100:.2f}%")
    print(f"  of those, voting CHANGED the label for "
          f"{(~agree_raw_v2[low]).mean() * 100:.1f}%")
    print(f"  of confident cells, voting changed "
          f"{(~agree_raw_v2[~low]).mean() * 100:.1f}%")

    # ---- figure: v2 labels, plus the watch populations highlighted ----
    xy = adata.obsm["X_umap_harmony"]
    fig, axes = plt.subplots(1, 2, figsize=(24, 11))

    ax = axes[0]
    vals = adata.obs["ct_fine_v2"].astype(str)
    cats = list(vals.value_counts().index)
    cmap = plt.get_cmap("tab20")
    for j, c in enumerate(cats):
        m = (vals == c).values
        ax.scatter(xy[m, 0], xy[m, 1], s=0.1, color=cmap(j % 20), rasterized=True)
        if m.sum() > 300:
            ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), c, fontsize=7,
                    ha="center", va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.8))
    ax.set_title(f"CellTypist v2 — fine over-clustering ({len(cats)} labels)")
    ax.set_xticks([]); ax.set_yticks([])

    ax = axes[1]
    ax.scatter(xy[:, 0], xy[:, 1], s=0.05, color="0.86", rasterized=True)
    hi = plt.get_cmap("tab10")
    for j, pop in enumerate(WATCH):
        m = (vals == pop).values
        if m.sum():
            ax.scatter(xy[m, 0], xy[m, 1], s=1.6, color=hi(j),
                       label=f"{pop} ({m.sum():,})", rasterized=True)
    ax.set_title("IL-17-axis populations")
    ax.legend(markerscale=8, fontsize=10, loc="best")
    ax.set_xticks([]); ax.set_yticks([])

    fig.tight_layout()
    fig.savefig(f"{OUT}/annotation_v2.png", dpi=110)

    adata.write_h5ad(FINAL, compression="gzip")
    vol.commit()
    print(f"\nwrote {FINAL}")

    with open(f"{OUT}/annotation_v2.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = annotate2.remote()
    with open("figures/annotation_v2.png", "wb") as f:
        f.write(png)
    print("saved figures/annotation_v2.png")
