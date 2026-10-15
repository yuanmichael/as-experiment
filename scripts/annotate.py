"""Reference-based cell type annotation (CellTypist) + a clustering resolution sweep.

Two things at once:
  1. Project our cells onto a published annotated immune atlas to inherit labels.
  2. Re-cluster at several resolutions, and measure how well each agrees with
     those labels -- i.e. let the reference tell us what resolution to use,
     instead of me picking 1.0 by convention.

    modal run scripts/annotate.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-annotate")
image = modal_app.image.add_local_python_source("modal_app")

IN = "/data/processed/as_processed.h5ad"
OUT = "/data/annotated"
ANNOTATED = f"{OUT}/as_annotated.h5ad"

RESOLUTIONS = [0.2, 0.4, 0.6, 1.0, 1.5, 2.0, 3.0]


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=8.0, timeout=10800)
def annotate():
    import os

    import anndata as ad
    import celltypist
    from celltypist import models
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    os.makedirs(OUT, exist_ok=True)
    sc.settings.n_jobs = 8
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)

    adata = ad.read_h5ad(IN)
    print(f"loaded {adata.n_obs:,} x {adata.n_vars:,}")

    # CellTypist expects log1p data normalised to 1e4 per cell -- exactly what
    # adata.X already holds, so no re-normalisation.
    models.download_models(force_update=False,
                           model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])

    for model_name, key in [("Immune_All_High.pkl", "ct_coarse"),
                            ("Immune_All_Low.pkl", "ct_fine")]:
        print(f"\n--- {model_name} ---")
        pred = celltypist.annotate(adata, model=model_name,
                                   majority_voting=True, over_clustering="leiden")
        adata.obs[key] = pred.predicted_labels["majority_voting"].values
        adata.obs[f"{key}_raw"] = pred.predicted_labels["predicted_labels"].values
        adata.obs[f"{key}_conf"] = pred.probability_matrix.max(axis=1).values
        n = adata.obs[key].nunique()
        print(f"  {n} distinct labels, median confidence "
              f"{adata.obs[f'{key}_conf'].median():.3f}")

    print("\n=== coarse labels ===")
    vc = adata.obs["ct_coarse"].value_counts()
    print((pd.DataFrame({"cells": vc, "pct": (100 * vc / adata.n_obs).round(2)})
           ).to_string())

    print("\n=== fine labels (top 40) ===")
    vc = adata.obs["ct_fine"].value_counts()
    print((pd.DataFrame({"cells": vc, "pct": (100 * vc / adata.n_obs).round(2)})
           ).head(40).to_string())

    # ---- resolution sweep ----
    # The neighbour graph from the Harmony embedding is already in the object.
    print("\n=== resolution sweep (agreement with reference labels) ===")
    rows = []
    for res in RESOLUTIONS:
        key = f"leiden_r{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        lab = adata.obs["ct_fine"].astype(str)
        clu = adata.obs[key].astype(str)
        ari = adjusted_rand_score(lab, clu)
        nmi = normalized_mutual_info_score(lab, clu)
        # purity: within a cluster, what fraction share the commonest label?
        pur = (adata.obs.groupby(key, observed=True)["ct_fine"]
                        .apply(lambda s: s.value_counts(normalize=True).iloc[0])
                        .median())
        rows.append({"resolution": res, "clusters": adata.obs[key].nunique(),
                     "ARI": round(ari, 4), "NMI": round(nmi, 4),
                     "median_label_purity": round(pur, 4)})
        print(f"  res={res:<4} clusters={rows[-1]['clusters']:>3}  "
              f"ARI={ari:.4f}  NMI={nmi:.4f}  purity={pur:.3f}")
    sweep = pd.DataFrame(rows)

    # ---- figure ----
    xy = adata.obsm["X_umap_harmony"]
    fig = plt.figure(figsize=(22, 15))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])

    for i, (key, title) in enumerate([("ct_coarse", "CellTypist — coarse"),
                                      ("ct_fine", "CellTypist — fine")]):
        ax = fig.add_subplot(gs[0, i])
        vals = adata.obs[key].astype(str)
        cats = list(vals.value_counts().index)
        cmap = plt.get_cmap("tab20")
        for j, c in enumerate(cats):
            m = (vals == c).values
            ax.scatter(xy[m, 0], xy[m, 1], s=0.1, color=cmap(j % 20),
                       rasterized=True)
            # label at the median position of the group
            if m.sum() > 400:
                ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), c,
                        fontsize=7.5, ha="center", va="center",
                        bbox=dict(fc="white", ec="none", alpha=0.65, pad=0.9))
        ax.set_title(f"{title}  ({len(cats)} labels)")
        ax.set_xticks([]); ax.set_yticks([])

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(sweep["resolution"], sweep["ARI"], "o-", label="ARI")
    ax.plot(sweep["resolution"], sweep["NMI"], "s-", label="NMI")
    ax.plot(sweep["resolution"], sweep["median_label_purity"], "^-",
            label="median label purity")
    ax.set_xlabel("Leiden resolution"); ax.set_ylabel("agreement with reference")
    ax.legend(); ax.grid(alpha=0.3)
    ax.set_title("which resolution best matches the reference taxonomy?")

    ax = fig.add_subplot(gs[1, 1])
    ax.plot(sweep["resolution"], sweep["clusters"], "o-", color="firebrick")
    ax.set_xlabel("Leiden resolution"); ax.set_ylabel("number of clusters")
    ax.grid(alpha=0.3); ax.set_title("cluster count vs resolution")

    fig.tight_layout()
    fig.savefig(f"{OUT}/annotation.png", dpi=110)

    sweep.to_csv(f"{OUT}/resolution_sweep.csv", index=False)
    adata.write_h5ad(ANNOTATED, compression="gzip")
    vol.commit()
    print(f"\nwrote {ANNOTATED}")

    with open(f"{OUT}/annotation.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = annotate.remote()
    with open("figures/annotation.png", "wb") as f:
        f.write(png)
    print(f"saved figures/annotation.png ({len(png) / 1e6:.1f} MB)")
