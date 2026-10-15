"""Filter -> normalize -> HVG -> PCA -> cluster, with and without batch correction.

Deliberately produces BOTH an uncorrected and a Harmony-corrected embedding, so
we can see whether cells cluster by cell type or by library before deciding to
trust the correction.

    modal run scripts/process.py
"""
import modal

import modal_app
from modal_app import H5AD, vol

app = modal.App("as-process")
image = modal_app.image.add_local_python_source("modal_app")

OUT = "/data/processed_clean"
PROCESSED = f"{OUT}/as_processed.h5ad"

# Agreed thresholds ("standard"): keep cells with 200-6000 genes and <=15% mito.
MIN_GENES, MAX_GENES, MAX_MT = 200, 6000, 15.0


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=8.0, timeout=10800)
def process():
    import os

    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc

    os.makedirs(OUT, exist_ok=True)
    sc.settings.n_jobs = 8
    pd.set_option("display.width", 250)

    ckpt = f"{OUT}/checkpoint_pass1.h5ad"
    if os.path.exists(ckpt):
        print(f"resuming from {ckpt}")
        adata = ad.read_h5ad(ckpt)
        return _pass2(adata)

    adata = ad.read_h5ad(H5AD)

    # CRITICAL: the matrix mixes 36,601 genes with 163 CITE-seq antibody-capture
    # features. Every non-responder sample was run with CITE-seq and most
    # responder samples were not, so antibody counts are ~7% of totals in one
    # group and ~1% in the other. Left in, they contaminate HVG selection and
    # skew CPM normalisation in a direction that mimics a biological result.
    # Filter by the feature_types column (by index -- several antibodies share a
    # name with a gene, so matching on names is unreliable).
    is_gene = (adata.var["feature_types"] == "Gene Expression").values
    print(f"features: {adata.n_vars:,} total, {int(is_gene.sum()):,} genes, "
          f"{int((~is_gene).sum()):,} antibody-capture (dropped)")
    adata = adata[:, is_gene].copy()
    adata.var_names_make_unique()

    adata = adata[adata.obs["Patient"] != "Multiplet"].copy()
    print(f"after dropping multiplets: {adata.n_obs:,} x {adata.n_vars:,}")

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)

    keep = (adata.obs["n_genes_by_counts"].between(MIN_GENES, MAX_GENES)
            & (adata.obs["pct_counts_mt"] <= MAX_MT))
    adata = adata[keep].copy()
    print(f"after cell QC: {adata.n_obs:,} cells")

    # A gene seen in almost no cells cannot inform clustering and only adds noise.
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"after gene filter (>=3 cells): {adata.n_vars:,} genes")

    # Keep the raw integer counts. Pseudobulk and differential expression must
    # start from counts, NOT from the log-normalised values used for clustering.
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)   # equalise sequencing depth
    sc.pp.log1p(adata)                             # ratios -> distances
    adata.raw = adata                              # log-norm all genes, for plotting

    # batch_key makes HVGs be chosen per library and pooled, so genes that are
    # variable only because of one bad library do not dominate the selection.
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="Sample_num")
    print(f"HVGs: {int(adata.var['highly_variable'].sum())}")

    hv = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=50, svd_solver="arpack")
    adata.obsm["X_pca"] = hv.obsm["X_pca"]
    del hv

    # ---- Pass 1: no batch correction ----
    print("\n[1/2] clustering WITHOUT batch correction...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden_raw",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)
    adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()
    print(f"  clusters: {adata.obs['leiden_raw'].nunique()}")

    adata.write_h5ad(ckpt, compression="gzip")
    vol.commit()
    print(f"  checkpointed to {ckpt}")
    return _pass2(adata)


def _pass2(adata):
    """Harmony correction, second clustering, diagnostics and figure."""
    import os

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc

    sc.settings.n_jobs = 8
    pd.set_option("display.width", 250)

    # ---- Pass 2: Harmony on library ----
    # Called directly rather than via sc.external: recent harmonypy returns a
    # tensor that scanpy's wrapper transposes into the wrong shape.
    print("\n[2/2] clustering WITH Harmony correction on Sample_num...")
    import harmonypy

    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["Sample_num"])
    Z = ho.Z_corr
    if hasattr(Z, "detach"):
        Z = Z.detach().cpu().numpy()
    Z = np.asarray(Z, dtype=np.float32)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    assert Z.shape == (adata.n_obs, adata.obsm["X_pca"].shape[1]), Z.shape
    adata.obsm["X_pca_harmony"] = np.ascontiguousarray(Z)

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep="X_pca_harmony")
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)
    adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
    print(f"  clusters: {adata.obs['leiden'].nunique()}")

    # ---- Did correction actually help? ----
    # A cluster made of cells from one library is a batch artifact. A cluster
    # drawing evenly from all 28 is a real cell type.
    print("\n=== how patient-pure is each cluster? (1.0 = one patient only) ===")
    for key in ["leiden_raw", "leiden"]:
        pur = (adata.obs.groupby(key, observed=True)["Patient"]
                        .apply(lambda s: s.value_counts(normalize=True).iloc[0]))
        n_bad = int((pur > 0.5).sum())
        print(f"  {key:<11} clusters={len(pur):>3}  median purity={pur.median():.3f}  "
              f"max={pur.max():.3f}  clusters >50% one patient: {n_bad}")

    print("\n=== singleplex-vs-pooled mixing per cluster ===")
    single = adata.obs["Sample_num"].astype(str).isin(
        [str(i) for i in range(1, 12)])
    adata.obs["batch_group"] = np.where(single, "singleplex", "pooled")
    overall = (adata.obs["batch_group"] == "singleplex").mean()
    print(f"  overall singleplex fraction: {overall:.3f}")
    for key in ["leiden_raw", "leiden"]:
        frac = adata.obs.groupby(key, observed=True)["batch_group"].apply(
            lambda s: (s == "singleplex").mean())
        print(f"  {key:<11} per-cluster singleplex fraction: "
              f"min={frac.min():.3f} median={frac.median():.3f} max={frac.max():.3f}")

    # ---- Figure ----
    fig, axes = plt.subplots(2, 3, figsize=(21, 13))
    for row, (emb, lab) in enumerate([("X_umap_raw", "leiden_raw"),
                                      ("X_umap_harmony", "leiden")]):
        xy = adata.obsm[emb]
        title = "uncorrected" if row == 0 else "Harmony-corrected"

        for col, (colour_by, name) in enumerate(
                [(lab, "cluster"), ("batch_group", "batch"), ("Patient", "patient")]):
            ax = axes[row, col]
            vals = adata.obs[colour_by].astype(str)
            cats = sorted(vals.unique())
            cmap = plt.get_cmap("tab20" if len(cats) <= 20 else "gist_ncar")
            for i, c in enumerate(cats):
                m = (vals == c).values
                ax.scatter(xy[m, 0], xy[m, 1], s=0.12,
                           color=cmap(i / max(len(cats) - 1, 1)), label=c,
                           rasterized=True)
            ax.set_title(f"{title} — by {name}")
            ax.set_xticks([]); ax.set_yticks([])
            if name == "batch":
                ax.legend(markerscale=30, fontsize=9, loc="best")
    fig.tight_layout()
    fig.savefig(f"{OUT}/umap_batch_check.png", dpi=110)

    adata.write_h5ad(PROCESSED, compression="gzip")
    vol.commit()
    print(f"\nwrote {PROCESSED} ({os.path.getsize(PROCESSED) / 1e9:.2f} GB)")

    with open(f"{OUT}/umap_batch_check.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = process.remote()
    with open("figures/umap_batch_check.png", "wb") as f:
        f.write(png)
    print(f"saved figures/umap_batch_check.png ({len(png) / 1e6:.1f} MB)")
