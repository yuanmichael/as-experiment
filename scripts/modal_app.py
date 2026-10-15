"""Modal app for the AS / TNF-inhibitor-response project (GSE277117).

The 5 GB processed AnnData object lives in a Modal Volume, not on the laptop.

    modal run scripts/modal_app.py::download   # one-time, ~5 GB into the Volume
    modal run scripts/modal_app.py::peek       # describe what is in it
"""
import modal

URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277117/"
       "suppl/GSE277117_All_cell.h5ad")
H5AD = "/data/GSE277117_All_cell.h5ad"

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("curl")
    .pip_install("scanpy==1.12.4", "anndata", "h5py", "pandas", "scikit-learn",
                 "pyarrow", "igraph", "leidenalg", "harmonypy", "celltypist")
)


@app.function(image=image, volumes={"/data": vol}, timeout=3600)
def download():
    """Stream the GEO file straight from NCBI into the Volume."""
    import os
    import subprocess

    if os.path.exists(H5AD):
        print(f"already present: {os.path.getsize(H5AD) / 1e9:.2f} GB")
        return

    subprocess.run(
        ["curl", "-sS", "--fail", "--retry", "5", "-C", "-", "-o", H5AD, URL],
        check=True,
    )
    vol.commit()
    print(f"downloaded {os.path.getsize(H5AD) / 1e9:.2f} GB -> {H5AD}")


@app.function(image=image, volumes={"/data": vol}, memory=4096, timeout=1800)
def peek():
    """Describe the object without loading the count matrix into memory."""
    import anndata as ad
    import pandas as pd

    pd.set_option("display.width", 200)

    # backed="r" leaves .X on disk; .obs/.var come into memory
    adata = ad.read_h5ad(H5AD, backed="r")
    print(f"shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"X on disk: {type(adata.X)}")

    print(f"\n--- obs columns ({len(adata.obs.columns)}) ---")
    for col in adata.obs.columns:
        s = adata.obs[col]
        n = s.nunique(dropna=True)
        head = list(pd.unique(s.dropna()))[:6]
        print(f"  {col:<28} {str(s.dtype):<10} {n:>6} unique  e.g. {head}")

    print(f"\n--- var columns ---")
    print(adata.var.head())

    print(f"\n--- obsm (embeddings) ---")
    for k, v in adata.obsm.items():
        print(f"  {k}: {v.shape}")

    print(f"\n--- layers / raw ---")
    print(f"  layers: {list(adata.layers.keys())}")
    print(f"  raw present: {adata.raw is not None}")

    print(f"\n--- uns keys ---")
    print(f"  {list(adata.uns.keys())}")


@app.function(image=image, volumes={"/data": vol}, memory=4096, timeout=1800)
def design():
    """Lay out the actual study design: who, which drug, which timepoint, what label."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 120)

    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs

    # Is .X raw integer counts, or already normalized?
    chunk = adata.X[:2000].toarray()
    print(f"X sample: min={chunk.min():.3f} max={chunk.max():.3f} dtype={chunk.dtype}")
    print(f"all integers? {np.allclose(chunk, np.round(chunk))}")
    print(f"median total per cell (first 2000): {np.median(chunk.sum(1)):,.0f}")

    print("\n=== cells per patient ===")
    print(obs["Patient"].value_counts().to_string())

    print("\n=== patient x (drug, stage, response) ===")
    tab = (obs.groupby(["Patient", "Drug", "Stage", "TNFi_Respond"], observed=True)
              .size().rename("cells").reset_index())
    print(tab.to_string(index=False))

    print("\n=== cohort definitions ===")
    base = obs[(obs["Patient"] != "Multiplet") & (obs["Stage"] == "base")]
    strict = base[base["Drug"] == "TNF"]
    print(f"  A (strict, TNF arm baseline):  {strict['Patient'].nunique()} patients")
    print(f"  B (expanded, all baselines):   {base['Patient'].nunique()} patients")

    print("\n=== the key count: patients with a usable TNFi label ===")
    real = obs[(obs["Patient"] != "Multiplet") & (obs["TNFi_Respond"] != "NA")]
    for stage in ["base", "1mo"]:
        sub = real[(real["Drug"] == "TNF") & (real["Stage"] == stage)]
        g = sub.groupby("TNFi_Respond", observed=True)["Patient"].nunique()
        print(f"  TNF / {stage:>4}: {dict(g)}  -> {sub['Patient'].nunique()} patients, "
              f"{len(sub):,} cells")
