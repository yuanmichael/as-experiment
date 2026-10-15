"""Annotate the clean (antibody-free) object.

Same approach as before: CellTypist reference mapping with its own fine-grained
over-clustering, plus Leiden at the two resolutions we agreed to keep.

    modal run --detach scripts/annotate_clean.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-annotate-clean")
image = modal_app.image.add_local_python_source("modal_app")

IN = "/data/processed_clean/as_processed.h5ad"
OUT = "/data/annotated_clean"
FINAL = f"{OUT}/as_annotated.h5ad"

WATCH = ["Type 17 helper T cells", "MAIT cells", "CRTAM+ gamma-delta T cells",
         "ILC3", "Regulatory T cells", "Tem/Effector helper T cells"]


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=8.0, timeout=10800)
def run():
    import os

    import anndata as ad
    import celltypist
    import pandas as pd
    import scanpy as sc

    os.makedirs(OUT, exist_ok=True)
    sc.settings.n_jobs = 8
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 250)

    adata = ad.read_h5ad(IN)
    print(f"loaded {adata.n_obs:,} x {adata.n_vars:,}")
    assert "feature_types" not in adata.var or \
        (adata.var["feature_types"] == "Gene Expression").all(), \
        "antibody features present -- rebuild from process.py"

    for model, pre in [("Immune_All_Low.pkl", "ct_fine"),
                       ("Immune_All_High.pkl", "ct_coarse")]:
        print(f"\n--- {model} ---")
        pred = celltypist.annotate(adata, model=model, majority_voting=True)
        adata.obs[pre] = pred.predicted_labels["majority_voting"].values
        adata.obs[f"{pre}_raw"] = pred.predicted_labels["predicted_labels"].values
        adata.obs[f"{pre}_conf"] = pred.probability_matrix.max(axis=1).values
        print(f"  {adata.obs[pre].nunique()} labels after voting, "
              f"{adata.obs[f'{pre}_raw'].nunique()} before")

    print("\n=== fine labels ===")
    vc = adata.obs["ct_fine"].value_counts()
    print(pd.DataFrame({"cells": vc, "pct": (100 * vc / adata.n_obs).round(3)})
          .to_string())

    print("\n=== IL-17-axis populations ===")
    raw = adata.obs["ct_fine_raw"].value_counts()
    fin = adata.obs["ct_fine"].value_counts()
    for pop in WATCH:
        print(f"  {pop:<32} raw={raw.get(pop, 0):>6,}  final={fin.get(pop, 0):>6,}")

    for res in [0.6, 1.0]:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_r{res}",
                     flavor="igraph", n_iterations=2, directed=False)
        print(f"  leiden r={res}: {adata.obs[f'leiden_r{res}'].nunique()} clusters")

    adata.write_h5ad(FINAL, compression="gzip")
    vol.commit()
    print(f"\nwrote {FINAL}")


@app.local_entrypoint()
def main():
    run.remote()
