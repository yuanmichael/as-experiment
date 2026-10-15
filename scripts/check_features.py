"""Are there non-gene features (CITE-seq antibody tags) in the feature list?

The genome-wide scan's top hits were CD44-1, CD162, HLA-DR, CD11a, CD1c --
antibody panel names, not gene symbols. If those are Antibody Capture features,
they are present only in the CITE-seq libraries and absent (zero) in the plain
scRNA ones, so they measure assay type, not biology.

    modal run scripts/check_features.py
"""
import modal

import modal_app
from modal_app import H5AD, vol

app = modal.App("as-featcheck")
image = modal_app.image.add_local_python_source("modal_app")

SUSPECTS = ["CD44-1", "CD162", "HLA-DR", "CD11a", "CD1c", "CD3-1", "CD19"]


@app.function(image=image, volumes={"/data": vol}, memory=16384, timeout=1800)
def check():
    import anndata as ad
    import numpy as np
    import pandas as pd

    pd.set_option("display.width", 220)

    raw = ad.read_h5ad(H5AD, backed="r")
    print("=== original file, var.feature_types ===")
    print(raw.var["feature_types"].value_counts().to_string())

    adt = raw.var_names[raw.var["feature_types"] != "Gene Expression"]
    print(f"\nnon-gene features: {len(adt)}")
    print(list(adt[:60]))

    print("\n=== are the scan's top hits in that list? ===")
    for s in SUSPECTS:
        print(f"  {s:<10} in ADT list: {s in set(adt)}")

    # How is an ADT feature distributed across libraries?
    pb = ad.read_h5ad("/data/pseudobulk/pseudobulk.h5ad")
    for s in ["CD44-1", "CD162"]:
        if s in pb.var_names:
            v = np.asarray(pb[:, s].X.todense()).ravel()
            df = pb.obs.assign(val=v)
            per = (df.groupby("sample_id", observed=True)
                     .agg(val=("val", "mean"),
                          resp=("TNFi_Respond", "first"),
                          drug=("Drug", "first")))
            nz = per["val"] > 0.5
            print(f"\n--- {s} ---")
            print(f"  samples with signal: {int(nz.sum())} of {len(per)}")
            print("  by response group:")
            print(per.assign(has=nz).groupby("resp", observed=True)["has"]
                     .agg(["size", "sum", "mean"]).round(3).to_string())
            print("  by drug arm:")
            print(per.assign(has=nz).groupby("drug", observed=True)["has"]
                     .agg(["size", "sum", "mean"]).round(3).to_string())


@app.local_entrypoint()
def main():
    check.remote()
