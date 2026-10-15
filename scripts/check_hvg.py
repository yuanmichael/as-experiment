"""Did the antibody-capture features contaminate the HVG selection, and so the
PCA, the clustering and the UMAP?

    modal run scripts/check_hvg.py
"""
import modal

import modal_app
from modal_app import H5AD, vol

app = modal.App("as-hvgcheck")
image = modal_app.image.add_local_python_source("modal_app")


@app.function(image=image, volumes={"/data": vol}, memory=16384, timeout=1800)
def check():
    import anndata as ad
    import numpy as np
    import pandas as pd

    raw = ad.read_h5ad(H5AD, backed="r")
    adt = set(raw.var_names[raw.var["feature_types"] == "Antibody Capture"])
    # var_names_make_unique appends -1 to the second copy of a duplicated name
    adt_variants = adt | {f"{a}-1" for a in adt}

    proc = ad.read_h5ad("/data/annotated/as_annotated_v2.h5ad", backed="r")
    hv = proc.var_names[proc.var["highly_variable"].values]
    hit = [g for g in hv if g in adt_variants]
    print(f"HVGs: {len(hv)}")
    print(f"HVGs that are antibody features: {len(hit)}")
    print(f"  {hit}")

    # How much signal do they carry? CPM share of ADT counts per group.
    pb = ad.read_h5ad("/data/pseudobulk/pseudobulk.h5ad")
    in_pb = [g for g in pb.var_names if g in adt_variants]
    print(f"\nantibody features present in pseudobulk: {len(in_pb)}")
    cnt = pb.layers["counts"]
    tot = np.asarray(cnt.sum(axis=1)).ravel()
    adt_tot = np.asarray(cnt[:, [pb.var_names.get_loc(g) for g in in_pb]]
                         .sum(axis=1)).ravel()
    frac = adt_tot / np.maximum(tot, 1)
    df = pb.obs.assign(adt_frac=frac)
    print("\nfraction of each group's total counts that are antibody tags:")
    print(df.groupby("TNFi_Respond", observed=True)["adt_frac"]
            .describe(percentiles=[.5, .9]).round(4).to_string())
    print("\nby drug arm:")
    print(df.groupby("Drug", observed=True)["adt_frac"]
            .describe(percentiles=[.5, .9]).round(4).to_string())


@app.local_entrypoint()
def main():
    check.remote()
