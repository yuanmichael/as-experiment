"""Collapse 314,368 cells into a patient-level table.

Produces three things, keyed by (Patient, Drug, Stage) x cell type:
  1. composition.parquet  -- fraction of each cell type per sample
  2. pseudobulk.h5ad      -- summed raw counts per sample x cell type, CPM+log
  3. scores.parquet       -- per-cell signature scores averaged per group

The signature is MSigDB Hallmark INTERFERON_ALPHA_RESPONSE: defined externally,
years before this study, so it cannot have been tuned on these patients.

    modal run --detach scripts/pseudobulk.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-pseudobulk")
image = (modal_app.image
         .add_local_python_source("modal_app")
         .add_local_file("data/hallmark.gmt", "/root/hallmark.gmt"))

IN = "/data/annotated_clean/as_annotated.h5ad"
OUT = "/data/pseudobulk_clean"

# A pseudobulk profile built from a handful of cells is noise, not a measurement.
MIN_CELLS = 20

HALLMARK_SETS = ["HALLMARK_INTERFERON_ALPHA_RESPONSE",
                 "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                 "HALLMARK_IL6_JAK_STAT3_SIGNALING"]

# No Hallmark set exists for Th17. This is hand-curated from the canonical
# markers and is therefore a weaker, more arguable definition than the above.
TH17_GENES = ["RORC", "RORA", "CCR6", "IL23R", "IL17A", "IL17F",
              "KLRB1", "CCL20", "IL1R1", "IL22", "AHR"]


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=8.0, timeout=7200)
def build():
    import os

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    os.makedirs(OUT, exist_ok=True)
    sc.settings.n_jobs = 8
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 300)

    adata = ad.read_h5ad(IN)
    print(f"loaded {adata.n_obs:,} x {adata.n_vars:,}")

    obs = adata.obs
    obs["sample_id"] = (obs["Patient"].astype(str) + "|" +
                        obs["Drug"].astype(str) + "|" + obs["Stage"].astype(str))
    obs["ct"] = obs["ct_fine"].astype(str)

    meta = (obs.groupby("sample_id", observed=True)
               .agg(Patient=("Patient", "first"), Drug=("Drug", "first"),
                    Stage=("Stage", "first"),
                    TNFi_Respond=("TNFi_Respond", "first"),
                    n_cells_sample=("ct", "size")))
    print(f"\n{len(meta)} samples (Patient x Drug x Stage)")

    # ---- 1. composition ----
    comp = pd.crosstab(obs["sample_id"], obs["ct"], normalize="index")
    comp.to_parquet(f"{OUT}/composition.parquet")
    print(f"composition: {comp.shape[0]} samples x {comp.shape[1]} cell types")

    print("\n=== cells per sample x cell type: how many groups clear the floor? ===")
    counts = pd.crosstab(obs["sample_id"], obs["ct"])
    ok = counts >= MIN_CELLS
    per_ct = pd.DataFrame({
        "total_cells": counts.sum(),
        "samples_with_data": ok.sum(),
        "pct_samples": (100 * ok.mean()).round(1),
        "median_cells": counts[counts > 0].median().round(0),
    }).sort_values("total_cells", ascending=False)
    print(per_ct.to_string())

    usable = per_ct.index[per_ct["pct_samples"] >= 90].tolist()
    print(f"\ncell types present in >=90% of samples ({len(usable)}): {usable}")

    # ---- 2. pseudobulk: sum RAW counts within each (sample, cell type) ----
    # Done as a sparse matrix product: (groups x cells) @ (cells x genes).
    key = obs["sample_id"].astype(str) + "||" + obs["ct"]
    groups = pd.Index(sorted(key.unique()))
    gidx = groups.get_indexer(key)
    ind = sp.csr_matrix(
        (np.ones(adata.n_obs, dtype=np.float32), (gidx, np.arange(adata.n_obs))),
        shape=(len(groups), adata.n_obs))

    counts_X = adata.layers["counts"]
    pb = ind @ counts_X                      # summed raw counts per group
    n_per_group = np.asarray(ind.sum(axis=1)).ravel()
    print(f"\npseudobulk matrix: {pb.shape[0]:,} groups x {pb.shape[1]:,} genes")

    pb_obs = pd.DataFrame({
        "group": groups,
        "sample_id": [g.split("||")[0] for g in groups],
        "ct": [g.split("||")[1] for g in groups],
        "n_cells": n_per_group,
    }).set_index("group")
    pb_obs = pb_obs.join(meta, on="sample_id")

    pbad = ad.AnnData(X=sp.csr_matrix(pb), obs=pb_obs, var=adata.var[[]].copy())
    pbad.layers["counts"] = pbad.X.copy()

    # Normalise for library size: a group of 5,000 cells sums to far more than
    # a group of 200, and that difference is not biology.
    sc.pp.normalize_total(pbad, target_sum=1e6)   # counts per million
    sc.pp.log1p(pbad)

    keep = pbad.obs["n_cells"] >= MIN_CELLS
    print(f"groups passing n_cells >= {MIN_CELLS}: {int(keep.sum()):,} "
          f"of {len(keep):,} ({100 * keep.mean():.1f}%)")
    pbad = pbad[keep].copy()
    pbad.write_h5ad(f"{OUT}/pseudobulk.h5ad", compression="gzip")

    # ---- 3. per-cell signature scores ----
    gmt = {}
    with open("/root/hallmark.gmt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gmt[parts[0]] = parts[2:]

    sigs = {name: gmt[name] for name in HALLMARK_SETS}
    sigs["TH17_CURATED"] = TH17_GENES

    print("\n=== signature scoring (per cell) ===")
    for name, genes in sigs.items():
        present = [g for g in genes if g in adata.var_names]
        # score_genes compares the set's mean expression against a reference
        # set of genes matched on expression level, so abundant genes do not
        # automatically produce high scores.
        sc.tl.score_genes(adata, present, score_name=name, ctrl_size=50)
        print(f"  {name:<38} {len(present):>3}/{len(genes):>3} genes found, "
              f"score range {adata.obs[name].min():.2f} to {adata.obs[name].max():.2f}")

    score_cols = list(sigs.keys())
    scores = (adata.obs.groupby(["sample_id", "ct"], observed=True)[score_cols]
                       .mean().reset_index())
    scores = scores.merge(counts.stack().rename("n_cells").reset_index()
                          .rename(columns={"col_0": "ct"}),
                          on=["sample_id", "ct"], how="left")
    scores = scores.merge(meta.reset_index(), on="sample_id", how="left")
    scores.to_parquet(f"{OUT}/scores.parquet")
    print(f"\nscores table: {scores.shape[0]:,} rows (sample x cell type)")

    # also keep a whole-sample (all cells pooled) version
    scores_all = (adata.obs.groupby("sample_id", observed=True)[score_cols]
                           .mean().join(meta))
    scores_all.to_parquet(f"{OUT}/scores_wholesample.parquet")

    print("\n=== first look: IFN-alpha score at BASELINE, by response ===")
    base = scores_all[scores_all["Stage"] == "base"]
    for cohort, sub in [("B expanded (all baselines)", base),
                        ("A strict (TNF arm only)", base[base["Drug"] == "TNF"])]:
        g = (sub.groupby("TNFi_Respond", observed=True)
                ["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
                .agg(["size", "mean", "std"]).round(4))
        print(f"\n  -- {cohort} --")
        print(g.to_string())

    vol.commit()
    print(f"\nwrote composition / pseudobulk / scores to {OUT}")


@app.local_entrypoint()
def main():
    build.remote()
