"""QC report for GSE277117 (run on Modal).

Deliberately does NOT filter anything. It measures the three standard quality
axes per cell, shows how they are distributed, and shows what each candidate
threshold would cost -- so the thresholds are a decision, not a default.

    modal run scripts/qc.py
"""
import modal

import modal_app
from modal_app import H5AD, vol

app = modal.App("as-qc")

# Modal only ships the entrypoint file into the container by default, so the
# sibling module has to be added to the image explicitly.
image = modal_app.image.add_local_python_source("modal_app")

OUT = "/data/qc"


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=4.0, timeout=3600)
def qc_report():
    import os

    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc

    os.makedirs(OUT, exist_ok=True)
    pd.set_option("display.width", 250)

    print("loading (this is the 6 GB step)...")
    adata = ad.read_h5ad(H5AD)
    adata.var_names_make_unique()
    print(f"  loaded: {adata.n_obs:,} x {adata.n_vars:,}")

    # Drop the droplets that caught two cells. Not a judgement call -- these
    # are not single cells, so no per-cell measurement of them is meaningful.
    n0 = adata.n_obs
    adata = adata[adata.obs["Patient"] != "Multiplet"].copy()
    print(f"  dropped Multiplet: {n0:,} -> {adata.n_obs:,} ({n0 - adata.n_obs:,} removed)")

    # The three standard QC axes.
    #   mt-  : mitochondrial genes. High fraction => the cell was dying/lysed
    #          and its cytoplasmic mRNA leaked out, leaving mitochondria behind.
    #   ribo : ribosomal protein genes. Very high fraction crowds out signal.
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]")
    adata.var["hb"] = adata.var_names.str.match(r"^HB[ABDEGQZ]\d?$")
    print(f"  mito genes: {adata.var['mt'].sum()}, "
          f"ribo: {adata.var['ribo'].sum()}, hemoglobin: {adata.var['hb'].sum()}")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=None,
        log1p=False, inplace=True,
    )

    obs = adata.obs
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]

    print("\n=== distribution across all cells ===")
    print(obs[metrics].describe(percentiles=[.01, .05, .25, .5, .75, .95, .99]).to_string())

    print("\n=== per-sample medians (looking for a bad library) ===")
    per = (obs.groupby("Sample_num", observed=True)
              .agg(cells=("total_counts", "size"),
                   med_genes=("n_genes_by_counts", "median"),
                   med_counts=("total_counts", "median"),
                   med_pct_mt=("pct_counts_mt", "median"),
                   p95_pct_mt=("pct_counts_mt", lambda s: np.percentile(s, 95)))
              .round(1))
    print(per.to_string())

    # What would each threshold cost us?
    print("\n=== cost of each candidate filter (applied alone) ===")
    cands = {
        "n_genes < 200":      obs["n_genes_by_counts"] < 200,
        "n_genes < 500":      obs["n_genes_by_counts"] < 500,
        "n_genes > 5000":     obs["n_genes_by_counts"] > 5000,
        "n_genes > 6000":     obs["n_genes_by_counts"] > 6000,
        "total_counts < 1000": obs["total_counts"] < 1000,
        "pct_mt > 5":         obs["pct_counts_mt"] > 5,
        "pct_mt > 10":        obs["pct_counts_mt"] > 10,
        "pct_mt > 15":        obs["pct_counts_mt"] > 15,
        "pct_mt > 20":        obs["pct_counts_mt"] > 20,
    }
    for name, mask in cands.items():
        n = int(mask.sum())
        print(f"  {name:<22} removes {n:>7,} cells ({100 * n / adata.n_obs:5.2f}%)")

    print("\n=== combined candidate filter sets ===")
    sets = {
        "lenient : genes>=200, mt<=20": (obs["n_genes_by_counts"] >= 200) & (obs["pct_counts_mt"] <= 20),
        "standard: genes 200-6000, mt<=15": (obs["n_genes_by_counts"].between(200, 6000)) & (obs["pct_counts_mt"] <= 15),
        "strict  : genes 500-5000, mt<=10": (obs["n_genes_by_counts"].between(500, 5000)) & (obs["pct_counts_mt"] <= 10),
    }
    for name, keep in sets.items():
        n = int(keep.sum())
        # does this filter hit some patients much harder than others?
        frac = obs.assign(k=keep).groupby("Patient", observed=True)["k"].mean()
        print(f"  {name:<34} keeps {n:>7,} ({100 * n / adata.n_obs:5.2f}%)  "
              f"worst-hit patient keeps {frac.min():.1%} ({frac.idxmin()})")

    # Figure: the classic QC triptych
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    axes[0].hist(obs["n_genes_by_counts"], bins=120)
    axes[0].set_xlabel("genes detected per cell"); axes[0].set_ylabel("cells")
    axes[0].axvline(200, color="r", ls="--"); axes[0].axvline(6000, color="r", ls="--")

    axes[1].hist(obs["pct_counts_mt"], bins=120)
    axes[1].set_xlabel("% mitochondrial counts"); axes[1].set_yscale("log")
    for t, c in [(10, "orange"), (15, "r")]:
        axes[1].axvline(t, color=c, ls="--")

    axes[2].scatter(obs["total_counts"], obs["n_genes_by_counts"],
                    c=obs["pct_counts_mt"], s=0.3, cmap="viridis", vmax=25)
    axes[2].set_xscale("log"); axes[2].set_yscale("log")
    axes[2].set_xlabel("total UMI counts"); axes[2].set_ylabel("genes detected")
    fig.colorbar(axes[2].collections[0], ax=axes[2], label="% mito")
    fig.tight_layout()
    fig.savefig(f"{OUT}/qc_overview.png", dpi=130)

    # Persist the metrics so we never recompute them
    obs[["Patient", "Drug", "Stage", "TNFi_Respond", "Sample_num"] + metrics
        + ["pct_counts_hb"]].to_parquet(f"{OUT}/qc_metrics.parquet")
    vol.commit()
    print(f"\nwrote {OUT}/qc_overview.png and {OUT}/qc_metrics.parquet")

    with open(f"{OUT}/qc_overview.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = qc_report.remote()
    out = "figures/qc_overview.png"
    with open(out, "wb") as f:
        f.write(png)
    print(f"saved {out} ({len(png) / 1e6:.1f} MB)")
