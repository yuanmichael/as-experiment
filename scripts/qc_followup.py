"""Follow-up on the QC report: is the quality variation tracking the outcome?

Reads only the saved per-cell metrics table (small), not the count matrix.

    modal run scripts/qc_followup.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-qc-followup")
image = modal_app.image.add_local_python_source("modal_app")


@app.function(image=image, volumes={"/data": vol}, memory=4096, timeout=900)
def followup():
    import numpy as np
    import pandas as pd

    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)

    df = pd.read_parquet("/data/qc/qc_metrics.parquet")

    print("=== sample -> patient / drug / stage / label ===")
    m = (df.groupby("Sample_num", observed=True)
           .agg(patients=("Patient", lambda s: ",".join(sorted(s.unique()))),
                drug=("Drug", lambda s: s.iloc[0]),
                stage=("Stage", lambda s: s.iloc[0]),
                label=("TNFi_Respond", lambda s: ",".join(sorted(s.unique()))),
                cells=("Patient", "size"),
                med_genes=("n_genes_by_counts", "median"),
                p95_mt=("pct_counts_mt", lambda s: round(np.percentile(s, 95), 1)),
                frac_mt_over_50=("pct_counts_mt", lambda s: round((s > 50).mean(), 4)),
                med_hb=("pct_counts_hb", "median"))
           .sort_values("p95_mt", ascending=False))
    print(m.to_string())

    print("\n=== is library size (pooled vs single) tied to the mito tail? ===")
    m["pooled"] = m["patients"].str.contains(",")
    print(m.groupby("pooled")[["cells", "med_genes", "p95_mt", "frac_mt_over_50"]]
           .median().round(3).to_string())

    print("\n=== THE CHECK: does cell quality differ by outcome group? ===")
    for cohort, sub in [
        ("B expanded (all baselines)", df[df["Stage"] == "base"]),
        ("A strict  (TNF arm baseline)",
         df[(df["Stage"] == "base") & (df["Drug"] == "TNF")]),
    ]:
        print(f"\n  -- {cohort} --")
        per_pt = (sub.groupby(["Patient", "TNFi_Respond"], observed=True)
                     .agg(cells=("n_genes_by_counts", "size"),
                          med_genes=("n_genes_by_counts", "median"),
                          med_mt=("pct_counts_mt", "median"),
                          frac_mt_over_20=("pct_counts_mt", lambda s: (s > 20).mean()))
                     .reset_index())
        print(per_pt.round(3).to_string(index=False))
        g = per_pt.groupby("TNFi_Respond", observed=True)[
            ["cells", "med_genes", "med_mt", "frac_mt_over_20"]].median().round(3)
        print("  group medians (per patient):")
        print(g.to_string())

    print("\n=== low-gene samples: platelet/RBC contamination? ===")
    low = m.nsmallest(5, "med_genes")[["patients", "med_genes", "med_hb", "cells"]]
    high = m.nlargest(5, "med_genes")[["patients", "med_genes", "med_hb", "cells"]]
    print("lowest med_genes:\n" + low.to_string())
    print("highest med_genes:\n" + high.to_string())
