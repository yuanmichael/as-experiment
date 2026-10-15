"""Redo the resolution sweep against INDEPENDENT labels.

The first sweep was circular: CellTypist's majority_voting used the res=1.0
Leiden clustering as its over-clustering, so `ct_fine` is constant within every
res=1.0 cluster by construction. Scoring that clustering against those labels
measures nothing.

`ct_fine_raw` holds the per-cell predictions made before majority voting. Those
are derived only from each cell's own expression, so they are independent of any
Leiden run and can referee the comparison fairly.

    modal run scripts/resolution_fix.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-resfix")
image = modal_app.image.add_local_python_source("modal_app")

IN = "/data/annotated/as_annotated.h5ad"
RESOLUTIONS = [0.2, 0.4, 0.6, 1.0, 1.5, 2.0, 3.0]


@app.function(image=image, volumes={"/data": vol}, memory=16384, timeout=1800)
def resweep():
    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    pd.set_option("display.width", 220)
    obs = ad.read_h5ad(IN, backed="r").obs

    print(f"ct_fine (majority-voted): {obs['ct_fine'].nunique()} labels")
    print(f"ct_fine_raw (per-cell):   {obs['ct_fine_raw'].nunique()} labels")
    print(f"median per-cell confidence: {obs['ct_fine_conf'].median():.3f}")
    print(f"cells below 0.5 confidence: "
          f"{(obs['ct_fine_conf'] < 0.5).mean() * 100:.2f}%")

    print("\n=== raw per-cell labels (independent of any clustering) ===")
    vc = obs["ct_fine_raw"].value_counts()
    print(pd.DataFrame({"cells": vc, "pct": (100 * vc / len(obs)).round(2)})
          .head(30).to_string())

    print("\n=== corrected sweep: clustering vs INDEPENDENT per-cell labels ===")
    lab = obs["ct_fine_raw"].astype(str)
    rows = []
    for res in RESOLUTIONS:
        key = f"leiden_r{res}"
        clu = obs[key].astype(str)
        ari = adjusted_rand_score(lab, clu)
        nmi = normalized_mutual_info_score(lab, clu)
        pur = (obs.groupby(key, observed=True)["ct_fine_raw"]
                  .apply(lambda s: s.value_counts(normalize=True).iloc[0]).median())
        rows.append({"resolution": res, "clusters": obs[key].nunique(),
                     "ARI": round(ari, 4), "NMI": round(nmi, 4),
                     "median_purity": round(pur, 4)})
        print(f"  res={res:<4} clusters={rows[-1]['clusters']:>3}  "
              f"ARI={ari:.4f}  NMI={nmi:.4f}  purity={pur:.3f}")
    sweep = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(sweep["resolution"], sweep["ARI"], "o-", label="ARI")
    ax.plot(sweep["resolution"], sweep["NMI"], "s-", label="NMI")
    ax.plot(sweep["resolution"], sweep["median_purity"], "^-", label="median purity")
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("agreement with independent per-cell labels")
    ax.set_title("corrected resolution sweep")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("/data/annotated/resolution_sweep_fixed.png", dpi=130)
    vol.commit()

    with open("/data/annotated/resolution_sweep_fixed.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = resweep.remote()
    with open("figures/resolution_sweep_fixed.png", "wb") as f:
        f.write(png)
    print("saved figures/resolution_sweep_fixed.png")
