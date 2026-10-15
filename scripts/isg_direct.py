"""Test the paper's Fig 2e / Fig 3a claim WITHOUT clustering or batch correction.

Their finding: the proportion of ISG-high non-classical monocytes at baseline is
higher in non-responders. Our sub-clustering did not reproduce it, but the
sub-clustering ran Harmony on Sample_num -- which removes between-patient
variation, i.e. the signal itself.

Here, instead: score each monocyte for the Hallmark IFN-alpha program, call a
cell ISG-high by threshold, and compute the per-patient fraction. No clustering,
no integration, nothing that can absorb a between-patient effect.

Also verifies that layers["counts"] really is raw integer counts.

    modal run --detach scripts/isg_direct.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-isg")
image = (modal_app.image
         .add_local_python_source("modal_app")
         .add_local_file("data/hallmark.gmt", "/root/hallmark.gmt"))

IN = "/data/annotated_clean/as_annotated.h5ad"
OUT = "/data/isg"

MONO = ["Classical monocytes", "Non-classical monocytes"]


def _one_per_patient(df):
    d = df[df["Stage"] == "base"].copy()
    d["_p"] = (d["Drug"] == "TNF").astype(int)
    return (d.sort_values("_p", ascending=False)
             .drop_duplicates("Patient", keep="first").drop(columns="_p"))


@app.function(image=image, volumes={"/data": vol}, memory=65536,
              cpu=8.0, timeout=7200)
def run():
    import os

    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import stats

    os.makedirs(OUT, exist_ok=True)
    sc.settings.n_jobs = 8
    pd.set_option("display.width", 260)
    pd.set_option("display.max_rows", 200)

    full = ad.read_h5ad(IN)

    # ---- sanity check the counts layer ----
    c = full.layers["counts"][:500].toarray()
    x = full.X[:500].toarray()
    print("=== layer sanity ===")
    print(f"  counts: min={c.min():.3f} max={c.max():.3f} "
          f"integers={np.allclose(c, np.round(c))} median_rowsum={np.median(c.sum(1)):,.0f}")
    print(f"  X     : min={x.min():.3f} max={x.max():.3f} "
          f"integers={np.allclose(x, np.round(x))} median_rowsum={np.median(x.sum(1)):,.1f}")
    assert np.allclose(c, np.round(c)), "counts layer is NOT raw counts"

    sub = full[full.obs["ct_fine"].isin(MONO)].copy()
    sub.obs["sample_id"] = (sub.obs["Patient"].astype(str) + "|" +
                            sub.obs["Drug"].astype(str) + "|" +
                            sub.obs["Stage"].astype(str))
    print(f"\nmonocytes: {sub.n_obs:,}")

    # normalise fresh from raw counts
    sub.X = sub.layers["counts"].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    gmt = {}
    with open("/root/hallmark.gmt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            gmt[p[0]] = p[2:]
    ifn = [g for g in gmt["HALLMARK_INTERFERON_ALPHA_RESPONSE"] if g in sub.var_names]
    sc.tl.score_genes(sub, ifn, score_name="IFN", ctrl_size=50)
    print(f"IFN score: {len(ifn)} genes, range "
          f"{sub.obs['IFN'].min():.3f} to {sub.obs['IFN'].max():.3f}")

    q = sub.obs["IFN"].quantile([.5, .75, .9, .95, .99]).round(3)
    print(f"score quantiles: {dict(q)}")

    meta = (sub.obs.groupby("sample_id", observed=True)
                   .agg(Patient=("Patient", "first"), Drug=("Drug", "first"),
                        Stage=("Stage", "first"),
                        TNFi_Respond=("TNFi_Respond", "first")))

    print("\n" + "=" * 74)
    print("ISG-HIGH FRACTION AT BASELINE, BY THRESHOLD AND COMPARTMENT")
    print("=" * 74)
    rows = []
    for pct in [75, 90, 95, 99]:
        thr = sub.obs["IFN"].quantile(pct / 100)
        sub.obs["hi"] = sub.obs["IFN"] > thr
        for comp, mask in [
                ("all monocytes", np.ones(sub.n_obs, bool)),
                ("classical", (sub.obs["ct_fine"] == "Classical monocytes").values),
                ("non-classical", (sub.obs["ct_fine"] == "Non-classical monocytes").values)]:
            o = sub.obs[mask]
            frac = (o.groupby("sample_id", observed=True)["hi"].mean()
                     .rename("isg_hi").to_frame().join(meta))
            frac["n"] = o.groupby("sample_id", observed=True)["hi"].size()
            frac = frac[frac["n"] >= 20]
            for cohort in ["A", "B"]:
                d = _one_per_patient(frac.reset_index())
                if cohort == "A":
                    d = d[d["Drug"] == "TNF"]
                bad = d[d["TNFi_Respond"] == "Bad"]["isg_hi"]
                good = d[d["TNFi_Respond"] == "Good"]["isg_hi"]
                if len(bad) < 3 or len(good) < 3:
                    continue
                _, p = stats.mannwhitneyu(bad, good, alternative="two-sided")
                rows.append({"threshold": f"top {100 - pct}%", "compartment": comp,
                             "cohort": cohort, "n_bad": len(bad), "n_good": len(good),
                             "pct_bad": 100 * bad.mean(), "pct_good": 100 * good.mean(),
                             "p": p})
    res = pd.DataFrame(rows)
    print(res.round(4).to_string(index=False))
    res.to_csv(f"{OUT}/isg_threshold_tests.csv", index=False)

    # ---- also: plain mean IFN score within non-classical monocytes ----
    print("\n=== mean IFN score per patient, by compartment ===")
    for comp, mask in [
            ("classical", (sub.obs["ct_fine"] == "Classical monocytes").values),
            ("non-classical", (sub.obs["ct_fine"] == "Non-classical monocytes").values)]:
        o = sub.obs[mask]
        m = (o.groupby("sample_id", observed=True)["IFN"].mean()
              .rename("IFN").to_frame().join(meta))
        d = _one_per_patient(m.reset_index())
        bad = d[d["TNFi_Respond"] == "Bad"]["IFN"]
        good = d[d["TNFi_Respond"] == "Good"]["IFN"]
        _, p = stats.mannwhitneyu(bad, good, alternative="two-sided")
        print(f"  {comp:<16} NR={bad.mean():.4f}  R={good.mean():.4f}  p={p:.4f}")

    # ---- AIF1 within non-classical monocytes, and its link to IFN ----
    print("\n=== AIF1 vs IFN score, per cell (the paper's Fig 3b logic) ===")
    aif1 = np.asarray(sub[:, "AIF1"].X.todense()).ravel()
    sub.obs["AIF1"] = aif1
    for comp, mask in [("classical", (sub.obs["ct_fine"] == "Classical monocytes").values),
                       ("non-classical", (sub.obs["ct_fine"] == "Non-classical monocytes").values)]:
        o = sub.obs[mask]
        r, p = stats.spearmanr(o["AIF1"], o["IFN"])
        hi = o["IFN"] > o["IFN"].quantile(0.9)
        print(f"  {comp:<16} spearman(AIF1, IFN) r={r:.3f} p={p:.2e}   "
              f"AIF1 in ISG-high={o.loc[hi, 'AIF1'].mean():.3f} "
              f"vs rest={o.loc[~hi, 'AIF1'].mean():.3f}")

    # ---- figure ----
    thr = sub.obs["IFN"].quantile(0.90)
    sub.obs["hi"] = sub.obs["IFN"] > thr
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))

    ax = axes[0]
    ax.hist(sub.obs["IFN"], bins=140, color="#4C72B0")
    ax.axvline(thr, color="r", ls="--", label=f"90th pct = {thr:.3f}")
    ax.set_yscale("log"); ax.set_xlabel("per-cell IFN-alpha score")
    ax.set_ylabel("monocytes"); ax.legend(); ax.set_title("is there an ISG-high mode?")

    for ax, comp in zip(axes[1:], ["Classical monocytes", "Non-classical monocytes"]):
        o = sub.obs[sub.obs["ct_fine"] == comp]
        frac = (o.groupby("sample_id", observed=True)["hi"].mean()
                 .rename("isg_hi").to_frame().join(meta))
        d = _one_per_patient(frac.reset_index())
        groups = [100 * d[d["TNFi_Respond"] == g]["isg_hi"].values
                  for g in ["Good", "Bad"]]
        bp = ax.boxplot(groups, tick_labels=["Responder", "Non-responder"],
                        widths=0.55, patch_artist=True, showfliers=False)
        for patch, col in zip(bp["boxes"], ["#4C72B0", "#C44E52"]):
            patch.set_facecolor(col); patch.set_alpha(0.4)
        for k, (v, col) in enumerate(zip(groups, ["#4C72B0", "#C44E52"])):
            ax.scatter(np.random.normal(k + 1, 0.06, len(v)), v, color=col,
                       s=40, zorder=3, edgecolor="white", lw=0.7)
        _, p = stats.mannwhitneyu(groups[1], groups[0], alternative="two-sided")
        ax.set_title(f"{comp}\nISG-high fraction, p = {p:.3f}")
        ax.set_ylabel("% ISG-high"); ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(f"{OUT}/isg_direct.png", dpi=125)
    vol.commit()
    with open(f"{OUT}/isg_direct.png", "rb") as f:
        return f.read()


@app.local_entrypoint()
def main():
    png = run.remote()
    with open("figures/isg_direct.png", "wb") as f:
        f.write(png)
    print("saved figures/isg_direct.png")
