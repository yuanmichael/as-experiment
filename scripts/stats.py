"""Patient-level statistics: does anything separate responders from non-responders?

Fixes the duplicate-baseline problem (KAS04/KAS05 each have a TNF and an IL17
baseline draw, so sample-level tests count them twice) by collapsing to exactly
one row per patient before any test is run.

Tests, at baseline only:
  - signature scores, per cell type
  - cell type composition
  - AIF1 (the paper's biomarker), per cell type
  - a genome-wide scan per cell type, with FDR

    modal run --detach scripts/stats.py
"""
import modal

import modal_app
from modal_app import vol

app = modal.App("as-stats")
image = modal_app.image.add_local_python_source("modal_app")

PB = "/data/pseudobulk_clean"
OUT = "/data/stats_clean"

SIGS = ["HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "TH17_CURATED"]

CORE_CT = ["Classical monocytes", "Tem/Temra cytotoxic T cells", "CD16+ NK cells",
           "Tcm/Naive helper T cells", "Naive B cells", "Tem/Trm cytotoxic T cells",
           "Memory B cells", "Non-classical monocytes", "Tem/Effector helper T cells",
           "DC2"]


def _one_row_per_patient(df, cohort):
    """Collapse to a single baseline sample per patient.

    cohort 'A': TNF-arm baselines only (12 Good vs 4 Bad).
    cohort 'B': all baselines; for the two patients with both a TNF and an IL17
                baseline draw, keep the TNF one so the arms stay comparable.
    """
    d = df[df["Stage"] == "base"].copy()
    if cohort == "A":
        d = d[d["Drug"] == "TNF"]
    d["_pref"] = (d["Drug"] == "TNF").astype(int)        # prefer TNF draw
    d = (d.sort_values("_pref", ascending=False)
           .drop_duplicates(subset=["Patient"], keep="first")
           .drop(columns="_pref"))
    return d


@app.function(image=image, volumes={"/data": vol}, memory=32768, timeout=3600)
def run():
    import os

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    os.makedirs(OUT, exist_ok=True)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)

    scores = pd.read_parquet(f"{PB}/scores.parquet")
    comp = pd.read_parquet(f"{PB}/composition.parquet")
    pbad = ad.read_h5ad(f"{PB}/pseudobulk.h5ad")
    print(f"scores {scores.shape}, composition {comp.shape}, pseudobulk {pbad.shape}")

    # ---------- sanity: how many patients per cohort after de-duplication ----------
    for cohort in ["A", "B"]:
        d = _one_row_per_patient(scores.drop_duplicates("sample_id"), cohort)
        g = d.groupby("TNFi_Respond", observed=True)["Patient"].nunique()
        print(f"  cohort {cohort}: {dict(g)}  ({d['Patient'].nunique()} patients, "
              f"{len(d)} rows)")

    def mw(a, b):
        """Mann-Whitney U: no normality assumption, right for n of this size."""
        if len(a) < 3 or len(b) < 3:
            return np.nan, np.nan
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        # rank-biserial correlation as the effect size
        return p, 2 * u / (len(a) * len(b)) - 1

    # ---------- 1. signature scores per cell type ----------
    print("\n" + "=" * 78)
    print("1. SIGNATURE SCORES AT BASELINE, BY CELL TYPE")
    print("=" * 78)
    rows = []
    for cohort in ["A", "B"]:
        for ct in CORE_CT:
            sub = _one_row_per_patient(scores[scores["ct"] == ct], cohort)
            bad = sub[sub["TNFi_Respond"] == "Bad"]
            good = sub[sub["TNFi_Respond"] == "Good"]
            for sig in SIGS:
                p, eff = mw(bad[sig].dropna(), good[sig].dropna())
                rows.append({"cohort": cohort, "ct": ct, "signature": sig,
                             "n_bad": len(bad), "n_good": len(good),
                             "mean_bad": bad[sig].mean(), "mean_good": good[sig].mean(),
                             "delta": bad[sig].mean() - good[sig].mean(),
                             "p": p, "effect": eff})
    sig_res = pd.DataFrame(rows)
    for cohort in ["A", "B"]:
        c = sig_res[(sig_res["cohort"] == cohort) & sig_res["p"].notna()].copy()
        c["q"] = multipletests(c["p"], method="fdr_bh")[1]
        print(f"\n--- cohort {cohort}: top 15 by p-value ---")
        print(c.nsmallest(15, "p")[["ct", "signature", "n_bad", "n_good",
                                    "mean_bad", "mean_good", "delta", "p", "q"]]
              .round(4).to_string(index=False))
        print(f"  significant at q<0.05: {int((c['q'] < 0.05).sum())} of {len(c)}")

    print("\n--- IFN-alpha specifically, every cell type, cohort B ---")
    ifn = sig_res[(sig_res["cohort"] == "B") &
                  (sig_res["signature"] == "HALLMARK_INTERFERON_ALPHA_RESPONSE")]
    print(ifn[["ct", "n_bad", "n_good", "mean_bad", "mean_good", "delta", "p"]]
          .sort_values("p").round(4).to_string(index=False))
    sig_res.to_csv(f"{OUT}/signature_tests.csv", index=False)

    # ---------- 2. composition ----------
    print("\n" + "=" * 78)
    print("2. CELL TYPE COMPOSITION AT BASELINE")
    print("=" * 78)
    meta = scores.drop_duplicates("sample_id")[
        ["sample_id", "Patient", "Drug", "Stage", "TNFi_Respond"]]
    cm = comp.reset_index().merge(meta, on="sample_id")
    rows = []
    for cohort in ["A", "B"]:
        sub = _one_row_per_patient(cm, cohort)
        bad = sub[sub["TNFi_Respond"] == "Bad"]
        good = sub[sub["TNFi_Respond"] == "Good"]
        for ct in comp.columns:
            p, eff = mw(bad[ct], good[ct])
            rows.append({"cohort": cohort, "ct": ct,
                         "mean_bad": bad[ct].mean(), "mean_good": good[ct].mean(),
                         "delta": bad[ct].mean() - good[ct].mean(),
                         "p": p, "effect": eff})
    comp_res = pd.DataFrame(rows)
    for cohort in ["A", "B"]:
        c = comp_res[(comp_res["cohort"] == cohort) & comp_res["p"].notna()].copy()
        c["q"] = multipletests(c["p"], method="fdr_bh")[1]
        print(f"\n--- cohort {cohort}: top 10 ---")
        print(c.nsmallest(10, "p").round(4).to_string(index=False))
        print(f"  significant at q<0.05: {int((c['q'] < 0.05).sum())} of {len(c)}")
    comp_res.to_csv(f"{OUT}/composition_tests.csv", index=False)

    # ---------- 3. AIF1, the paper's biomarker ----------
    print("\n" + "=" * 78)
    print("3. AIF1 PSEUDOBULK EXPRESSION AT BASELINE")
    print("=" * 78)
    pb = pbad.obs.copy()
    pb["AIF1"] = np.asarray(pbad[:, "AIF1"].X.todense()).ravel()
    rows = []
    for cohort in ["A", "B"]:
        for ct in CORE_CT:
            sub = _one_row_per_patient(pb[pb["ct"] == ct], cohort)
            bad = sub[sub["TNFi_Respond"] == "Bad"]["AIF1"]
            good = sub[sub["TNFi_Respond"] == "Good"]["AIF1"]
            p, eff = mw(bad, good)
            rows.append({"cohort": cohort, "ct": ct, "n_bad": len(bad),
                         "n_good": len(good), "mean_bad": bad.mean(),
                         "mean_good": good.mean(),
                         "delta": bad.mean() - good.mean(), "p": p})
    aif1 = pd.DataFrame(rows)
    print(aif1.sort_values(["cohort", "p"]).round(4).to_string(index=False))
    aif1.to_csv(f"{OUT}/aif1_tests.csv", index=False)

    # ---------- 4. genome-wide scan ----------
    print("\n" + "=" * 78)
    print("4. GENOME-WIDE SCAN PER CELL TYPE (cohort B, baseline)")
    print("=" * 78)
    X = pbad.X.toarray() if hasattr(pbad.X, "toarray") else np.asarray(pbad.X)
    genes = np.array(pbad.var_names)
    allres = []
    for ct in CORE_CT:
        mask = (pbad.obs["ct"] == ct).values
        sub_obs = pbad.obs[mask].copy()
        sub_obs["_i"] = np.arange(mask.sum())
        sub_obs = _one_row_per_patient(sub_obs, "B")
        idx = sub_obs["_i"].values
        Xi = X[mask][idx]
        y = (sub_obs["TNFi_Respond"] == "Bad").values
        if y.sum() < 3 or (~y).sum() < 3:
            continue
        # drop genes that are essentially silent in this cell type
        expressed = (Xi > 0).mean(axis=0) > 0.5
        u, p = stats.mannwhitneyu(Xi[y][:, expressed], Xi[~y][:, expressed],
                                  alternative="two-sided", axis=0)
        q = multipletests(p, method="fdr_bh")[1]
        d = pd.DataFrame({"ct": ct, "gene": genes[expressed], "p": p, "q": q,
                          "mean_bad": Xi[y][:, expressed].mean(axis=0),
                          "mean_good": Xi[~y][:, expressed].mean(axis=0)})
        d["delta"] = d["mean_bad"] - d["mean_good"]
        allres.append(d)
        print(f"  {ct:<30} {int(expressed.sum()):>6,} genes tested, "
              f"n={int(y.sum())} bad vs {int((~y).sum())} good, "
              f"q<0.05: {int((q < 0.05).sum())}")
    gw = pd.concat(allres)
    gw.to_csv(f"{OUT}/genomewide.csv.gz", index=False, compression="gzip")
    print(f"\n  total gene x cell-type tests: {len(gw):,}")
    print(f"  passing q<0.05 anywhere:       {int((gw['q'] < 0.05).sum()):,}")
    print("\n  top 25 by raw p-value:")
    print(gw.nsmallest(25, "p").round(4).to_string(index=False))

    vol.commit()
    print(f"\nwrote results to {OUT}")


@app.local_entrypoint()
def main():
    run.remote()
