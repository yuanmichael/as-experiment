"""Could GSE277117 have found the TRBV9 story? And can we see the CD8 phenotype?

Two separate questions:

  1. TCR V-gene usage. TRBV9 clonotypes are the AS finding (Faham: TRBV9-SVGLYST
     in 30-39% of cases, 0% of controls). Reading V-gene usage needs either a
     dedicated VDJ library or 5' GEX chemistry. With 3' chemistry the reads sit
     at the 3' end of the transcript -- over the CONSTANT region (TRBC1/TRBC2),
     not the V segment 1-2 kb upstream. So this is a chemistry question, and it
     is answerable from the data: if TRBC2 is abundant while every TRBV gene is
     near zero, the library cannot see V segments.

  2. The CD8 phenotype (KLRB1/CD161+, CCR6+, exhaustion-resistant) is ordinary
     gene expression and IS readable from 3' data. Testable here.

     Confound to handle: MAIT cells are KLRB1-high and CCR6+ by definition, so
     a naive KLRB1+CCR6+ gate just re-finds MAIT cells. They are separated out
     via SLC4A10 / TRAV1-2 and reported both ways.

    modal run --detach scripts/trbv_check.py::tcr_feasibility
    modal run --detach scripts/trbv_check.py::cd8_phenotype
"""
import modal

H5AD = "/data/annotated_clean/as_annotated.h5ad"

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("scanpy==1.12.4", "anndata", "h5py", "pandas", "scikit-learn",
                 "pyarrow", "numpy", "scipy", "statsmodels")
)

EXHAUSTION = ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1",
              "BTLA", "CD244", "VSIR"]
CYTOTOXIC = ["GZMB", "GZMA", "GZMK", "PRF1", "GNLY", "NKG7", "IFNG", "TNF"]
TC17 = ["KLRB1", "CCR6", "RORC", "IL23R", "IL17A", "CXCR6", "ZBTB16"]
MAIT = ["SLC4A10", "TRAV1-2", "KLRB1", "ZBTB16", "NCR3"]

CD8_TYPES = ["Tem/Temra cytotoxic T cells", "Tem/Trm cytotoxic T cells",
             "Tcm/Naive cytotoxic T cells", "MAIT cells"]


@app.function(image=image, volumes={"/data": vol}, memory=32768, timeout=3600)
def tcr_feasibility():
    """Is V-gene usage readable at all from this library?"""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    pd.set_option("display.width", 200)
    a = ad.read_h5ad(H5AD)
    X = a.layers["counts"] if "counts" in a.layers else a.X
    X = sp.csr_matrix(X).tocsc()
    genes = np.asarray(a.var_names)

    def stats(name):
        if name not in set(genes):
            return None
        j = int(np.where(genes == name)[0][0])
        col = X[:, j]
        return dict(gene=name, total_counts=float(col.sum()),
                    pct_cells=100.0 * col.nnz / X.shape[0],
                    mean_per_cell=float(col.sum() / X.shape[0]))

    print("=" * 88)
    print("Q1  Is TRBV9 even in the feature list, and is it measurable?")
    print("=" * 88)

    trbv = sorted(g for g in genes if g.startswith("TRBV"))
    trav = sorted(g for g in genes if g.startswith("TRAV"))
    const = [g for g in ("TRBC1", "TRBC2", "TRAC", "TRDC", "TRGC1", "TRGC2")
             if g in set(genes)]
    print(f"TRBV genes in matrix: {len(trbv)}")
    print(f"TRAV genes in matrix: {len(trav)}")
    print(f"constant-region genes present: {const}")

    rows = [stats(g) for g in const]
    rows += [stats(g) for g in ["TRBV9"] if g in set(genes)]
    rows += [stats(g) for g in trbv[:12]]
    rows = [r for r in rows if r]
    df = pd.DataFrame(rows).sort_values("total_counts", ascending=False)
    print("\n" + df.to_string(index=False))

    if trbv:
        tot_v = sum(stats(g)["total_counts"] for g in trbv)
        tot_c = sum(stats(g)["total_counts"] for g in const) or 1
        print(f"\nTotal TRBV counts across all {len(trbv)} V genes: {tot_v:,.0f}")
        print(f"Total TCR constant-region counts:                {tot_c:,.0f}")
        print(f"V : C ratio = {tot_v / tot_c:.4f}")
        print("\n  A 5' library gives V:C on the order of 0.1-1.")
        print("  A 3' library gives V:C near zero -- V segments are ~1-2 kb")
        print("  upstream of where the reads land, so they are invisible.")

    print("\n" + "=" * 88)
    print("Q2  Does the object carry a VDJ/clonotype annotation?")
    print("=" * 88)
    cand = [c for c in a.obs.columns
            if any(k in c.lower() for k in
                   ("tcr", "clono", "vdj", "cdr3", "trb", "tra_", "chain"))]
    print(f"obs columns suggesting VDJ data: {cand or 'NONE'}")
    print(f"uns keys: {list(a.uns.keys())}")


@app.function(image=image, volumes={"/data": vol}, memory=65536, timeout=7200)
def cd8_phenotype():
    """The KLRB1+/CCR6+ exhaustion-resistant CD8 population, at patient level."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from scipy import stats as st

    pd.set_option("display.width", 220)
    a = ad.read_h5ad(H5AD)
    print(f"{a.n_obs:,} cells")

    # normalise from raw counts; the stored .X may already be processed
    X = a.layers["counts"] if "counts" in a.layers else a.X
    X = sp.csr_matrix(X)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    X = sp.diags(1e4 / tot) @ X
    X.data = np.log1p(X.data)
    genes = pd.Index(np.asarray(a.var_names))

    obs = a.obs
    is_cd8 = obs["ct_fine"].astype(str).isin(CD8_TYPES).to_numpy()
    print(f"CD8/cytotoxic-lineage cells: {is_cd8.sum():,}")
    print(obs.loc[is_cd8, "ct_fine"].value_counts().to_string())

    def vec(g):
        if g not in genes:
            return None
        return np.asarray(X[:, genes.get_loc(g)].todense()).ravel()

    def score(gs):
        have = [g for g in gs if g in genes]
        if not have:
            return None, []
        m = np.vstack([vec(g) for g in have])
        z = (m - m.mean(axis=1, keepdims=True)) / (m.std(axis=1, keepdims=True) + 1e-9)
        return z.mean(axis=0), have

    print("\n" + "=" * 88)
    print("Marker availability")
    print("=" * 88)
    for nm, gs in (("Tc17/CD161", TC17), ("exhaustion", EXHAUSTION),
                   ("cytotoxic", CYTOTOXIC), ("MAIT", MAIT)):
        print(f"  {nm:12s} {[g for g in gs if g in genes]}")
        missing = [g for g in gs if g not in genes]
        if missing:
            print(f"               missing: {missing}")

    klrb1, ccr6 = vec("KLRB1"), vec("CCR6")
    if klrb1 is None or ccr6 is None:
        print("KLRB1 or CCR6 absent; cannot proceed")
        return

    # Gate on detection, not an arbitrary quantile -- these are sparse in 3' data
    gate = is_cd8 & (klrb1 > 0) & (ccr6 > 0)
    print(f"\nKLRB1+ CCR6+ CD8 cells: {gate.sum():,} "
          f"({100 * gate.sum() / max(is_cd8.sum(), 1):.1f}% of CD8 lineage)")

    is_mait = obs["ct_fine"].astype(str).eq("MAIT cells").to_numpy()
    slc = vec("SLC4A10")
    mait_like = is_mait | (slc is not None and (slc > 0))
    print(f"  of which MAIT or SLC4A10+: {int((gate & mait_like).sum()):,}")
    print(f"  non-MAIT KLRB1+CCR6+ CD8:  {int((gate & ~mait_like).sum()):,}")

    exh, exh_g = score(EXHAUSTION)
    cyt, _ = score(CYTOTOXIC)

    print("\n" + "=" * 88)
    print("Are KLRB1+CCR6+ CD8 cells less exhausted than other CD8 cells?")
    print("=" * 88)
    other = is_cd8 & ~gate
    for label, mask in (("KLRB1+CCR6+ (all)", gate),
                        ("KLRB1+CCR6+ (non-MAIT)", gate & ~mait_like),
                        ("other CD8", other)):
        n = int(mask.sum())
        if n < 50:
            print(f"  {label:26s} n={n:,} -- too few")
            continue
        print(f"  {label:26s} n={n:7,}  exhaustion {exh[mask].mean():+.3f}  "
              f"cytotoxic {cyt[mask].mean():+.3f}")
    u, p = st.mannwhitneyu(exh[gate & ~mait_like], exh[other])
    print(f"\n  non-MAIT KLRB1+CCR6+ vs other CD8, exhaustion score: p = {p:.3g}")
    print("  (cell-level test -- NOT a patient-level claim; cells are not independent)")

    # Patient level: the only unit that can carry a response label
    print("\n" + "=" * 88)
    print("Patient-level: does this population track TNFi response?")
    print("=" * 88)
    base = obs["Stage"].astype(str).eq("base").to_numpy()
    valid = base & (obs["Patient"].astype(str) != "Multiplet").to_numpy()
    rows = []
    pats = pd.unique(obs.loc[valid, "Patient"].astype(str))
    for p_ in pats:
        m = valid & (obs["Patient"].astype(str) == p_).to_numpy()
        cd8 = m & is_cd8
        if cd8.sum() < 50:
            continue
        g = m & gate & ~mait_like
        rows.append(dict(
            patient=p_,
            drug=obs.loc[m, "Drug"].astype(str).iloc[0],
            resp=obs.loc[m, "TNFi_Respond"].astype(str).iloc[0],
            n_cd8=int(cd8.sum()),
            pct_klrb1_ccr6=100.0 * g.sum() / cd8.sum(),
            exh_in_cd8=float(exh[cd8].mean()),
            exh_in_gate=float(exh[g].mean()) if g.sum() >= 10 else np.nan,
        ))
    df = pd.DataFrame(rows)
    print(df.round(3).to_string(index=False))

    for cohort, sub in (("B (all baselines)", df),
                        ("A (TNF arm only)", df[df["drug"] == "TNF"])):
        good = sub[sub["resp"] == "Good"]["pct_klrb1_ccr6"]
        bad = sub[sub["resp"] == "Bad"]["pct_klrb1_ccr6"]
        if len(good) < 2 or len(bad) < 2:
            print(f"\n  {cohort}: too few patients ({len(good)} vs {len(bad)})")
            continue
        u, p = st.mannwhitneyu(bad, good)
        print(f"\n  {cohort}: responders n={len(good)} median {good.median():.2f}%  |  "
              f"non-responders n={len(bad)} median {bad.median():.2f}%  |  p = {p:.3f}")


@app.local_entrypoint()
def main():
    tcr_feasibility.remote()
    cd8_phenotype.remote()
