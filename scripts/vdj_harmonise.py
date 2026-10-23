"""Put 291 libraries on one gene axis, and pick highly variable genes ACROSS studies.

Each library was ingested keeping its own gene axis on purpose -- 32,285 to
36,601 genes depending on the reference the submitters used -- so that a bad
intersection could not silently destroy data at ingest. This is where they are
reconciled.

Two things have to be right.

1. Duplicate gene symbols. A 10x features.tsv maps several Ensembl IDs to the
   same symbol, so a naive lookup silently keeps whichever column came first and
   throws away real counts. Columns sharing a symbol are summed.

2. HVG selection has to be done PER STUDY and then combined, never pooled.
   Pooling makes the top variable genes the ones that differ most BETWEEN
   studies -- which is batch, not biology. Selecting within each study and
   keeping genes that recur ranks on reproducible variation instead. The
   overlap between studies' own HVG lists is reported, because if that overlap
   is poor the batch problem is worse than the gene axis problem and no
   intersection fixes it.

Reads what the build already wrote to the volume. Downloads nothing.

    modal run scripts/vdj_harmonise.py::stats_all
    modal run scripts/vdj_harmonise.py::harmonise
"""

import json

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
VDJ = "/data/vdj"
HARM = "/data/vdj/_harmonised"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("pandas", "numpy", "scipy", "pyarrow")
)

# Dropped after measurement: mostly not T cells (1.1% and 4.6% median pairing).
DROP = {"GSE256141", "GSE292177"}
# Never contributes to gene selection either -- choosing genes using the test
# set is leakage, even unsupervised leakage.
TEST = {"GSE288581"}

N_HVG = 2000
TOP_PER_STUDY = 4000      # candidate pool within each study before combining

# Minimum fraction of cells in which a gene must be detected to be eligible.
#
# Without this the selection is pathological. var/mean rewards a gene detected
# in 0.3% of cells that happens to spike in a few: the first run returned 2,000
# HVGs whose MEDIAN detection rate was 0.4%, so about half the feature budget
# went to genes absent from 99.6% of cells. Binning dispersion by mean
# expression is supposed to control this and does not control it enough.
MIN_DETECT_PCT = 1.0

# Genes guaranteed a place in the output space regardless of dispersion rank.
#
# These are not used to DEFINE the label -- the target stays an unsupervised
# latent, so no hand-written signature caps the ceiling. They are here because
# the latent has to be interpretable afterwards: if a programme's marker genes
# are not in the output space, that programme cannot be read out at all.
# HAVCR2 (TIM-3), detected in 2.4% of cells, missed the dispersion cut on the
# first pass and would have taken the exhaustion readout with it.
MARKER_PANEL = [
    # exhaustion / inhibitory
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "ENTPD1", "CD160",
    "BTLA", "EOMES",
    # cytotoxicity
    "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7", "IFNG", "KLRG1",
    "FGFBP2",
    # naive / memory / stemness
    "TCF7", "SELL", "CCR7", "LEF1", "IL7R", "CD27", "CD28", "BACH2", "FAS",
    # regulatory
    "FOXP3", "IL2RA", "IKZF2", "TNFRSF18",
    # lineage and unconventional
    "CD3E", "CD3D", "CD4", "CD8A", "CD8B", "TRAC", "TRBC2",
    "KLRB1", "SLC4A10", "ZBTB16", "TRDC", "TRGC1",
    # activation / proliferation / stress
    "MKI67", "HLA-DRA", "IL2", "TNF", "NR4A1", "HSPA1A",
]


@app.function(image=image, volumes={"/data": vol}, cpu=2.0, memory=8192,
              timeout=1800, max_containers=24, retries=2)
def stats_one(lib: str) -> dict:
    """Per-gene detection and log-normalised moments for one library."""
    import os

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    gse = lib.split("/")[0]
    d = f"{VDJ}/{lib}"
    out = f"{d}/gene_stats.parquet"
    if os.path.exists(out):
        return {"lib": lib, "gse": gse, "status": "cached"}

    with open(f"{d}/genes.txt") as fh:
        genes = [g.strip() for g in fh.read().splitlines()]
    X = sp.load_npz(f"{d}/expr.npz").tocsc()
    if X.shape[1] != len(genes):
        return {"lib": lib, "gse": gse, "status":
                f"shape mismatch {X.shape[1]} vs {len(genes)}"}

    # Sum columns that share a gene symbol, rather than keeping the first.
    sym = pd.Index(genes)
    n_dup = int(len(sym) - sym.nunique())
    if n_dup:
        codes, uniq = pd.factorize(sym)
        agg = sp.csr_matrix(
            (np.ones(len(codes)), (np.arange(len(codes)), codes)),
            shape=(len(codes), len(uniq)))
        X = (X @ agg).tocsc()
        genes = list(uniq)

    # CPM to 10k, log1p -- the standard scale for dispersion.
    lib_sz = np.asarray(X.sum(axis=1)).ravel()
    lib_sz[lib_sz == 0] = 1
    Xn = sp.diags(1e4 / lib_sz) @ X.tocsr()
    Xn = Xn.tocoo()
    Xn.data = np.log1p(Xn.data)
    Xn = Xn.tocsc()

    n_cells = X.shape[0]
    s1 = np.asarray(Xn.sum(axis=0)).ravel()
    s2 = np.asarray(Xn.multiply(Xn).sum(axis=0)).ravel()
    n_det = np.diff(X.tocsc().indptr)

    pd.DataFrame({"gene": genes, "n_det": n_det, "s1": s1, "s2": s2,
                  "n_cells": n_cells, "gse": gse, "lib": lib}).to_parquet(
        out, index=False)
    vol.commit()
    return {"lib": lib, "gse": gse, "status": "ok", "n_genes": len(genes),
            "n_cells": int(n_cells), "n_dup_symbols": n_dup}


@app.function(image=image, volumes={"/data": vol}, timeout=7200)
def stats_all():
    import glob
    import os

    libs = []
    for p in sorted(glob.glob(f"{VDJ}/*/*/genes.txt")):
        gse, name = p.split("/")[-3], p.split("/")[-2]
        if gse in DROP:
            continue
        libs.append(f"{gse}/{name}")
    print(f"{len(libs)} libraries ({len(DROP)} studies excluded)")

    res = list(stats_one.map(libs, return_exceptions=True))
    ok = [r for r in res if isinstance(r, dict) and r.get("status") in ("ok", "cached")]
    print(f"{len(ok)}/{len(res)} done")
    dups = sum(r.get("n_dup_symbols", 0) for r in res if isinstance(r, dict))
    print(f"duplicate gene symbols summed: {dups:,} columns across all libraries")
    for r in res:
        if not isinstance(r, dict):
            print(f"  EXCEPTION {r}")
        elif r.get("status") not in ("ok", "cached"):
            print(f"  FAILED {r.get('lib')}: {r.get('status')}")
    vol.reload()
    return len(ok)


@app.function(image=image, volumes={"/data": vol}, cpu=4.0, memory=32768,
              timeout=7200)
def harmonise():
    import glob
    import os

    import numpy as np
    import pandas as pd

    vol.reload()
    files = [p for p in sorted(glob.glob(f"{VDJ}/*/*/gene_stats.parquet"))
             if p.split("/")[-3] not in DROP]
    print(f"reading {len(files)} per-library stats")
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)

    # ---------------------------------------------------------- the gene axis
    per_lib = df.groupby("lib")["gene"].nunique()
    sets = df.groupby("lib")["gene"].apply(set)
    common = set.intersection(*sets.tolist())
    union = set.union(*sets.tolist())
    print("\n" + "=" * 88)
    print("GENE AXIS")
    print("=" * 88)
    print(f"  genes per library   min {per_lib.min():,}  median "
          f"{int(per_lib.median()):,}  max {per_lib.max():,}")
    print(f"  union               {len(union):,}")
    print(f"  intersection        {len(common):,}  "
          f"({100*len(common)/len(union):.1f}% of union)")

    by_gse = df.groupby("gse")["gene"].apply(set)
    for gse, s in by_gse.items():
        print(f"    {gse:12s} {len(s):6,} genes, "
              f"{100*len(s & common)/max(len(common),1):5.1f}% of the intersection")

    # ------------------------------------------------- per-study HVG, then vote
    #
    # Selecting HVGs on pooled data would rank genes by BETWEEN-study variance,
    # i.e. batch. Rank within each study, then keep genes that recur.
    train = df[~df["gse"].isin(TEST)]
    print(f"\n  gene selection uses {train['gse'].nunique()} training studies; "
          f"{', '.join(sorted(TEST))} excluded (selecting on the test set is "
          f"leakage)")

    votes, per_study_sets = {}, {}
    for gse, g in train.groupby("gse"):
        a = (g[g["gene"].isin(common)]
             .groupby("gene")
             .agg(s1=("s1", "sum"), s2=("s2", "sum"),
                  n=("n_cells", "sum"), det=("n_det", "sum")))
        a = a[a["n"] > 0]
        mean = a["s1"] / a["n"]
        var = (a["s2"] / a["n"] - mean ** 2).clip(lower=0)
        det_pct = 100 * a["det"] / a["n"]
        keep = (mean > 0) & (det_pct >= MIN_DETECT_PCT)
        print(f"    {gse:12s} {int(keep.sum()):6,} of {len(a):,} genes pass "
              f"the {MIN_DETECT_PCT}% detection floor")
        mean, var = mean[keep], var[keep]
        # Seurat-style: dispersion normalised within mean-expression bins, so
        # highly expressed genes do not dominate purely by having more variance.
        disp = var / mean
        bins = pd.qcut(mean, min(20, max(2, mean.nunique() // 50)),
                       duplicates="drop")
        z = disp.groupby(bins, observed=True).transform(
            lambda s: (s - s.mean()) / (s.std() if s.std() > 0 else 1.0))
        top = set(z.sort_values(ascending=False).head(TOP_PER_STUDY).index)
        per_study_sets[gse] = top
        for gene in top:
            votes[gene] = votes.get(gene, 0) + 1

    n_studies = len(per_study_sets)
    vote_s = pd.Series(votes).sort_values(ascending=False)
    print("\n" + "=" * 88)
    print(f"HVG AGREEMENT ACROSS {n_studies} STUDIES "
          f"(top {TOP_PER_STUDY:,} within each)")
    print("=" * 88)
    for k in range(n_studies, 0, -1):
        print(f"  in {k}/{n_studies} studies' top list: "
              f"{int((vote_s == k).sum()):6,} genes")
    pairwise = []
    keys = sorted(per_study_sets)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = per_study_sets[keys[i]], per_study_sets[keys[j]]
            pairwise.append(len(a & b) / TOP_PER_STUDY)
    if pairwise:
        print(f"\n  median pairwise overlap of study HVG lists: "
              f"{100*np.median(pairwise):.1f}%")
        print("  (low overlap would mean batch dominates biology, and no gene")
        print("   axis fixes that -- it would have to be handled in the model)")

    hvg = list(vote_s.head(N_HVG).index)
    forced = [g for g in MARKER_PANEL if g in common and g not in hvg]
    missing = [g for g in MARKER_PANEL if g not in common]
    hvg = hvg + forced
    print(f"\n  {len(forced)} marker-panel genes added that dispersion ranking "
          f"missed: {', '.join(forced) if forced else '(none)'}")
    if missing:
        print(f"  NOT ON THE COMMON GENE AXIS AT ALL: {', '.join(missing)}")

    # ------------------------------------------------------------------ write
    os.makedirs(HARM, exist_ok=True)
    genes_sorted = sorted(common)
    with open(f"{HARM}/common_genes.txt", "w") as fh:
        fh.write("\n".join(genes_sorted))
    with open(f"{HARM}/hvg_{N_HVG}.txt", "w") as fh:
        fh.write("\n".join(hvg))
    vote_s.rename("n_studies").to_frame().to_parquet(f"{HARM}/hvg_votes.parquet")

    agg = (df[df["gene"].isin(common)].groupby("gene")
             .agg(s1=("s1", "sum"), s2=("s2", "sum"), n=("n_cells", "sum"),
                  det=("n_det", "sum")))
    agg["mean"] = agg["s1"] / agg["n"]
    agg["var"] = (agg["s2"] / agg["n"] - agg["mean"] ** 2).clip(lower=0)
    agg["pct_cells_detected"] = 100 * agg["det"] / agg["n"]
    agg["is_hvg"] = agg.index.isin(hvg)
    agg.to_parquet(f"{HARM}/gene_stats.parquet")

    meta = dict(n_libraries=len(files), n_studies=int(df["gse"].nunique()),
                n_common_genes=len(common), n_union_genes=len(union),
                n_hvg=len(hvg), top_per_study=TOP_PER_STUDY,
                min_detect_pct=MIN_DETECT_PCT,
                median_hvg_detect_pct=round(
                    float(agg.loc[hvg, "pct_cells_detected"].median()), 2),
                selection_studies=sorted(per_study_sets),
                excluded_from_selection=sorted(TEST),
                dropped_studies=sorted(DROP),
                median_pairwise_hvg_overlap=round(float(np.median(pairwise)), 4)
                if pairwise else None)
    with open(f"{HARM}/meta.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    vol.commit()

    print("\n" + "=" * 88)
    print("SELECTED")
    print("=" * 88)
    h = agg[agg["is_hvg"]]
    print(f"  {len(hvg):,} HVGs from {len(common):,} common genes")
    print(f"  median % of cells detecting an HVG: "
          f"{h['pct_cells_detected'].median():.1f}%  "
          f"(all genes: {agg['pct_cells_detected'].median():.1f}%)")
    print(f"\n  sanity -- are known T cell programme genes in?")
    for g in ["CD3E", "CD8A", "CD4", "GZMB", "PRF1", "IFNG", "PDCD1", "TOX",
              "LAG3", "HAVCR2", "TIGIT", "FOXP3", "IL2RA", "TCF7", "SELL",
              "CCR7", "IL7R", "KLRB1", "SLC4A10", "NKG7", "GNLY"]:
        if g in agg.index:
            r = agg.loc[g]
            print(f"    {g:9s} {'HVG' if r['is_hvg'] else '   '}  "
                  f"detected in {r['pct_cells_detected']:5.1f}% of cells")
        else:
            print(f"    {g:9s} NOT IN THE COMMON GENE AXIS")
    print(f"\nwrote {HARM}/")
    return meta


@app.local_entrypoint()
def main():
    n = stats_all.remote()
    print(f"\nper-library stats done for {n} libraries; harmonising\n")
    meta = harmonise.remote()
    print(json.dumps(meta, indent=1))
