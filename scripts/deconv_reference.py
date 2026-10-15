"""Build a cell-type signature matrix from the annotated GSE277117 object.

Project 2a, step 1. The single-cell object is 314k cells and needs far more RAM
than the laptop has, so the heavy step runs on Modal and exports two small
tables that the local deconvolution reads:

  signature_matrix.tsv   marker genes x 23 cell types, mean CPM profile
  reference_props.tsv    true cell-type proportions per patient (for validation)

It also exports held-out pseudobulk: patients are split in half, the signature
is built from one half and pseudobulk from the other, so the local step can be
validated against proportions that are actually known before it is trusted on
bulk cohorts where they are not.

    modal run scripts/deconv_reference.py::peek
    modal run --detach scripts/deconv_reference.py::build
"""
import modal

H5AD = "/data/annotated_clean/as_annotated.h5ad"
OUT = "/data/deconv"

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("scanpy==1.12.4", "anndata", "h5py", "pandas", "scikit-learn",
                 "pyarrow", "numpy", "scipy")
)

# Candidate names for the two obs columns we need, most specific first.
CELLTYPE_KEYS = [
    "ct_fine", "majority_voting", "cell_type", "celltype", "predicted_labels",
]
PATIENT_KEYS = ["Patient", "patient", "Sample_num", "sample", "donor"]


def _pick(obs, candidates, what):
    for k in candidates:
        if k in obs.columns:
            return k
    raise KeyError(f"no {what} column found; have {list(obs.columns)}")


@app.function(image=image, volumes={"/data": vol}, memory=8192, timeout=1800)
def peek():
    """What is actually in the annotated object."""
    import anndata as ad
    import pandas as pd

    pd.set_option("display.width", 220)
    a = ad.read_h5ad(H5AD, backed="r")
    print(f"shape: {a.n_obs:,} cells x {a.n_vars:,} genes")
    print(f"layers: {list(a.layers.keys())}")
    print(f"\nobs columns:")
    for c in a.obs.columns:
        s = a.obs[c]
        print(f"  {c:<30} {str(s.dtype):<10} {s.nunique():>5} unique  "
              f"e.g. {list(pd.unique(s.dropna()))[:4]}")
    ct = _pick(a.obs, CELLTYPE_KEYS, "cell type")
    pt = _pick(a.obs, PATIENT_KEYS, "patient")
    print(f"\nusing cell type = {ct!r}, patient = {pt!r}")
    print(f"\ncells per cell type:\n{a.obs[ct].value_counts().to_string()}")
    print(f"\npatients: {a.obs[pt].nunique()}")


@app.function(image=image, volumes={"/data": vol}, memory=65536, timeout=7200)
def build(n_markers_per_type: int = 120, seed: int = 0):
    """Signature matrix + held-out pseudobulk, written to the Volume."""
    import os

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    os.makedirs(OUT, exist_ok=True)

    a = ad.read_h5ad(H5AD)
    ct_key = _pick(a.obs, CELLTYPE_KEYS, "cell type")
    pt_key = _pick(a.obs, PATIENT_KEYS, "patient")
    print(f"{a.n_obs:,} cells x {a.n_vars:,} genes; "
          f"cell type = {ct_key!r}, patient = {pt_key!r}")

    # Raw counts, not the normalised matrix -- a signature matrix must be built
    # on the same scale the mixture is measured on.
    X = a.layers["counts"] if "counts" in a.layers else a.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr()
    print(f"counts matrix: {X.shape}, nnz={X.nnz:,}")

    ct = a.obs[ct_key].astype(str).to_numpy()
    pt = a.obs[pt_key].astype(str).to_numpy()
    genes = np.asarray(a.var_names)

    # Drop any residual non-gene features by index, never by name -- several
    # antibody-capture tags share a symbol with a real gene.
    if "feature_types" in a.var.columns:
        keep = (a.var["feature_types"].to_numpy() == "Gene Expression")
        print(f"feature_types filter: {keep.sum():,}/{len(keep):,} kept")
        X, genes = X[:, keep], genes[keep]

    patients = np.array(sorted(set(pt)))
    patients = patients[patients != "Multiplet"]
    rng.shuffle(patients)
    half = len(patients) // 2
    ref_pats, test_pats = set(patients[:half]), set(patients[half:])
    print(f"{len(patients)} patients -> signature from {len(ref_pats)}, "
          f"held-out pseudobulk from {len(test_pats)}")

    types = np.array(sorted(set(ct)))
    print(f"{len(types)} cell types")

    def profile(mask_cells):
        """Mean CPM profile over the selected cells."""
        sub = X[mask_cells]
        tot = np.asarray(sub.sum(axis=1)).ravel()
        tot[tot == 0] = 1
        # normalise each cell to CPM first so deep cells do not dominate
        inv = sp.diags(1e6 / tot)
        return np.asarray((inv @ sub).mean(axis=0)).ravel()

    in_ref = np.isin(pt, list(ref_pats))

    sig = pd.DataFrame(index=genes, columns=types, dtype=float)
    n_cells = {}
    for t in types:
        m = in_ref & (ct == t)
        n_cells[t] = int(m.sum())
        sig[t] = profile(m) if m.sum() else 0.0
        print(f"  {t:<34} {n_cells[t]:>7,} cells")

    # Marker selection: per cell type, the genes with the highest ratio of this
    # type's mean to the next-highest type's mean. A signature matrix wants
    # genes that separate types, not genes that are merely abundant.
    logsig = np.log2(sig + 1.0)
    chosen: set[str] = set()
    for t in types:
        others = logsig.drop(columns=t).max(axis=1)
        score = logsig[t] - others
        # require the gene to be meaningfully expressed in its own type
        ok = sig[t] > 1.0
        chosen |= set(score[ok].sort_values(ascending=False).head(n_markers_per_type).index)
    markers = sorted(chosen)
    print(f"\nmarker genes: {len(markers)}")

    sig_m = sig.loc[markers]
    sig_m.to_csv(f"{OUT}/signature_matrix.tsv", sep="\t")

    pd.Series(n_cells).to_csv(f"{OUT}/signature_ncells.tsv", sep="\t", header=["n_cells"])

    # Held-out pseudobulk: sum raw counts per patient, plus that patient's true
    # proportions. The whole transcriptome is kept so the local step can subset.
    rows, props = {}, {}
    for p in sorted(test_pats):
        m = pt == p
        if m.sum() < 500:
            continue
        rows[p] = np.asarray(X[m].sum(axis=0)).ravel()
        vc = pd.Series(ct[m]).value_counts(normalize=True)
        props[p] = vc.reindex(types).fillna(0.0)
    pb = pd.DataFrame(rows, index=genes)
    pb.to_csv(f"{OUT}/heldout_pseudobulk.tsv.gz", sep="\t", compression="gzip")
    pd.DataFrame(props).to_csv(f"{OUT}/heldout_true_props.tsv", sep="\t")
    print(f"held-out pseudobulk: {pb.shape[1]} patients x {pb.shape[0]:,} genes")

    # True proportions for every patient, for reference
    full = (pd.crosstab(pd.Series(pt, name="patient"), pd.Series(ct, name="type"),
                        normalize="index"))
    full.to_csv(f"{OUT}/reference_props.tsv", sep="\t")

    vol.commit()
    print(f"\nwrote to {OUT}: signature_matrix.tsv ({sig_m.shape}), "
          f"heldout_pseudobulk.tsv.gz, heldout_true_props.tsv, reference_props.tsv")


@app.function(image=image, volumes={"/data": vol}, memory=4096, timeout=900)
def fetch(name: str):
    """Return one exported file's bytes so it can be saved locally."""
    with open(f"{OUT}/{name}", "rb") as fh:
        return fh.read()


@app.local_entrypoint()
def main():
    """Build on Modal, then pull the small tables down to data/deconv/."""
    import pathlib

    build.remote()
    out = pathlib.Path(__file__).resolve().parents[1] / "data" / "deconv"
    out.mkdir(parents=True, exist_ok=True)
    for name in ("signature_matrix.tsv", "signature_ncells.tsv",
                 "heldout_true_props.tsv", "reference_props.tsv",
                 "heldout_pseudobulk.tsv.gz"):
        (out / name).write_bytes(fetch.remote(name))
        print(f"saved {out / name}")
