"""Consolidate 291 libraries into one training array on the harmonised gene axis.

Turns per-library output into the arrays a model can actually be fed:

    X.npz        sparse counts, cells x 2,009 harmonised genes
    cells.parquet  cdr3_aa, v_gene, j_gene, clone_size, gse, library, subject
    vocab.json   amino acid, V gene and J gene vocabularies

Three things done here rather than in the model, because getting them wrong is
invisible once training starts:

1. Counts are kept as COUNTS. The decoder models them with a negative binomial,
   which needs integers and the per-cell library size. Log-normalising here
   would throw away the information the likelihood depends on.

2. Clone size is recomputed GLOBALLY per study, not per library. A clone split
   across two libraries of the same subject would otherwise look like two
   separate singletons, which would corrupt the one evaluation stratum we
   actually trust.

3. A subject id is derived, and it is a HEURISTIC. GEO does not reliably record
   which libraries came from the same person, and a cell-level split would let
   the model recognise a donor and look up the answer. Where the naming makes
   replicates obvious (GSE301133 deposits ...1058a and ...1058b) they are
   merged; otherwise subject falls back to library, which is the conservative
   direction for leakage within a study but cannot catch a donor appearing
   under two unrelated names. For that reason the HEADLINE split is
   leave-one-study-out, which is immune to the heuristic being wrong.

    modal run scripts/vdj_dataset.py::build
"""

import json
import re

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
VDJ = "/data/vdj"
HARM = "/data/vdj/_harmonised"
DS = "/data/vdj/_dataset"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("pandas", "numpy", "scipy", "pyarrow")
)

DROP = {"GSE256141", "GSE292177"}

AA = "ACDEFGHIKLMNPQRSTVWY"
MAX_CDR3 = 30          # covers essentially all TRB CDR3s; longer are truncated

# Trailing tokens that mark a technical replicate or lane of the SAME subject
# rather than a different one. Deliberately narrow -- merging two real subjects
# is worse than failing to merge two replicates, because the former inflates
# the score and the latter only wastes a little data.
REPLICATE_SUFFIX = re.compile(r"(?<=\d)[ab]$|_(?:r|rep|lane|run)?[12]$", re.I)


def subject_of(gse: str, library: str) -> str:
    """Best guess at which person a library came from. See module docstring."""
    s = library
    if gse == "GSE301133":
        # tbt04301058a / tbt04301058b are the same subject
        s = re.sub(r"[ab]$", "", s)
    elif gse == "GSE271413":
        # 1313_5, 1324_tcb, 1490_1 -> leading numeric id is the subject
        m = re.match(r"(\d+)", s)
        if m:
            s = m.group(1)
    elif gse == "GSE272993":
        # 2021_04_29_1_4 -> date_cohort identifies the run; the trailing index
        # is the lane. Keep everything but the final token.
        parts = s.split("_")
        if len(parts) > 2:
            s = "_".join(parts[:-1])
    elif gse == "GSE288581":
        # SF73-4, SF76, HC1788 -> strip a trailing fraction index
        s = re.sub(r"[-_]\d$", "", s)
    return f"{gse}:{s}"


@app.function(image=image, volumes={"/data": vol}, cpu=8.0, memory=65536,
              timeout=10800)
def build():
    import glob
    import os

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    vol.reload()
    hvg_file = sorted(glob.glob(f"{HARM}/hvg_*.txt"))[0]
    with open(hvg_file) as fh:
        hvg = [g.strip() for g in fh.read().splitlines() if g.strip()]
    hvg_idx = {g: i for i, g in enumerate(hvg)}
    print(f"{len(hvg)} target genes from {os.path.basename(hvg_file)}")

    libs = []
    for p in sorted(glob.glob(f"{VDJ}/*/*/cells.parquet")):
        gse, name = p.split("/")[-3], p.split("/")[-2]
        if gse not in DROP:
            libs.append((gse, name, os.path.dirname(p)))
    print(f"{len(libs)} libraries")

    blocks, frames, skipped = [], [], []
    for n, (gse, name, d) in enumerate(libs, 1):
        try:
            cells = pd.read_parquet(f"{d}/cells.parquet")
            with open(f"{d}/genes.txt") as fh:
                genes = [g.strip() for g in fh.read().splitlines()]
            X = sp.load_npz(f"{d}/expr.npz").tocsc()
        except Exception as exc:
            skipped.append((f"{gse}/{name}", f"read: {type(exc).__name__}"))
            continue
        if X.shape[0] != len(cells):
            skipped.append((f"{gse}/{name}",
                            f"rows {X.shape[0]} vs cells {len(cells)}"))
            continue

        # Sum duplicate symbols, exactly as the gene-stats pass did, so the
        # columns mean the same thing here as they did during selection.
        sym = pd.Index(genes)
        if len(sym) != sym.nunique():
            codes, uniq = pd.factorize(sym)
            agg = sp.csr_matrix(
                (np.ones(len(codes)), (np.arange(len(codes)), codes)),
                shape=(len(codes), len(uniq)))
            X = (X @ agg).tocsc()
            genes = list(uniq)

        pos = {g: i for i, g in enumerate(genes)}
        take = [pos[g] for g in hvg if g in pos]
        slot = [hvg_idx[g] for g in hvg if g in pos]
        if len(take) < 0.9 * len(hvg):
            skipped.append((f"{gse}/{name}",
                            f"only {len(take)}/{len(hvg)} target genes present"))
            continue

        sub = X[:, take].tocsc()
        # Place into the canonical column order, leaving absent genes at zero.
        M = sp.lil_matrix((X.shape[0], len(hvg)), dtype=np.float32)
        M[:, slot] = sub
        blocks.append(M.tocsr())

        cells = cells.copy()
        cells["library"] = name
        cells["gse"] = gse
        cells["subject"] = [subject_of(gse, name)] * len(cells)
        # Keep each cell's total counts across ALL genes: the NB decoder needs
        # the true library size, not the total over 2,009 selected genes.
        cells["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
        frames.append(cells)

        if n % 40 == 0:
            print(f"  {n}/{len(libs)} libraries, "
                  f"{sum(b.shape[0] for b in blocks):,} cells so far")

    if not blocks:
        print("nothing assembled")
        return {}

    X = sp.vstack(blocks, format="csr")
    cells = pd.concat(frames, ignore_index=True)
    assert X.shape[0] == len(cells), (X.shape, len(cells))

    # ------------------------------------------------- global clone size
    # Per-library counts miss a clone spread over several libraries of the same
    # study, turning one expanded clone into several fake singletons.
    key = cells["gse"] + "|" + cells["cdr3_aa"] + "|" + cells["v_gene"]
    cells["clone_size_global"] = key.map(key.value_counts())

    # ------------------------------------------------------------ vocab
    aa_vocab = {c: i + 1 for i, c in enumerate(AA)}      # 0 = pad
    v_vocab = {g: i + 1 for i, g in enumerate(
        sorted(cells["v_gene"].astype(str).unique()))}   # 0 = unknown
    j_vocab = {g: i + 1 for i, g in enumerate(
        sorted(cells["j_gene"].astype(str).unique()))}

    os.makedirs(DS, exist_ok=True)
    sp.save_npz(f"{DS}/X.npz", X)
    cells.to_parquet(f"{DS}/cells.parquet", index=False)
    with open(f"{DS}/genes.txt", "w") as fh:
        fh.write("\n".join(hvg))
    with open(f"{DS}/vocab.json", "w") as fh:
        json.dump({"aa": aa_vocab, "v_gene": v_vocab, "j_gene": j_vocab,
                   "max_cdr3": MAX_CDR3}, fh, indent=1)

    meta = dict(
        n_cells=int(X.shape[0]), n_genes=int(X.shape[1]),
        nnz=int(X.nnz), density=round(X.nnz / (X.shape[0] * X.shape[1]), 4),
        n_studies=int(cells["gse"].nunique()),
        n_subjects=int(cells["subject"].nunique()),
        n_libraries=int(cells["library"].nunique()),
        n_clonotypes=int(cells.groupby(["gse", "cdr3_aa", "v_gene"]).ngroups),
        singletons=int((cells["clone_size_global"] == 1).sum()),
        n_v_genes=len(v_vocab), n_j_genes=len(j_vocab),
        median_cdr3_len=int(cells["cdr3_aa"].str.len().median()),
        pct_cdr3_over_max=round(
            100 * float((cells["cdr3_aa"].str.len() > MAX_CDR3).mean()), 3),
        skipped=skipped,
    )
    with open(f"{DS}/meta.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    vol.commit()

    print("\n" + "=" * 88)
    print("DATASET")
    print("=" * 88)
    for k, v in meta.items():
        if k != "skipped":
            print(f"  {k:22s} {v}")
    print(f"\n  cells per study:")
    print(cells["gse"].value_counts().to_string())
    print(f"\n  cells by clone size:")
    cs = cells["clone_size_global"]
    for lo, hi in [(1, 1), (2, 2), (3, 5), (6, 20), (21, 100), (101, 10**9)]:
        n = int(((cs >= lo) & (cs <= hi)).sum())
        label = f"{lo}" if lo == hi else f"{lo}-{hi if hi < 10**9 else '+'}"
        print(f"    clone size {label:>8s}: {n:>10,} cells "
              f"({100*n/len(cs):5.1f}%)")
    if skipped:
        print(f"\n  {len(skipped)} libraries skipped:")
        for lib, why in skipped[:20]:
            print(f"    {lib}: {why}")
    return meta


@app.local_entrypoint()
def main():
    print(json.dumps(build.remote(), indent=1)[:2000])
