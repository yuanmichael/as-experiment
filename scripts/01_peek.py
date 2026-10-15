"""Open one 10x .h5 count matrix and describe what is actually inside it.

Usage: python -I scripts/01_peek.py <path-to-filtered_feature_bc_matrix.h5>
"""
import sys

import numpy as np
import scanpy as sc

path = sys.argv[1]
adata = sc.read_10x_h5(path)
adata.var_names_make_unique()

print(f"file: {path}")
print(f"shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"\nfirst 3 cell barcodes: {list(adata.obs_names[:3])}")
print(f"first 5 gene names:   {list(adata.var_names[:5])}")
print(f"\nmatrix type: {type(adata.X)}  dtype: {adata.X.dtype}")

# How empty is it?
nnz = adata.X.nnz
total = adata.n_obs * adata.n_vars
print(f"non-zero entries: {nnz:,} of {total:,} ({100 * nnz / total:.1f}% filled)")

# Per-cell sequencing depth: the reason raw counts are not comparable across cells
depth = np.asarray(adata.X.sum(axis=1)).ravel()
genes_per_cell = np.asarray((adata.X > 0).sum(axis=1)).ravel()
print("\ntotal UMI counts per cell:")
for q in (0, 10, 50, 90, 100):
    print(f"  {q:3d}th percentile: {np.percentile(depth, q):>10,.0f}")
print(f"  ratio deepest/shallowest: {depth.max() / depth.min():.0f}x")
print(f"\ngenes detected per cell: median {np.median(genes_per_cell):,.0f} "
      f"(of {adata.n_vars:,} measured)")

# The gene the paper cares about
for gene in ("AIF1", "IFI6", "ISG15", "MX1"):
    if gene in adata.var_names:
        col = adata[:, gene].X
        vals = np.asarray(col.todense()).ravel()
        print(f"\n{gene}: detected in {(vals > 0).sum():,}/{adata.n_obs:,} cells "
              f"({100 * (vals > 0).mean():.1f}%), max count in one cell = {vals.max():.0f}")
