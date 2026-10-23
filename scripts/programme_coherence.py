"""Are my marker "programmes" actually single programmes?

Each programme score is the mean of N normalised gene values. That only measures
one thing if those genes move together. If they do not, the average blends
unrelated signals, and no compact representation can predict it -- not because
the representation is poor but because the target is not a quantity.

This is the suspicion raised by the round-trip diagnostic: exhaustion is
expressed in 41% of cells and recovers at r = 0.60, while regulatory is
expressed in 34% and recovers at 0.70. More widely expressed and harder to
recover, so sparsity does not explain it. A plausible reason is that my
exhaustion set mixes contexts -- CTLA4 is a core regulatory-T-cell gene, TIGIT
sits on NK cells, BTLA on B cells.

Reports per programme:
  mean pairwise correlation between its genes  (does it cohere?)
  first-eigenvalue share                       (is there one dominant axis?)
  each gene's loading on that axis             (which genes belong, which don't)
  cross-programme score correlation            (are the programmes distinct?)

Local, CPU, no cost -- reads the dataset already on the volume via a small
Modal function because the matrix lives there.

    modal run scripts/programme_coherence.py::check
"""

import json

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
DS = "/data/vdj/_dataset"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("pandas==2.2.3", "numpy==2.1.3", "scipy==1.14.1",
                 "pyarrow==18.1.0")
)

PROGRAMMES = {
    "exhaustion": ["PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
                   "ENTPD1", "CD160", "BTLA"],
    "cytotoxicity": ["GZMA", "GZMB", "GZMH", "PRF1", "GNLY", "NKG7", "IFNG",
                     "KLRG1", "FGFBP2"],
    "naive_memory": ["TCF7", "SELL", "CCR7", "LEF1", "IL7R", "CD27", "CD28"],
    "regulatory": ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "CTLA4"],
}


@app.function(image=image, volumes={"/data": vol}, cpu=8.0, memory=65536,
              timeout=3600)
def check(n_cells: int = 150000):
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    vol.reload()
    X = sp.load_npz(f"{DS}/X.npz").tocsr()
    genes = open(f"{DS}/genes.txt").read().split()
    gidx = {g: i for i, g in enumerate(genes)}

    rng = np.random.default_rng(0)
    idx = np.sort(rng.choice(X.shape[0], min(n_cells, X.shape[0]),
                             replace=False))
    Y = X[idx].toarray()
    d = Y.sum(1, keepdims=True)
    d[d == 0] = 1
    Yn = np.log1p(1e4 * Y / d)
    print(f"{len(idx):,} cells\n")

    out, score_cols = {}, {}
    for name, gs in PROGRAMMES.items():
        present = [g for g in gs if g in gidx]
        M = Yn[:, [gidx[g] for g in present]]
        C = np.corrcoef(M, rowvar=False)
        off = C[~np.eye(len(present), dtype=bool)]

        # One dominant axis? Eigen-decompose the correlation matrix.
        w, v = np.linalg.eigh(C)
        order = np.argsort(w)[::-1]
        w, v = w[order], v[:, order]
        pc1_share = float(w[0] / w.sum())
        load = v[:, 0]
        if load.sum() < 0:
            load = -load

        score_cols[name] = M.mean(1)
        out[name] = {
            "n_genes": len(present),
            "mean_pairwise_r": round(float(off.mean()), 4),
            "min_pairwise_r": round(float(off.min()), 4),
            "max_pairwise_r": round(float(off.max()), 4),
            "pc1_share": round(pc1_share, 4),
            "loadings": {g: round(float(l), 3)
                         for g, l in zip(present, load)},
        }

        print(f"{name.upper():16s} {len(present)} genes   "
              f"mean pairwise r = {off.mean():.3f}   "
              f"PC1 explains {100*pc1_share:.1f}%")
        for g, l in sorted(zip(present, load), key=lambda t: -t[1]):
            flag = "  <-- barely belongs" if l < 0.20 else ""
            print(f"    {g:9s} loading {l:+.3f}{flag}")
        print(f"    weakest gene pair r = {off.min():+.3f}\n")

    print("=" * 72)
    print("CROSS-PROGRAMME SCORE CORRELATION")
    print("(high values mean the scores are not measuring distinct things)")
    print("=" * 72)
    names = list(score_cols)
    print("            " + "".join(f"{n[:9]:>11s}" for n in names))
    cross = {}
    for a in names:
        row = []
        for b in names:
            r = float(np.corrcoef(score_cols[a], score_cols[b])[0, 1])
            row.append(r)
            cross[f"{a}|{b}"] = round(r, 4)
        print(f"  {a[:10]:10s}" + "".join(f"{x:>11.3f}" for x in row))
    out["cross_programme"] = cross

    print("\nmean pairwise r near 0 means the genes are not co-regulated, so")
    print("their average is a blend of unrelated signals rather than one axis.")
    print("A low PC1 share says the same thing: no single dominant direction.")

    with open(f"{DS}/programme_coherence.json", "w") as fh:
        json.dump(out, fh, indent=1)
    vol.commit()
    return out


@app.local_entrypoint()
def main():
    check.remote()
