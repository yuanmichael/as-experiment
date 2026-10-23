"""Did the compression lose the exhaustion programme, or can my probe just not find it?

A linear ridge probe on the stage-1 latent recovered the exhaustion score at
r = 0.50, against 0.88 for naive/memory. I read that as the 40x bottleneck
discarding a sparse programme. That conclusion was not established, because two
other explanations fit the same number:

  * the probe is LINEAR. If exhaustion is encoded as a curved combination of
    the latent dimensions, ridge misses it and reports a low r whatever is
    actually in there.
  * the TARGET is crude. The score is the mean of 9 log-normalised genes, and
    most cells have all 9 at zero, so it is a spike at 0 with a thin tail.
    Correlating against a near-constant target has a low ceiling no matter what
    you predict from, and that ceiling would be mine rather than the model's.

Two tests separate them, neither requiring retraining:

  1. ROUND TRIP, no probe at all. Push cells through the frozen encoder AND
     decoder, compute the programme score from the RECONSTRUCTED genes, and
     correlate with the score from the real genes. This asks what the latent
     preserved using the model's own decoder, with no linearity assumption and
     no probe of mine in the path. It is the better test.

  2. NONLINEAR PROBE. Same inputs as the ridge, a small MLP instead. If the MLP
     is far ahead, the information was there and linearity was the limit.

A third reference is included for context: a probe restricted to cells where
the programme is actually expressed. If r rises sharply there, the low overall
number was being driven by the mass of all-zero cells.

    modal run scripts/vdj_latent_probe.py::diagnose --z-dim 100
"""

import json

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
DS = "/data/vdj/_dataset"
LAT = "/data/vdj/_expr_latent"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("torch==2.5.1", "pandas==2.2.3", "numpy==2.1.3",
                 "scipy==1.14.1", "pyarrow==18.1.0", "scikit-learn==1.5.2")
)

TEST_STUDY = "GSE288581"
STUDY_HOLDOUT = "GSE329527"
VAL_SUBJECT_FRAC = 0.15

PROGRAMMES = {
    "exhaustion": ["PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
                   "ENTPD1", "CD160", "BTLA"],
    "cytotoxicity": ["GZMA", "GZMB", "GZMH", "PRF1", "GNLY", "NKG7", "IFNG",
                     "KLRG1", "FGFBP2"],
    "naive_memory": ["TCF7", "SELL", "CCR7", "LEF1", "IL7R", "CD27", "CD28"],
    "regulatory": ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "CTLA4"],
}


@app.function(image=image, volumes={"/data": vol}, gpu="L4", cpu=8.0,
              memory=65536, timeout=7200)
def diagnose(z_dim: int = 100, lib_emb: int = 16, n_cells: int = 60000):
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import torch
    import torch.nn as nn
    from sklearn.linear_model import Ridge
    from sklearn.neural_network import MLPRegressor

    vol.reload()
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    X = sp.load_npz(f"{DS}/X.npz").tocsr()
    cells = pd.read_parquet(f"{DS}/cells.parquet")
    genes = open(f"{DS}/genes.txt").read().split()
    gidx = {g: i for i, g in enumerate(genes)}
    n_genes = X.shape[1]

    lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
    lib[lib == 0] = 1.0
    lib_codes, lib_uniq = pd.factorize(cells["library"])
    n_libs = len(lib_uniq)

    # Same split as training, so this is measured on cells the VAE never saw.
    gse, subj = cells["gse"].values, cells["subject"].values
    pool = ~((gse == TEST_STUDY) | (gse == STUDY_HOLDOUT))
    tsub = pd.unique(subj[pool])
    rng = np.random.default_rng(1234)
    val_s = set(rng.choice(tsub, max(1, int(VAL_SUBJECT_FRAC * len(tsub))),
                           replace=False))
    is_val = pool & np.isin(subj, list(val_s))
    idx = np.where(is_val)[0]
    if len(idx) > n_cells:
        idx = np.random.default_rng(0).choice(idx, n_cells, replace=False)
        idx.sort()
    print(f"z_dim={z_dim}  {len(idx):,} held-out cells")

    # Architecture must match vdj_expr_vae.py exactly or the weights will not load.
    class ExprVAE(nn.Module):
        def __init__(self):
            super().__init__()
            self.lib = nn.Embedding(n_libs, lib_emb)
            self.enc = nn.Sequential(
                nn.Linear(n_genes + lib_emb, 512), nn.LayerNorm(512), nn.GELU(),
                nn.Linear(512, 256), nn.LayerNorm(256), nn.GELU())
            self.mu = nn.Linear(256, z_dim)
            self.logvar = nn.Linear(256, z_dim)
            self.dec = nn.Sequential(
                nn.Linear(z_dim + lib_emb, 256), nn.LayerNorm(256), nn.GELU(),
                nn.Linear(256, 512), nn.LayerNorm(512), nn.GELU(),
                nn.Linear(512, n_genes))
            self.log_theta = nn.Parameter(torch.zeros(n_genes))

    model = ExprVAE().to(dev)
    sd = torch.load(f"{LAT}/expr_vae_{z_dim}.pt", map_location=dev,
                    weights_only=True)
    model.load_state_dict(sd)
    model.eval()

    # ---- one pass: latent, and the decoder's reconstructed expression ----
    Zs, Rs = [], []
    with torch.no_grad():
        for i in range(0, len(idx), 4096):
            b = idx[i:i + 4096]
            raw = torch.tensor(X[b].toarray(), dtype=torch.float32, device=dev)
            l = torch.tensor(lib[b], device=dev)
            bc = torch.tensor(lib_codes[b], device=dev)
            e = model.lib(bc)
            h = model.enc(torch.cat([torch.log1p(1e4 * raw / l.unsqueeze(-1)),
                                     e], -1))
            mu = model.mu(h)
            rate = torch.softmax(model.dec(torch.cat([mu, e], -1)), -1)
            Zs.append(mu.cpu().numpy())
            # rate is the predicted fraction of transcripts per gene, so
            # 1e4 * rate is already CPM-to-10k -- directly comparable to the
            # normalisation applied to the observed counts below.
            Rs.append(torch.log1p(1e4 * rate).cpu().numpy())
    Z = np.concatenate(Zs)
    Rec = np.concatenate(Rs)

    Y = X[idx].toarray()
    d = Y.sum(1, keepdims=True)
    d[d == 0] = 1
    Obs = np.log1p(1e4 * Y / d)

    def corr(a, b):
        if a.std() < 1e-12 or b.std() < 1e-12:
            return None
        return round(float(np.corrcoef(a, b)[0, 1]), 4)

    out = {"z_dim": z_dim, "n_cells": len(idx)}
    print(f"\n{'programme':16s} {'ridge':>8s} {'MLP':>8s} {'ROUND TRIP':>11s} "
          f"{'expressed':>10s} {'pct>0':>7s}")
    print("-" * 68)

    for name, gs in PROGRAMMES.items():
        cols = [gidx[g] for g in gs if g in gidx]
        obs_score = Obs[:, cols].mean(1)
        rec_score = Rec[:, cols].mean(1)
        half = len(Z) // 2

        ridge = corr(
            Ridge(alpha=1.0).fit(Z[:half], obs_score[:half]).predict(Z[half:]),
            obs_score[half:])

        mlp = MLPRegressor(hidden_layer_sizes=(128, 64), max_iter=60,
                           early_stopping=True, random_state=0)
        mlp.fit(Z[:half], obs_score[:half])
        mlp_r = corr(mlp.predict(Z[half:]), obs_score[half:])

        # No probe at all: the model's own encode-decode round trip.
        rt = corr(rec_score, obs_score)

        # Restricted to cells that actually express the programme, to see how
        # much of the low number is the mass of all-zero cells.
        nz = obs_score > 0
        pct_nz = 100 * float(nz.mean())
        expressed = None
        if nz.sum() > 2000:
            Zn, sn = Z[nz], obs_score[nz]
            h = len(Zn) // 2
            expressed = corr(
                Ridge(alpha=1.0).fit(Zn[:h], sn[:h]).predict(Zn[h:]), sn[h:])

        out[name] = {"ridge": ridge, "mlp": mlp_r, "round_trip": rt,
                     "ridge_expressed_only": expressed,
                     "pct_cells_expressing": round(pct_nz, 1)}
        print(f"{name:16s} {ridge if ridge is not None else 0:8.3f} "
              f"{mlp_r if mlp_r is not None else 0:8.3f} "
              f"{rt if rt is not None else 0:11.3f} "
              f"{expressed if expressed is not None else 0:10.3f} "
              f"{pct_nz:6.1f}%")

    print("\nROUND TRIP is the test that matters: no probe of mine in the path,")
    print("no linearity assumption. High round trip + low ridge means the latent")
    print("kept the programme and my probe could not read it. Both low means the")
    print("compression really did discard it.")

    with open(f"{LAT}/probe_diagnostics_{z_dim}.json", "w") as fh:
        json.dump(out, fh, indent=1)
    vol.commit()
    return out


@app.local_entrypoint()
def main(z_dim: int = 100):
    print(json.dumps(diagnose.remote(z_dim=z_dim), indent=1))
