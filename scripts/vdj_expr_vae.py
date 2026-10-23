"""Stage 1: compress 2,009 genes into a learned latent, using expression ALONE.

This model never sees a TCR, a gene list, a cell type label or a diagnosis. Its
only objective is to rebuild each cell's counts from a narrow bottleneck. The
50 dimensions are not defined by anyone; they are whatever the compression
needs them to be, which is coordinated gene programmes, because genes that move
together are cheaper to store together.

Why a separate stage at all. The stage-2 experiment measures how much of a cell
state the receptor explains, and that is only a clean measurement if the target
space is fixed BEFORE the TCR is involved. Trained jointly, the shared latent
would be dominated by the expression branch -- expression predicts expression
far better than a receptor does -- and the TCR branch would wither into an
input the model ignores while still appearing to be half the architecture.
That is PLAPT's prot_zero pathology: a branch contributing +0.005 that looked
load-bearing until it was deleted.

Two things are measured rather than assumed:

  * EFFECTIVE RANK of the latent. Asking for 50 dimensions does not mean
    getting 50. plapt-ablation found ProtBERT's pooler_output had nominal
    dimension 1024 and effective rank 1.00, with every pair of unrelated
    proteins at cosine similarity 1.0000. If this latent asks for 50 and
    delivers 8, then 8 is the real answer.

  * Whether the latent is just encoding BATCH. A library embedding is given to
    both encoder and decoder so the model can account for batch without
    spending latent capacity on it, and the latent is then probed for library
    identity directly. A latent that predicts which library a cell came from
    is a batch detector, and stage 2 would score well by predicting batch.

    modal run scripts/vdj_expr_vae.py::train_expr --z-dim 50
    modal run scripts/vdj_expr_vae.py::sweep
"""

import json

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
DS = "/data/vdj/_dataset"
OUT = "/data/vdj/_expr_latent"

# Versions pinned. The first attempt at this run died on
# LogisticRegression(multi_class=...), an argument that existed when the code
# was written and had been removed by the time it ran.
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


def nb_nll(x, mu, theta, eps=1e-8):
    import torch
    mu = mu.clamp(min=eps)
    theta = theta.clamp(min=eps).unsqueeze(0)
    log_tm = (theta + mu).log()
    ll = (torch.lgamma(x + theta) - torch.lgamma(theta) - torch.lgamma(x + 1)
          + theta * (theta.log() - log_tm) + x * (mu.log() - log_tm))
    return -ll.sum(-1)


def effective_rank(Z):
    """exp of the entropy of the normalised singular value spectrum.

    Reports how many dimensions the representation ACTUALLY uses, as opposed to
    how many it nominally has. A latent spread evenly over 50 directions scores
    ~50; one that has collapsed onto a single direction scores ~1.
    """
    import numpy as np
    Zc = Z - Z.mean(0, keepdims=True)
    s = np.linalg.svd(Zc, compute_uv=False)
    p = s / max(s.sum(), 1e-12)
    p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


@app.function(image=image, volumes={"/data": vol}, gpu="L4", cpu=8.0,
              memory=65536, timeout=14400)
def train_expr(z_dim: int = 50, epochs: int = 12, beta: float = 1.0,
               batch: int = 512, lr: float = 1e-3, lib_emb: int = 16,
               save_latent: bool = True, tag: str = ""):
    import os
    import time

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import torch
    import torch.nn as nn
    from sklearn.linear_model import LogisticRegression, Ridge

    vol.reload()
    torch.manual_seed(0)
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    X = sp.load_npz(f"{DS}/X.npz").tocsr()
    cells = pd.read_parquet(f"{DS}/cells.parquet")
    genes = open(f"{DS}/genes.txt").read().split()
    gidx = {g: i for i, g in enumerate(genes)}
    n_genes = X.shape[1]
    print(f"z_dim={z_dim} beta={beta} | {X.shape[0]:,} cells x {n_genes:,} genes")

    lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
    lib[lib == 0] = 1.0
    lib_codes, lib_uniq = pd.factorize(cells["library"])
    n_libs = len(lib_uniq)

    gse, subj = cells["gse"].values, cells["subject"].values
    is_test, is_study = gse == TEST_STUDY, gse == STUDY_HOLDOUT
    pool = ~(is_test | is_study)
    tsub = pd.unique(subj[pool])
    rng = np.random.default_rng(1234)
    val_s = set(rng.choice(tsub, max(1, int(VAL_SUBJECT_FRAC * len(tsub))),
                           replace=False))
    is_val = pool & np.isin(subj, list(val_s))
    is_train = pool & ~is_val
    print(f"  train {int(is_train.sum()):,}  val {int(is_val.sum()):,}  "
          f"study {int(is_study.sum()):,}  test {int(is_test.sum()):,}")

    class ExprVAE(nn.Module):
        def __init__(self):
            super().__init__()
            # The library embedding lets the model account for batch WITHOUT
            # spending latent capacity on it. Without this the first few latent
            # dimensions would just encode which sample a cell came from.
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

        def encode(self, x, b):
            # log1p on CPM-normalised counts for the ENCODER input only; the
            # decoder still models raw counts, so no information is lost from
            # the likelihood.
            h = self.enc(torch.cat([x, self.lib(b)], -1))
            return self.mu(h), self.logvar(h).clamp(-8, 8)

        def forward(self, x_norm, x_raw, b, l):
            mu, logvar = self.encode(x_norm, b)
            z = mu + torch.randn_like(mu) * (0.5 * logvar).exp() \
                if self.training else mu
            rate = torch.softmax(self.dec(torch.cat([z, self.lib(b)], -1)), -1)
            return rate * l.unsqueeze(-1), self.log_theta.exp(), mu, logvar

    model = ExprVAE().to(dev)
    print(f"  {sum(p.numel() for p in model.parameters()):,} parameters")
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    tr_idx = np.where(is_train)[0]

    def batches(idx, size):
        for i in range(0, len(idx), size):
            b = idx[i:i + size]
            raw = torch.tensor(X[b].toarray(), dtype=torch.float32, device=dev)
            l = torch.tensor(lib[b], device=dev)
            norm = torch.log1p(1e4 * raw / l.unsqueeze(-1))
            yield b, raw, norm, l, torch.tensor(lib_codes[b], device=dev)

    def encode_all(mask):
        model.eval()
        idx = np.where(mask)[0]
        zs = []
        with torch.no_grad():
            for _, raw, norm, l, bc in batches(idx, 4096):
                mu, _ = model.encode(norm, bc)
                zs.append(mu.cpu().numpy())
        return idx, (np.concatenate(zs) if zs else np.zeros((0, z_dim)))

    def eval_nll(mask):
        model.eval()
        tot, n = 0.0, 0
        with torch.no_grad():
            for _, raw, norm, l, bc in batches(np.where(mask)[0], 4096):
                m, th, _, _ = model(norm, raw, bc, l)
                tot += nb_nll(raw, m, th).sum().item()
                n += len(raw)
        return tot / max(n, 1)

    hist, best, t0 = [], None, time.time()
    for ep in range(epochs):
        model.train()
        rng.shuffle(tr_idx)
        run, nb = 0.0, 0
        for _, raw, norm, l, bc in batches(tr_idx, batch):
            m, th, mu, logvar = model(norm, raw, bc, l)
            rec = nb_nll(raw, m, th).mean()
            kl = (-0.5 * (1 + logvar - mu.pow(2) - logvar.exp()).sum(-1)).mean()
            loss = rec + beta * kl
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            run += loss.item()
            nb += 1
        v = eval_nll(is_val)
        hist.append({"epoch": ep, "train": run / max(nb, 1), "val_nll": v,
                     "kl": float(kl.item()),
                     "mins": round((time.time() - t0) / 60, 1)})
        print(f"  ep{ep:02d} train {run/max(nb,1):9.1f}  val_nll {v:9.2f}  "
              f"kl {kl.item():7.2f}  {hist[-1]['mins']:5.1f}m")
        if best is None or v < best:
            best = v
            # To the VOLUME, not /tmp: a container that dies or is preempted
            # takes /tmp with it, and the point of checkpointing is to survive
            # exactly that.
            os.makedirs(OUT, exist_ok=True)
            torch.save(model.state_dict(), f"{OUT}/ckpt_{z_dim}.pt")
            vol.commit()

    model.load_state_dict(torch.load(f"{OUT}/ckpt_{z_dim}.pt",
                                     weights_only=True))

    res = {"z_dim": z_dim, "beta": beta, "epochs": epochs, "history": hist,
           "val_nll": eval_nll(is_val), "study_nll": eval_nll(is_study),
           "test_nll": eval_nll(is_test)}

    # ------------------------------------------------------------ SAVE FIRST
    # The first attempt at this ran all 12 epochs and then threw in the
    # diagnostics below, discarding the entire run. Weights and latents are
    # the expensive artefact; persist them before computing anything optional.
    if save_latent:
        os.makedirs(OUT, exist_ok=True)
        allz = np.zeros((X.shape[0], z_dim), dtype=np.float32)
        for mask in (is_train, is_val, is_study, is_test):
            idx, Z = encode_all(mask)
            allz[idx] = Z
        np.save(f"{OUT}/z_expr_{z_dim}.npy", allz)
        torch.save(model.state_dict(), f"{OUT}/expr_vae_{z_dim}.pt")
        vol.commit()
        print(f"  saved latent for all {X.shape[0]:,} cells")

    # ---------------------------------------------------- is the latent real?
    vi, Zv = encode_all(is_val)
    res["effective_rank"] = round(effective_rank(Zv), 2)
    res["nominal_dim"] = z_dim
    Zn = Zv / np.maximum(np.linalg.norm(Zv, axis=1, keepdims=True), 1e-9)
    samp = Zn[np.random.default_rng(0).choice(len(Zn), min(3000, len(Zn)),
                                              replace=False)]
    cos = samp @ samp.T
    res["mean_pairwise_cosine"] = round(
        float((cos.sum() - np.trace(cos)) / (len(cos) ** 2 - len(cos))), 4)

    # Does the latent carry real biology, and does it carry batch?
    Yv = X[vi].toarray()
    depth = Yv.sum(1, keepdims=True)
    depth[depth == 0] = 1
    Yn = np.log1p(1e4 * Yv / depth)
    half = len(Zv) // 2
    for name, gs in PROGRAMMES.items():
        cols = [gidx[g] for g in gs if g in gidx]
        score = Yn[:, cols].mean(1)
        pred = Ridge(alpha=1.0).fit(Zv[:half], score[:half]).predict(Zv[half:])
        res[f"latent_{name}_r"] = round(
            float(np.corrcoef(pred, score[half:])[0, 1]), 4)

    # A latent that identifies the library is a batch detector. Chance is
    # 1/n_libraries, so compare against that, not against zero.
    # Wrapped: a diagnostic that fails must not take the run down with it.
    try:
        labs = lib_codes[vi]
        keep = pd.Series(labs).map(pd.Series(labs).value_counts()) >= 50
        if keep.sum() > 500:
            Zk, lk = Zv[keep.values], labs[keep.values]
            h = len(Zk) // 2
            clf = LogisticRegression(max_iter=300).fit(Zk[:h], lk[:h])
            acc = float((clf.predict(Zk[h:]) == lk[h:]).mean())
            res["library_probe_acc"] = round(acc, 4)
            res["library_probe_chance"] = round(1.0 / len(np.unique(lk)), 4)
    except Exception as exc:
        res["library_probe_error"] = f"{type(exc).__name__}: {exc}"
        print(f"  library probe failed: {res['library_probe_error']}")

    if save_latent:
        with open(f"{OUT}/meta_{z_dim}.json", "w") as fh:
            json.dump({k: v for k, v in res.items() if k != "history"}, fh,
                      indent=1)
        vol.commit()

    print("\n" + json.dumps({k: v for k, v in res.items()
                             if k != "history"}, indent=1))
    return res


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def report():
    import glob
    vol.reload()
    rows = []
    for p in sorted(glob.glob(f"{OUT}/meta_*.json")):
        with open(p) as fh:
            rows.append(json.load(fh))
    hdr = (f"{'z':>4s} {'eff rank':>9s} {'val nll':>10s} {'cos':>7s} "
           f"{'lib acc':>8s} {'chance':>7s} {'exh':>6s} {'cyt':>6s} {'naive':>6s}")
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda d: d["z_dim"]):
        print(f"{r['z_dim']:4d} {r.get('effective_rank',0):9.2f} "
              f"{r.get('val_nll',0):10.2f} {r.get('mean_pairwise_cosine',0):7.3f} "
              f"{r.get('library_probe_acc',0):8.3f} "
              f"{r.get('library_probe_chance',0):7.4f} "
              f"{r.get('latent_exhaustion_r',0):6.3f} "
              f"{r.get('latent_cytotoxicity_r',0):6.3f} "
              f"{r.get('latent_naive_memory_r',0):6.3f}")
    print("\neff rank << z means the latent is narrower than it looks.")
    print("lib acc >> chance means it is partly a batch detector.")
    return rows
