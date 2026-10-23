"""Does a T cell's receptor sequence predict its own transcriptome?

A conditional VAE: CDR3 + V + J in, that cell's 2,009-gene count vector out.

    CDR3 (<=30 aa) ---> transformer ---\
    V gene          ---> embedding -----+--> z (50-d) --> decoder --> NB(mu, theta)
    J gene          ---> embedding ----/

Loss = negative-binomial negative log-likelihood + beta * KL. The NB is not
decoration: 81% of entries are zero, variance exceeds the mean (overdispersion),
and sequencing depth varies several-fold per cell. Squared error on
log-normalised counts -- what mvTCR uses -- assumes none of that.

THE POINT OF THIS SCRIPT IS THE SHUFFLED ARM. A reconstruction loss always
falls; predicting the global mean profile for every cell already scores well
because "mostly zero" is an easy guess. So the only question that matters is
whether the TCR input moves the number, and the only way to know is to destroy
the TCR-cell correspondence and retrain. PLAPT looked like it used its protein
branch until the branch was zeroed and the score moved 0.005.

Three splits of increasing difficulty:
    val      15% of subjects, drawn from every training study  (in-distribution)
    study    GSE329527 entirely                                (new cohort)
    test     GSE288581 -- ankylosing spondylitis               (never seen, OOD)

    modal run scripts/vdj_vae.py::train --arm real
    modal run scripts/vdj_vae.py::train --arm shuffled
    modal run scripts/vdj_vae.py::baseline
"""

import json

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
DS = "/data/vdj/_dataset"
RUNS = "/data/vdj/_runs"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("torch==2.5.1", "pandas", "numpy", "scipy", "pyarrow",
                 "scikit-learn")
)

STUDY_HOLDOUT = "GSE329527"     # a whole cohort withheld
TEST_STUDY = "GSE288581"        # AS; out of distribution by construction
VAL_SUBJECT_FRAC = 0.15

# Marker programmes, used ONLY to interpret the latent after the fact via a
# linear probe. They do not enter the loss and do not define the target.
PROGRAMMES = {
    "exhaustion": ["PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
                   "ENTPD1", "CD160", "BTLA"],
    "cytotoxicity": ["GZMA", "GZMB", "GZMH", "PRF1", "GNLY", "NKG7", "IFNG",
                     "KLRG1", "FGFBP2"],
    "naive_memory": ["TCF7", "SELL", "CCR7", "LEF1", "IL7R", "CD27", "CD28"],
    "regulatory": ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "CTLA4"],
}


def build_model(n_v, n_j, n_genes, max_len, d_model=128, n_layers=2,
                n_heads=4, z_dim=50):
    import torch
    import torch.nn as nn

    class Encoder(nn.Module):
        def __init__(self):
            super().__init__()
            self.aa = nn.Embedding(21, d_model, padding_idx=0)
            self.pos = nn.Embedding(max_len, d_model)
            layer = nn.TransformerEncoderLayer(
                d_model, n_heads, d_model * 2, dropout=0.1,
                batch_first=True, norm_first=True)
            self.tr = nn.TransformerEncoder(layer, n_layers)
            self.v = nn.Embedding(n_v + 1, 32, padding_idx=0)
            self.j = nn.Embedding(n_j + 1, 16, padding_idx=0)
            self.head = nn.Sequential(
                nn.Linear(d_model + 48, 256), nn.LayerNorm(256), nn.GELU())
            self.mu = nn.Linear(256, z_dim)
            self.logvar = nn.Linear(256, z_dim)

        def forward(self, seq, v, j):
            pad = seq == 0
            pos = torch.arange(seq.shape[1], device=seq.device)
            h = self.aa(seq) + self.pos(pos)[None]
            h = self.tr(h, src_key_padding_mask=pad)
            # mean over real residues only; padding must not dilute the pool
            keep = (~pad).float().unsqueeze(-1)
            h = (h * keep).sum(1) / keep.sum(1).clamp(min=1)
            h = self.head(torch.cat([h, self.v(v), self.j(j)], dim=-1))
            return self.mu(h), self.logvar(h).clamp(-8, 8)

    class Decoder(nn.Module):
        def __init__(self):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(z_dim, 256), nn.LayerNorm(256), nn.GELU(),
                nn.Linear(256, 512), nn.LayerNorm(512), nn.GELU(),
                nn.Linear(512, n_genes))
            # One dispersion per gene, learned. theta = exp(log_theta).
            self.log_theta = nn.Parameter(torch.zeros(n_genes))

        def forward(self, z, lib):
            # Predict a PROPORTION per gene, then scale by this cell's depth.
            # Predicting counts directly would make the model spend capacity
            # learning sequencing depth, which is not biology.
            rate = torch.softmax(self.net(z), dim=-1)
            return rate * lib.unsqueeze(-1), self.log_theta.exp()

    class VAE(nn.Module):
        def __init__(self):
            super().__init__()
            self.enc, self.dec = Encoder(), Decoder()

        def forward(self, seq, v, j, lib):
            mu, logvar = self.enc(seq, v, j)
            z = mu + torch.randn_like(mu) * (0.5 * logvar).exp() \
                if self.training else mu
            m, th = self.dec(z, lib)
            return m, th, mu, logvar

    return VAE()


def nb_nll(x, mu, theta, eps=1e-8):
    """Negative binomial negative log-likelihood, summed over genes."""
    import torch
    mu = mu.clamp(min=eps)
    theta = theta.clamp(min=eps).unsqueeze(0)
    log_tm = (theta + mu).log()
    ll = (torch.lgamma(x + theta) - torch.lgamma(theta) - torch.lgamma(x + 1)
          + theta * (theta.log() - log_tm) + x * (mu.log() - log_tm))
    return -ll.sum(-1)


@app.function(image=image, volumes={"/data": vol}, gpu="L4", cpu=8.0,
              memory=65536, timeout=14400)
def train(arm: str = "real", epochs: int = 12, beta: float = 1.0,
          z_dim: int = 50, batch: int = 512, lr: float = 1e-3,
          seed: int = 0, tag: str = ""):
    import os
    import time

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import torch
    from sklearn.linear_model import Ridge

    vol.reload()
    torch.manual_seed(seed)
    np.random.seed(seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"arm={arm} beta={beta} z={z_dim} epochs={epochs} device={dev}")

    X = sp.load_npz(f"{DS}/X.npz").tocsr()
    cells = pd.read_parquet(f"{DS}/cells.parquet")
    with open(f"{DS}/vocab.json") as fh:
        vocab = json.load(fh)
    genes = open(f"{DS}/genes.txt").read().split()
    gidx = {g: i for i, g in enumerate(genes)}
    max_len = vocab["max_cdr3"]
    print(f"{X.shape[0]:,} cells x {X.shape[1]:,} genes, {X.nnz:,} nnz")

    # ------------------------------------------------------------- tokenise
    aa = vocab["aa"]
    seqs = np.zeros((len(cells), max_len), dtype=np.int64)
    for i, s in enumerate(cells["cdr3_aa"].astype(str).values):
        for k, c in enumerate(s[:max_len]):
            seqs[i, k] = aa.get(c, 0)
    vg = cells["v_gene"].astype(str).map(vocab["v_gene"]).fillna(0).astype(np.int64).values
    jg = cells["j_gene"].astype(str).map(vocab["j_gene"]).fillna(0).astype(np.int64).values

    # Depth over the MODELLED genes, not all genes -- the decoder's softmax is
    # over 2,009 columns, so scaling by the full library size would bias mu up.
    lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
    lib[lib == 0] = 1.0

    # ---------------------------------------------------------------- splits
    gse = cells["gse"].values
    subj = cells["subject"].values
    is_test = gse == TEST_STUDY
    is_study = gse == STUDY_HOLDOUT
    pool = ~(is_test | is_study)
    train_subjects = pd.unique(subj[pool])
    rng = np.random.default_rng(1234)                 # fixed: split is not a
    val_subjects = set(rng.choice(                    # hyperparameter
        train_subjects, max(1, int(VAL_SUBJECT_FRAC * len(train_subjects))),
        replace=False))
    is_val = pool & np.isin(subj, list(val_subjects))
    is_train = pool & ~is_val
    for name, m in [("train", is_train), ("val", is_val),
                    ("study", is_study), ("test", is_test)]:
        print(f"  {name:6s} {int(m.sum()):>9,} cells  "
              f"{len(pd.unique(subj[m])):>4} subjects  "
              f"{len(pd.unique(gse[m])):>2} studies")

    # ------------------------------------------- the arm that decides anything
    if arm == "shuffled":
        # Permute receptors WITHIN subject, so study and donor structure are
        # preserved and only the TCR-to-cell correspondence is destroyed. A
        # global shuffle would also scramble batch, which would confound the
        # comparison.
        order = np.arange(len(cells))
        for s in pd.unique(subj):
            idx = np.where(subj == s)[0]
            perm = idx.copy()
            rng.shuffle(perm)
            order[idx] = perm
        seqs, vg, jg = seqs[order], vg[order], jg[order]
        print("  TCRs permuted within subject")

    model = build_model(len(vocab["v_gene"]), len(vocab["j_gene"]),
                        X.shape[1], max_len, z_dim=z_dim).to(dev)
    n_par = sum(p.numel() for p in model.parameters())
    print(f"  {n_par:,} parameters")
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)

    tr_idx = np.where(is_train)[0]

    libs_all = cells["library"].values

    def evaluate(mask, probe=False):
        model.eval()
        idx = np.where(mask)[0]
        tot, n = 0.0, 0
        zs, ys, lbs = [], [], []
        with torch.no_grad():
            for i in range(0, len(idx), 2048):
                b = idx[i:i + 2048]
                xb = torch.tensor(X[b].toarray(), dtype=torch.float32, device=dev)
                lb = torch.tensor(lib[b], device=dev)
                m, th, mu, _ = model(
                    torch.tensor(seqs[b], device=dev),
                    torch.tensor(vg[b], device=dev),
                    torch.tensor(jg[b], device=dev), lb)
                tot += nb_nll(xb, m, th).sum().item()
                n += len(b)
                if probe:
                    zs.append(mu.cpu().numpy())
                    ys.append(X[b].toarray())
                    lbs.append(libs_all[b])
        out = {"nll": tot / max(n, 1), "n": n}
        if probe and zs:
            Z = np.concatenate(zs)
            Y = np.concatenate(ys)
            L = np.concatenate(lbs)
            depth = Y.sum(1, keepdims=True)
            depth[depth == 0] = 1
            Yn = np.log1p(1e4 * Y / depth)

            def probe_r(Zs, score):
                half = len(Zs) // 2
                if half < 50:
                    return None
                pred = Ridge(alpha=1.0).fit(Zs[:half], score[:half]).predict(Zs[half:])
                tgt = score[half:]
                if pred.std() == 0 or tgt.std() == 0:
                    return None
                return round(float(np.corrcoef(pred, tgt)[0, 1]), 4)

            # The whole-split probe is confounded wherever the latent can
            # encode WHICH SAMPLE a cell came from: sample identity tracks
            # compartment (synovial fluid vs blood in GSE288581) and
            # compartment tracks exhaustion. Repeating the probe inside single
            # libraries removes that path -- if the correlation only exists
            # across libraries, it was never sequence -> state.
            big = [l for l, c in pd.Series(L).value_counts().items() if c >= 400]
            for name, gs in PROGRAMMES.items():
                cols = [gidx[g] for g in gs if g in gidx]
                if not cols:
                    continue
                score = Yn[:, cols].mean(1)
                r = probe_r(Z, score)
                if r is not None:
                    out[f"probe_{name}_r"] = r
                within = [probe_r(Z[L == l], score[L == l]) for l in big]
                within = [w for w in within if w is not None]
                if within:
                    out[f"probe_{name}_within_lib_r"] = round(
                        float(np.median(within)), 4)
                    out[f"probe_{name}_n_libs"] = len(within)
        return out

    best, hist, t0 = None, [], time.time()
    for ep in range(epochs):
        model.train()
        rng.shuffle(tr_idx)
        run, nb = 0.0, 0
        for i in range(0, len(tr_idx), batch):
            b = tr_idx[i:i + batch]
            xb = torch.tensor(X[b].toarray(), dtype=torch.float32, device=dev)
            lb = torch.tensor(lib[b], device=dev)
            m, th, mu, logvar = model(
                torch.tensor(seqs[b], device=dev),
                torch.tensor(vg[b], device=dev),
                torch.tensor(jg[b], device=dev), lb)
            rec = nb_nll(xb, m, th).mean()
            kl = (-0.5 * (1 + logvar - mu.pow(2) - logvar.exp()).sum(-1)).mean()
            loss = rec + beta * kl
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            run += loss.item()
            nb += 1
        v = evaluate(is_val)
        # Posterior collapse watch: if the latent carries nothing, KL -> 0 and
        # the model is just predicting an average cell.
        hist.append({"epoch": ep, "train_loss": run / max(nb, 1),
                     "val_nll": v["nll"], "kl": float(kl.item()),
                     "mins": round((time.time() - t0) / 60, 1)})
        print(f"  ep{ep:02d} train {run/max(nb,1):9.1f}  val_nll {v['nll']:9.2f}"
              f"  kl {kl.item():7.3f}  {hist[-1]['mins']:5.1f}m")
        if best is None or v["nll"] < best:
            best = v["nll"]
            torch.save(model.state_dict(), f"/tmp/best_{arm}.pt")

    model.load_state_dict(torch.load(f"/tmp/best_{arm}.pt"))
    res = {"arm": arm, "beta": beta, "z_dim": z_dim, "epochs": epochs,
           "seed": seed, "n_params": n_par, "history": hist,
           "val": evaluate(is_val, probe=True),
           "study_holdout": evaluate(is_study, probe=True),
           "test_AS": evaluate(is_test, probe=True)}

    # Stratify the AS test by clone size. Expanded clones are the easy cases;
    # singletons are the honest stratum.
    cs = cells["clone_size_global"].values
    for lo, hi, label in [(1, 1, "singleton"), (2, 5, "2-5"),
                          (6, 100, "6-100"), (101, 10**9, "101+")]:
        m = is_test & (cs >= lo) & (cs <= hi)
        if m.sum() > 200:
            res[f"test_clone_{label}"] = evaluate(m)

    os.makedirs(RUNS, exist_ok=True)
    name = tag or f"{arm}_b{beta}_z{z_dim}_s{seed}"
    with open(f"{RUNS}/{name}.json", "w") as fh:
        json.dump(res, fh, indent=1)
    vol.commit()
    print("\n" + json.dumps({k: v for k, v in res.items()
                             if k != "history"}, indent=1))
    return res


@app.function(image=image, volumes={"/data": vol}, cpu=8.0, memory=65536,
              timeout=3600)
def baseline():
    """The floor: predict the same profile for every cell. No TCR, no training.

    Fitted on the training split only, then scored on each held-out split with
    the identical NB likelihood, so the numbers are directly comparable.
    """
    import os

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import torch

    vol.reload()
    X = sp.load_npz(f"{DS}/X.npz").tocsr()
    cells = pd.read_parquet(f"{DS}/cells.parquet")
    lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
    lib[lib == 0] = 1.0

    gse, subj = cells["gse"].values, cells["subject"].values
    is_test, is_study = gse == TEST_STUDY, gse == STUDY_HOLDOUT
    pool = ~(is_test | is_study)
    tsub = pd.unique(subj[pool])
    rng = np.random.default_rng(1234)
    val_s = set(rng.choice(tsub, max(1, int(VAL_SUBJECT_FRAC * len(tsub))),
                           replace=False))
    is_val = pool & np.isin(subj, list(val_s))
    is_train = pool & ~is_val

    tot = np.asarray(X[np.where(is_train)[0]].sum(axis=0)).ravel()
    rate = torch.tensor(tot / max(tot.sum(), 1), dtype=torch.float32)

    # Fit one dispersion per gene by moment matching on the training split.
    tr = np.where(is_train)[0]
    sub = X[tr]
    mu_hat = np.outer(lib[tr], rate.numpy())
    var = np.asarray(sub.power(2).mean(axis=0)).ravel() - \
        np.asarray(sub.mean(axis=0)).ravel() ** 2
    mean = np.asarray(sub.mean(axis=0)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        theta = np.where(var > mean, mean ** 2 / (var - mean), 1e4)
    theta = np.clip(np.nan_to_num(theta, nan=1e4), 1e-2, 1e4)
    th = torch.tensor(theta, dtype=torch.float32)

    out = {"arm": "global_mean"}
    for name, m in [("val", is_val), ("study_holdout", is_study),
                    ("test_AS", is_test)]:
        idx = np.where(m)[0]
        s, n = 0.0, 0
        for i in range(0, len(idx), 4096):
            b = idx[i:i + 4096]
            xb = torch.tensor(X[b].toarray(), dtype=torch.float32)
            mu = rate.unsqueeze(0) * torch.tensor(lib[b]).unsqueeze(-1)
            s += nb_nll(xb, mu, th).sum().item()
            n += len(b)
        out[name] = {"nll": s / max(n, 1), "n": n}
        print(f"  {name:14s} nll {out[name]['nll']:10.2f}  n={n:,}")

    os.makedirs(RUNS, exist_ok=True)
    with open(f"{RUNS}/global_mean.json", "w") as fh:
        json.dump(out, fh, indent=1)
    vol.commit()
    return out


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def summary():
    import glob

    vol.reload()
    rows = []
    for p in sorted(glob.glob(f"{RUNS}/*.json")):
        with open(p) as fh:
            r = json.load(fh)
        rows.append(r)
    hdr = (f"{'run':26s} {'val':>10s} {'study':>10s} {'AS test':>10s}"
           f" {'epochs':>7s}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        name = r.get("arm", "?")
        if "beta" in r:
            name += f" b{r['beta']} z{r['z_dim']}"
        print(f"{name:26s} {r.get('val',{}).get('nll',float('nan')):10.2f} "
              f"{r.get('study_holdout',{}).get('nll',float('nan')):10.2f} "
              f"{r.get('test_AS',{}).get('nll',float('nan')):10.2f} "
              f"{r.get('epochs','-'):>7}")

    print("\n" + "=" * 78)
    print("LINEAR PROBE r  (latent -> marker programme score, held-out half)")
    print("=" * 78)
    for r in rows:
        if "beta" not in r:
            continue
        name = f"{r['arm']} b{r['beta']} z{r['z_dim']}"
        for split in ("val", "study_holdout", "test_AS"):
            d = r.get(split, {})
            ps = {k.replace("probe_", "").replace("_r", ""): v
                  for k, v in d.items() if k.startswith("probe_")}
            if ps:
                print(f"  {name:22s} {split:14s} " +
                      "  ".join(f"{k}={v:+.3f}" for k, v in ps.items()))

    print("\n" + "=" * 78)
    print("AS TEST BY CLONE SIZE  (singletons are the honest stratum)")
    print("=" * 78)
    for r in rows:
        if "beta" not in r:
            continue
        name = f"{r['arm']} b{r['beta']} z{r['z_dim']}"
        strata = {k.replace("test_clone_", ""): v["nll"]
                  for k, v in r.items() if k.startswith("test_clone_")}
        if strata:
            print(f"  {name:22s} " +
                  "  ".join(f"{k}={v:.1f}" for k, v in strata.items()))

    print("\n" + "=" * 78)
    print("KL TRAJECTORY  (KL -> 0 means posterior collapse: the latent is")
    print("carrying nothing and the decoder is predicting an average cell)")
    print("=" * 78)
    for r in rows:
        h = r.get("history") or []
        if h:
            name = f"{r['arm']} b{r.get('beta')} z{r.get('z_dim')}"
            kls = [f"{e['kl']:.2f}" for e in h]
            print(f"  {name:22s} {' '.join(kls)}")
            print(f"  {'':22s} total {h[-1]['mins']:.1f} min")

    print("\nLOWER NLL IS BETTER. The number that matters is real vs shuffled:")
    print("if they match, the receptor sequence carries no information about")
    print("the cell's transcriptome and the project stops here.")
    return rows
