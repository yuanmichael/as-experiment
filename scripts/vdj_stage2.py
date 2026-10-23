"""Stage 2: predict the frozen expression latent from the receptor alone.

    CDR3 + V + J  ->  transformer  ->  predict z_expr (frozen, from stage 1)

Stage 1 is frozen on purpose. It defined the target space without ever seeing a
TCR, so stage 2 cannot quietly redefine what it is being graded on. Trained
jointly, the shared latent would be dominated by the expression branch --
expression predicts expression far better than a receptor does -- and the TCR
branch would wither while still looking load-bearing. That is plapt-ablation's
prot_zero result: a branch worth +0.005 that appeared to be half the model.

What makes this different from the raw-count attempt: the target is ~100 dense
dimensions instead of 2,009 counts that are 81% zeros, so the objective is no
longer dominated by sequencing depth. In the raw-count run the receptor moved
the loss by 2.5% of the achievable range -- invisible -- while a targeted probe
showed a clear real-vs-shuffled gap. Same signal, an objective that can see it.

Three reference points, because a loss value alone means nothing:
    mean-z      predict the training-mean latent for every cell. The floor.
    shuffled    identical model, receptors permuted within subject. The control.
    real        the model under test.

    modal run scripts/vdj_stage2.py::train2 --arm real --z-dim 100
    modal run scripts/vdj_stage2.py::train2 --arm shuffled --z-dim 100
    modal run scripts/vdj_stage2.py::report
"""

import json

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
DS = "/data/vdj/_dataset"
LAT = "/data/vdj/_expr_latent"
RUNS = "/data/vdj/_stage2"

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
              memory=65536, timeout=14400)
def train2(arm: str = "real", z_dim: int = 100, epochs: int = 15,
           batch: int = 1024, lr: float = 1e-3, d_model: int = 128,
           n_layers: int = 3, seed: int = 0, tag: str = ""):
    import os
    import time

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import torch
    import torch.nn as nn
    from sklearn.linear_model import Ridge

    vol.reload()
    torch.manual_seed(seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    Z = np.load(f"{LAT}/z_expr_{z_dim}.npy")
    cells = pd.read_parquet(f"{DS}/cells.parquet")
    X = sp.load_npz(f"{DS}/X.npz").tocsr()
    genes = open(f"{DS}/genes.txt").read().split()
    gidx = {g: i for i, g in enumerate(genes)}
    with open(f"{DS}/vocab.json") as fh:
        vocab = json.load(fh)
    max_len = vocab["max_cdr3"]
    assert len(Z) == len(cells), (Z.shape, len(cells))
    print(f"arm={arm} z_dim={z_dim} | {len(cells):,} cells")

    aa = vocab["aa"]
    seqs = np.zeros((len(cells), max_len), dtype=np.int64)
    for i, s in enumerate(cells["cdr3_aa"].astype(str).values):
        for k, c in enumerate(s[:max_len]):
            seqs[i, k] = aa.get(c, 0)
    vg = cells["v_gene"].astype(str).map(vocab["v_gene"]).fillna(0).astype(np.int64).values
    jg = cells["j_gene"].astype(str).map(vocab["j_gene"]).fillna(0).astype(np.int64).values

    gse, subj = cells["gse"].values, cells["subject"].values
    libs_all = cells["library"].values
    is_test, is_study = gse == TEST_STUDY, gse == STUDY_HOLDOUT
    pool = ~(is_test | is_study)
    tsub = pd.unique(subj[pool])
    rng = np.random.default_rng(1234)
    val_s = set(rng.choice(tsub, max(1, int(VAL_SUBJECT_FRAC * len(tsub))),
                           replace=False))
    is_val = pool & np.isin(subj, list(val_s))
    is_train = pool & ~is_val

    # Standardise the target on the TRAINING split only, so that the reported
    # R^2 is against the training mean and a held-out split cannot flatter the
    # model by happening to be centred differently.
    mu_z = Z[is_train].mean(0, keepdims=True)
    sd_z = Z[is_train].std(0, keepdims=True)
    sd_z[sd_z < 1e-6] = 1.0
    Zs = (Z - mu_z) / sd_z

    if arm == "shuffled":
        order = np.arange(len(cells))
        for s in pd.unique(subj):
            idx = np.where(subj == s)[0]
            perm = idx.copy()
            rng.shuffle(perm)
            order[idx] = perm
        seqs, vg, jg = seqs[order], vg[order], jg[order]
        print("  receptors permuted within subject")

    class TCREnc(nn.Module):
        def __init__(self):
            super().__init__()
            self.aa = nn.Embedding(21, d_model, padding_idx=0)
            self.pos = nn.Embedding(max_len, d_model)
            layer = nn.TransformerEncoderLayer(
                d_model, 4, d_model * 2, dropout=0.1, batch_first=True,
                norm_first=True)
            self.tr = nn.TransformerEncoder(layer, n_layers)
            self.v = nn.Embedding(len(vocab["v_gene"]) + 1, 32, padding_idx=0)
            self.j = nn.Embedding(len(vocab["j_gene"]) + 1, 16, padding_idx=0)
            self.head = nn.Sequential(
                nn.Linear(d_model + 48, 256), nn.LayerNorm(256), nn.GELU(),
                nn.Linear(256, 256), nn.LayerNorm(256), nn.GELU(),
                nn.Linear(256, z_dim))

        def forward(self, seq, v, j):
            pad = seq == 0
            pos = torch.arange(seq.shape[1], device=seq.device)
            h = self.tr(self.aa(seq) + self.pos(pos)[None],
                        src_key_padding_mask=pad)
            keep = (~pad).float().unsqueeze(-1)
            h = (h * keep).sum(1) / keep.sum(1).clamp(min=1)
            return self.head(torch.cat([h, self.v(v), self.j(j)], -1))

    model = TCREnc().to(dev)
    print(f"  {sum(p.numel() for p in model.parameters()):,} parameters")
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    tr_idx = np.where(is_train)[0]

    def predict(mask):
        model.eval()
        idx = np.where(mask)[0]
        out = []
        with torch.no_grad():
            for i in range(0, len(idx), 4096):
                b = idx[i:i + 4096]
                out.append(model(torch.tensor(seqs[b], device=dev),
                                 torch.tensor(vg[b], device=dev),
                                 torch.tensor(jg[b], device=dev)).cpu().numpy())
        return idx, np.concatenate(out) if out else np.zeros((0, z_dim))

    def scores(mask, probe=False):
        idx, P = predict(mask)
        T = Zs[idx]
        mse = float(((P - T) ** 2).mean())
        # R^2 against predicting the training mean (which is 0 after
        # standardising). Negative means worse than the constant predictor.
        r2 = float(1 - ((P - T) ** 2).sum() / max((T ** 2).sum(), 1e-9))
        # Correlation per latent dimension, then its median: a fairer summary
        # than overall MSE when some dimensions are intrinsically noisier.
        rs = []
        for k in range(z_dim):
            if P[:, k].std() > 1e-9 and T[:, k].std() > 1e-9:
                rs.append(float(np.corrcoef(P[:, k], T[:, k])[0, 1]))
        out = {"mse": round(mse, 5), "r2": round(r2, 5), "n": len(idx),
               "median_dim_r": round(float(np.median(rs)), 4) if rs else None,
               "max_dim_r": round(float(np.max(rs)), 4) if rs else None}

        if probe:
            # The same measurement the raw-count model reported, so the two are
            # directly comparable: can a linear read of the PREDICTED latent
            # recover a marker programme score, within a single library?
            Y = X[idx].toarray()
            d = Y.sum(1, keepdims=True)
            d[d == 0] = 1
            Yn = np.log1p(1e4 * Y / d)
            L = libs_all[idx]
            big = [l for l, c in pd.Series(L).value_counts().items() if c >= 400]

            def probe_r(Zp, y):
                h = len(Zp) // 2
                if h < 50:
                    return None
                p = Ridge(alpha=1.0).fit(Zp[:h], y[:h]).predict(Zp[h:])
                t = y[h:]
                if p.std() == 0 or t.std() == 0:
                    return None
                return float(np.corrcoef(p, t)[0, 1])

            for name, gs in PROGRAMMES.items():
                cols = [gidx[g] for g in gs if g in gidx]
                y = Yn[:, cols].mean(1)
                r = probe_r(P, y)
                if r is not None:
                    out[f"probe_{name}_r"] = round(r, 4)
                w = [probe_r(P[L == l], y[L == l]) for l in big]
                w = [x for x in w if x is not None]
                if w:
                    out[f"probe_{name}_within_lib_r"] = round(
                        float(np.median(w)), 4)
                    # Spread matters: a median over 13 libraries can bounce.
                    out[f"probe_{name}_within_lib_iqr"] = [
                        round(float(np.percentile(w, 25)), 4),
                        round(float(np.percentile(w, 75)), 4)]
                    out[f"probe_{name}_n_libs"] = len(w)
        return out

    hist, best, t0 = [], None, time.time()
    for ep in range(epochs):
        model.train()
        rng.shuffle(tr_idx)
        run, nb = 0.0, 0
        for i in range(0, len(tr_idx), batch):
            b = tr_idx[i:i + batch]
            pred = model(torch.tensor(seqs[b], device=dev),
                         torch.tensor(vg[b], device=dev),
                         torch.tensor(jg[b], device=dev))
            loss = ((pred - torch.tensor(Zs[b], device=dev)) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            run += loss.item()
            nb += 1
        v = scores(is_val)
        hist.append({"epoch": ep, "train_mse": run / max(nb, 1),
                     "val_mse": v["mse"], "val_r2": v["r2"],
                     "mins": round((time.time() - t0) / 60, 1)})
        print(f"  ep{ep:02d} train {run/max(nb,1):7.4f}  val_mse {v['mse']:7.4f}"
              f"  val_r2 {v['r2']:+7.4f}  {hist[-1]['mins']:5.1f}m")
        if best is None or v["mse"] < best:
            best = v["mse"]
            os.makedirs(RUNS, exist_ok=True)
            torch.save(model.state_dict(), f"{RUNS}/ckpt_{arm}_{z_dim}.pt")
            vol.commit()

    model.load_state_dict(torch.load(f"{RUNS}/ckpt_{arm}_{z_dim}.pt",
                                     weights_only=True))

    res = {"arm": arm, "z_dim": z_dim, "epochs": epochs, "history": hist,
           "n_params": sum(p.numel() for p in model.parameters())}
    # Save before diagnostics: the weights are the expensive artefact.
    name = tag or f"{arm}_z{z_dim}_s{seed}"
    with open(f"{RUNS}/{name}.json", "w") as fh:
        json.dump(res, fh, indent=1)
    vol.commit()

    for split, m in [("val", is_val), ("study_holdout", is_study),
                     ("test_AS", is_test)]:
        try:
            res[split] = scores(m, probe=True)
        except Exception as exc:
            res[split] = {"error": f"{type(exc).__name__}: {exc}"}
            print(f"  {split} scoring failed: {res[split]['error']}")

    # The floor: predict the training mean for every cell.
    for split, m in [("val", is_val), ("study_holdout", is_study),
                     ("test_AS", is_test)]:
        T = Zs[np.where(m)[0]]
        res[f"floor_{split}_mse"] = round(float((T ** 2).mean()), 5)

    with open(f"{RUNS}/{name}.json", "w") as fh:
        json.dump(res, fh, indent=1)
    vol.commit()
    print("\n" + json.dumps({k: v for k, v in res.items()
                             if k != "history"}, indent=1))
    return res


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def report():
    import glob
    vol.reload()
    rows = []
    for p in sorted(glob.glob(f"{RUNS}/*.json")):
        with open(p) as fh:
            r = json.load(fh)
        if "arm" in r:
            rows.append(r)

    hdr = (f"{'run':22s} {'val r2':>8s} {'val medr':>9s} "
           f"{'study r2':>9s} {'AS r2':>8s} {'AS medr':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda d: (d["z_dim"], d["arm"])):
        n = f"{r['arm']} z{r['z_dim']}"
        print(f"{n:22s} {r.get('val',{}).get('r2',float('nan')):8.4f} "
              f"{r.get('val',{}).get('median_dim_r',float('nan')) or float('nan'):9.4f} "
              f"{r.get('study_holdout',{}).get('r2',float('nan')):9.4f} "
              f"{r.get('test_AS',{}).get('r2',float('nan')):8.4f} "
              f"{r.get('test_AS',{}).get('median_dim_r',float('nan')) or float('nan'):8.4f}")

    print("\n" + "=" * 92)
    print("WITHIN-LIBRARY PROGRAMME PROBE  (median [IQR] across libraries)")
    print("directly comparable to the raw-count model: real 0.244 / shuffled")
    print("0.007 for exhaustion on the AS test")
    print("=" * 92)
    for r in sorted(rows, key=lambda d: (d["z_dim"], d["arm"])):
        n = f"{r['arm']} z{r['z_dim']}"
        for split in ("val", "study_holdout", "test_AS"):
            d = r.get(split, {})
            parts = []
            for prog in PROGRAMMES:
                k = f"probe_{prog}_within_lib_r"
                if k in d:
                    iqr = d.get(f"probe_{prog}_within_lib_iqr", [None, None])
                    parts.append(f"{prog[:4]}={d[k]:+.3f}[{iqr[0]:+.2f},{iqr[1]:+.2f}]")
            if parts:
                print(f"  {n:18s} {split:14s} " + "  ".join(parts))
    print("\nr2 is against predicting the training-mean latent. Negative means")
    print("worse than that constant. real vs shuffled is the only comparison")
    print("that establishes the receptor carries information.")
    return rows
