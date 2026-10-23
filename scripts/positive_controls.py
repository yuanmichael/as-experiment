"""Positive controls: can this pipeline find signal that is definitely there?

Every result so far is a negative measured against a shuffled baseline. That
establishes "not much signal". It cannot establish that the machinery is
capable of finding signal, because a bug in tokenisation, the splits, the probe
or the latent would produce identical numbers. Until something known to be
receptor-encoded comes back positive, every negative is uninterpretable.

Three targets, in order of how certain they are:

  1. PLUMBING -- CDR3 length, and V gene identity. These are the model's own
     INPUT. Recovery must be near-perfect. If the predicted latent cannot
     recover the length of the sequence it was computed from, the encoder or
     the tokenisation is broken and nothing else means anything.

  2. CD4 vs CD8 -- a real, published receptor->phenotype relationship:
     thymic selection on different MHC classes leaves different V-gene usage
     and CDR3 length distributions. Expect modest, NOT near-perfect, because
     we hold only the beta chain; beta-only classifiers in the literature sit
     around AUC 0.6-0.7.

  3. MAIT identity -- the weakest of the three for the same reason. MAIT cells
     are defined by a near-invariant ALPHA chain (TRAV1-2), which we do not
     have. Their beta chains show biased TRBV20/TRBV6 usage and little else,
     so a modest result here is expected and a null is not damning.

Each target is scored two ways:

  MODEL     the stage-2 predicted latent (receptor -> z), probed for the label.
  DIRECT    plain logistic regression on hand-made receptor features: V and J
            one-hot, CDR3 length, amino acid composition.

DIRECT is the load-bearing comparison. It asks whether the information exists
in the receptor at all, independent of anything we built. If DIRECT clearly
beats MODEL, the information is there and our model is failing to use it --
which would mean the whole chain of negatives was about our machinery, not the
biology.

Labels come from expression, so they are noisy but independent of the receptor.

    modal run scripts/positive_controls.py::controls
"""

import json

import modal

app = modal.App("as-experiment")
vol = modal.Volume.from_name("as-data", create_if_missing=True)
DS = "/data/vdj/_dataset"
RUNS = "/data/vdj/_stage2"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("torch==2.5.1", "pandas==2.2.3", "numpy==2.1.3",
                 "scipy==1.14.1", "pyarrow==18.1.0", "scikit-learn==1.5.2")
)

TEST_STUDY = "GSE288581"
STUDY_HOLDOUT = "GSE329527"
VAL_SUBJECT_FRAC = 0.15
AA = "ACDEFGHIKLMNPQRSTVWY"


@app.function(image=image, volumes={"/data": vol}, gpu="L4", cpu=8.0,
              memory=65536, timeout=10800)
def controls(z_dim: int = 100, d_model: int = 128, n_layers: int = 3,
             n_cells: int = 200000):
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import torch
    import torch.nn as nn
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.metrics import roc_auc_score

    vol.reload()
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    X = sp.load_npz(f"{DS}/X.npz").tocsr()
    cells = pd.read_parquet(f"{DS}/cells.parquet")
    genes = open(f"{DS}/genes.txt").read().split()
    gidx = {g: i for i, g in enumerate(genes)}
    with open(f"{DS}/vocab.json") as fh:
        vocab = json.load(fh)
    max_len = vocab["max_cdr3"]

    gse, subj = cells["gse"].values, cells["subject"].values
    pool = ~((gse == TEST_STUDY) | (gse == STUDY_HOLDOUT))
    tsub = pd.unique(subj[pool])
    vr = np.random.default_rng(1234)
    val_s = set(vr.choice(tsub, max(1, int(VAL_SUBJECT_FRAC * len(tsub))),
                          replace=False))
    is_val = pool & np.isin(subj, list(val_s))
    is_train = pool & ~is_val

    # ------------------------------------------------------------- receptors
    aa_map = vocab["aa"]
    cdr3 = cells["cdr3_aa"].astype(str).values
    seqs = np.zeros((len(cells), max_len), dtype=np.int64)
    for i, s in enumerate(cdr3):
        for k, c in enumerate(s[:max_len]):
            seqs[i, k] = aa_map.get(c, 0)
    vg = cells["v_gene"].astype(str).map(vocab["v_gene"]).fillna(0).astype(np.int64).values
    jg = cells["j_gene"].astype(str).map(vocab["j_gene"]).fillna(0).astype(np.int64).values
    cdr3_len = np.array([len(s) for s in cdr3], dtype=np.float32)

    # DIRECT features: nothing learned, just the receptor written out.
    nV, nJ = len(vocab["v_gene"]) + 1, len(vocab["j_gene"]) + 1
    comp = np.zeros((len(cells), len(AA)), dtype=np.float32)
    for i, s in enumerate(cdr3):
        if not s:
            continue
        for c in s:
            k = AA.find(c)
            if k >= 0:
                comp[i, k] += 1
        comp[i] /= max(len(s), 1)
    Voh = np.zeros((len(cells), nV), dtype=np.float32)
    Voh[np.arange(len(cells)), vg] = 1
    Joh = np.zeros((len(cells), nJ), dtype=np.float32)
    Joh[np.arange(len(cells)), jg] = 1
    DIRECT = np.hstack([Voh, Joh, cdr3_len[:, None] / 30.0, comp])
    print(f"DIRECT features: {DIRECT.shape[1]} "
          f"(V one-hot {nV}, J one-hot {nJ}, length, 20 aa composition)")

    # ------------------------------------------------------ model prediction
    class TCREnc(nn.Module):
        def __init__(self):
            super().__init__()
            self.aa = nn.Embedding(21, d_model, padding_idx=0)
            self.pos = nn.Embedding(max_len, d_model)
            layer = nn.TransformerEncoderLayer(
                d_model, 4, d_model * 2, dropout=0.1, batch_first=True,
                norm_first=True)
            self.tr = nn.TransformerEncoder(layer, n_layers)
            self.v = nn.Embedding(nV, 32, padding_idx=0)
            self.j = nn.Embedding(nJ, 16, padding_idx=0)
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

    m = TCREnc().to(dev)
    m.load_state_dict(torch.load(f"{RUNS}/ckpt_real_{z_dim}.pt",
                                 map_location=dev, weights_only=True))
    m.eval()
    MODEL = np.zeros((len(cells), z_dim), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, len(cells), 8192):
            b = np.arange(i, min(i + 8192, len(cells)))
            MODEL[b] = m(torch.tensor(seqs[b], device=dev),
                         torch.tensor(vg[b], device=dev),
                         torch.tensor(jg[b], device=dev)).cpu().numpy()

    # --------------------------------------------- labels, from expression
    def col(g):
        return X[:, gidx[g]].toarray().ravel() if g in gidx else None

    cd4, cd8a, cd8b = col("CD4"), col("CD8A"), col("CD8B")
    cd8 = (cd8a if cd8a is not None else 0) + (cd8b if cd8b is not None else 0)
    # Unambiguous cells only: one lineage marker present, the other absent.
    lab_cd48 = np.full(len(cells), -1, dtype=np.int64)
    lab_cd48[(cd8 > 0) & (cd4 == 0)] = 1
    lab_cd48[(cd4 > 0) & (cd8 == 0)] = 0

    mait_g = col("SLC4A10")
    lab_mait = (mait_g > 0).astype(np.int64) if mait_g is not None else None

    rng = np.random.default_rng(0)

    def auc(feats, label, train_mask, test_mask, tag):
        tr = np.where(train_mask & (label >= 0))[0]
        te = np.where(test_mask & (label >= 0))[0]
        if len(tr) > n_cells:
            tr = rng.choice(tr, n_cells, replace=False)
        if len(te) > 60000:
            te = rng.choice(te, 60000, replace=False)
        if len(np.unique(label[tr])) < 2 or len(np.unique(label[te])) < 2:
            return None
        clf = LogisticRegression(max_iter=400, C=1.0).fit(feats[tr], label[tr])
        p = clf.predict_proba(feats[te])[:, 1]
        return {"auc": round(float(roc_auc_score(label[te], p)), 4),
                "n_train": len(tr), "n_test": len(te),
                "pct_positive": round(float(label[te].mean()), 3)}

    def r2(feats, target, train_mask, test_mask):
        tr = np.where(train_mask)[0]
        te = np.where(test_mask)[0]
        if len(tr) > n_cells:
            tr = rng.choice(tr, n_cells, replace=False)
        if len(te) > 60000:
            te = rng.choice(te, 60000, replace=False)
        pred = Ridge(alpha=1.0).fit(feats[tr], target[tr]).predict(feats[te])
        t = target[te]
        if pred.std() < 1e-12 or t.std() < 1e-12:
            return None
        return round(float(np.corrcoef(pred, t)[0, 1]), 4)

    out = {}
    splits = {"val": is_val, "study_holdout": gse == STUDY_HOLDOUT,
              "test_AS": gse == TEST_STUDY}

    print("\n" + "=" * 78)
    print("1. PLUMBING -- recovering the model's OWN INPUT from its output")
    print("   must be near 1.0; anything else means the encoder is broken")
    print("=" * 78)
    for split, mask in splits.items():
        rl = r2(MODEL, cdr3_len, is_train, mask)
        out[f"plumbing_cdr3_len_{split}"] = rl
        print(f"  {split:14s} CDR3 length from latent, r = {rl}")
    # V gene is categorical: can the latent tell the two commonest V genes apart?
    vc = pd.Series(vg).value_counts()
    if len(vc) >= 2:
        a, b = vc.index[0], vc.index[1]
        lab_v = np.full(len(cells), -1, dtype=np.int64)
        lab_v[vg == a] = 1
        lab_v[vg == b] = 0
        for split, mask in splits.items():
            s = auc(MODEL, lab_v, is_train, mask, "v")
            out[f"plumbing_vgene_{split}"] = s
            if s:
                print(f"  {split:14s} V gene {a} vs {b}, AUC = {s['auc']}")

    print("\n" + "=" * 78)
    print("2. CD4 vs CD8 -- real published relationship; beta chain only, so")
    print("   expect roughly 0.6-0.7, not 0.95")
    print("=" * 78)
    print(f"  labelled cells: {int((lab_cd48 >= 0).sum()):,} of {len(cells):,} "
          f"({100*(lab_cd48 >= 0).mean():.1f}%)")
    for split, mask in splits.items():
        mo = auc(MODEL, lab_cd48, is_train, mask, "cd48")
        di = auc(DIRECT, lab_cd48, is_train, mask, "cd48")
        out[f"cd4cd8_model_{split}"], out[f"cd4cd8_direct_{split}"] = mo, di
        if mo and di:
            print(f"  {split:14s} MODEL {mo['auc']:.4f}   DIRECT {di['auc']:.4f}"
                  f"   (n_test {mo['n_test']:,}, "
                  f"{100*mo['pct_positive']:.0f}% CD8)")

    if lab_mait is not None:
        print("\n" + "=" * 78)
        print("3. MAIT (SLC4A10+) -- weakest control: MAIT identity is carried")
        print("   by the ALPHA chain, which this dataset does not contain")
        print("=" * 78)
        print(f"  MAIT cells: {int(lab_mait.sum()):,} "
              f"({100*lab_mait.mean():.2f}%)")
        for split, mask in splits.items():
            mo = auc(MODEL, lab_mait, is_train, mask, "mait")
            di = auc(DIRECT, lab_mait, is_train, mask, "mait")
            out[f"mait_model_{split}"], out[f"mait_direct_{split}"] = mo, di
            if mo and di:
                print(f"  {split:14s} MODEL {mo['auc']:.4f}   "
                      f"DIRECT {di['auc']:.4f}")

    print("\n" + "=" * 78)
    print("READING THIS")
    print("=" * 78)
    print("  plumbing near 1.0, CD4/CD8 DIRECT ~0.6-0.7, MODEL similar")
    print("     -> the pipeline works; the weak programme signal is real biology")
    print("  DIRECT clearly above MODEL")
    print("     -> the information is in the receptor and our model is not")
    print("        using it. The negatives were about our machinery.")
    print("  plumbing well below 1.0")
    print("     -> encoder or tokenisation is broken; re-measure everything")
    print("  CD4/CD8 at 0.5 for BOTH")
    print("     -> labels are wrong, or the beta chain really carries nothing")

    with open(f"{DS}/positive_controls.json", "w") as fh:
        json.dump(out, fh, indent=1)
    vol.commit()
    return out


@app.local_entrypoint()
def main():
    controls.remote()
