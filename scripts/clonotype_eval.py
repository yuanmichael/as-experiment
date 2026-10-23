"""Re-measure the receptor -> state signal at the CLONOTYPE level, by clone size.

Nothing is retrained. The stage-2 checkpoints are reused and only the unit of
evaluation changes, because two things were wrong with measuring per cell.

PSEUDO-REPLICATION. Cells in one expanded clone share a receptor AND a history,
so they are not independent observations. Correlating across cells counts one
biological event hundreds of times. Same error as the KLRB1 result that read
p = 0.0014 per cell and p = 0.289 per patient.

DILUTION. Most T cells are naive -- they have never met their antigen, so the
receptor has not acted and there is no receptor->state relationship to find in
them. 85.8% of our clonotypes are singletons. If the mechanism only operates in
antigen-experienced cells, most of the dataset is noise for this question.

Clone size is never an input to the model, so stratifying by it is not leakage;
it selects the subpopulation where the biology can exist.

Two effects are confounded in "bigger clones score better": averaging more
cells reduces target noise, and expanded clones are antigen-experienced. They
are separated by scoring expanded clones twice -- once against the mean of all
their cells, once against a single randomly chosen cell:

    mean-target minus single-cell-target  =  the averaging benefit
    single-cell-target vs singletons      =  the biology

Programme scores are rebuilt as PC1 of each gene set rather than its mean.
The mean weighted CD160 (loading 0.089) and BTLA (0.075) equally with TOX and
PDCD1, and CTLA4 sat in both the exhaustion and regulatory sets, which made
those two scores correlate at 0.451. PC1 is derived from the data and
downweights the passengers automatically.

    modal run scripts/clonotype_eval.py::evaluate
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

# CTLA4 removed from exhaustion: it is a canonical regulatory-T-cell gene and
# having it in both sets made the two scores correlate at 0.451.
PROGRAMMES = {
    "exhaustion": ["PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "ENTPD1",
                   "CD160", "BTLA"],
    "cytotoxicity": ["GZMA", "GZMB", "GZMH", "PRF1", "GNLY", "NKG7", "IFNG",
                     "KLRG1", "FGFBP2"],
    "naive_memory": ["TCF7", "SELL", "CCR7", "LEF1", "IL7R", "CD27", "CD28"],
    "regulatory": ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "CTLA4"],
}

STRATA = [(1, 1, "singleton"), (2, 2, "2"), (3, 5, "3-5"),
          (6, 20, "6-20"), (21, 10**9, "21+")]


@app.function(image=image, volumes={"/data": vol}, gpu="L4", cpu=8.0,
              memory=65536, timeout=10800)
def evaluate(z_dim: int = 100, d_model: int = 128, n_layers: int = 3,
             min_clones: int = 150):
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import torch
    import torch.nn as nn
    from sklearn.linear_model import Ridge

    vol.reload()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    rng = np.random.default_rng(0)

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
    splits = {"val": is_val, "study_holdout": gse == STUDY_HOLDOUT,
              "test_AS": gse == TEST_STUDY}

    # ---- PC1-weighted programme scores, axes fitted on TRAINING cells only --
    fit_idx = np.sort(rng.choice(np.where(is_train)[0],
                                 min(120000, int(is_train.sum())),
                                 replace=False))
    Yf = X[fit_idx].toarray()
    df_ = Yf.sum(1, keepdims=True)
    df_[df_ == 0] = 1
    Yf = np.log1p(1e4 * Yf / df_)
    axes = {}
    for name, gs in PROGRAMMES.items():
        cols = [gidx[g] for g in gs if g in gidx]
        M = Yf[:, cols]
        C = np.corrcoef(M, rowvar=False)
        w, v = np.linalg.eigh(C)
        load = v[:, int(np.argmax(w))]
        if load.sum() < 0:
            load = -load
        axes[name] = (cols, load / max(np.abs(load).sum(), 1e-9),
                      M.mean(0), M.std(0) + 1e-9)
        print(f"  {name:14s} PC1 share "
              f"{w.max()/w.sum():.3f}  loadings "
              + " ".join(f"{g}:{l:+.2f}" for g, l in zip(
                  [g for g in gs if g in gidx], load)))

    def programme_scores(idx):
        Y = X[idx].toarray()
        d = Y.sum(1, keepdims=True)
        d[d == 0] = 1
        Yn = np.log1p(1e4 * Y / d)
        out = {}
        for name, (cols, load, mu, sd) in axes.items():
            out[name] = ((Yn[:, cols] - mu) / sd) @ load
        return out

    # ------------------------------------------------- tokenise the receptors
    aa = vocab["aa"]
    seqs = np.zeros((len(cells), max_len), dtype=np.int64)
    for i, s in enumerate(cells["cdr3_aa"].astype(str).values):
        for k, c in enumerate(s[:max_len]):
            seqs[i, k] = aa.get(c, 0)
    vg = cells["v_gene"].astype(str).map(vocab["v_gene"]).fillna(0).astype(np.int64).values
    jg = cells["j_gene"].astype(str).map(vocab["j_gene"]).fillna(0).astype(np.int64).values

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

    def predictions(arm):
        import os
        ck = f"{RUNS}/ckpt_{arm}_{z_dim}.pt"
        if not os.path.exists(ck):
            return None
        m = TCREnc().to(dev)
        m.load_state_dict(torch.load(ck, map_location=dev, weights_only=True))
        m.eval()
        out = np.zeros((len(cells), z_dim), dtype=np.float32)
        with torch.no_grad():
            for i in range(0, len(cells), 8192):
                b = np.arange(i, min(i + 8192, len(cells)))
                out[b] = m(torch.tensor(seqs[b], device=dev),
                           torch.tensor(vg[b], device=dev),
                           torch.tensor(jg[b], device=dev)).cpu().numpy()
        return out

    def probe(P, y, groups):
        """Ridge within each group, median r across groups with >= min_clones."""
        rs = []
        for g in pd.unique(groups):
            m = groups == g
            if m.sum() < min_clones:
                continue
            Pg, yg = P[m], y[m]
            h = len(Pg) // 2
            if h < 60 or yg.std() < 1e-9:
                continue
            pred = Ridge(alpha=1.0).fit(Pg[:h], yg[:h]).predict(Pg[h:])
            t = yg[h:]
            if pred.std() > 1e-12 and t.std() > 1e-12:
                rs.append(float(np.corrcoef(pred, t)[0, 1]))
        if not rs:
            return None
        return {"median_r": round(float(np.median(rs)), 4),
                "iqr": [round(float(np.percentile(rs, 25)), 4),
                        round(float(np.percentile(rs, 75)), 4)],
                "n_groups": len(rs)}

    results = {}
    for arm in ("real", "shuffled"):
        P_all = predictions(arm)
        if P_all is None:
            print(f"\n{arm}: no checkpoint, skipped")
            continue
        print(f"\n{'='*86}\n{arm.upper()}\n{'='*86}")
        results[arm] = {}

        for split, mask in splits.items():
            idx = np.where(mask)[0]
            sc = programme_scores(idx)
            sub = cells.iloc[idx]
            clone_key = (sub["gse"].astype(str) + "|"
                         + sub["cdr3_aa"].astype(str) + "|"
                         + sub["v_gene"].astype(str)).values
            lib = sub["library"].values

            frame = pd.DataFrame({"clone": clone_key, "lib": lib,
                                  "row": np.arange(len(idx))})
            for name in PROGRAMMES:
                frame[name] = sc[name]

            # One row per clonotype: mean target, and a single random cell's
            # target, so averaging and biology can be told apart.
            grp = frame.groupby("clone", sort=False)
            size = grp["row"].size()
            mean_rows = grp["row"].first().values          # receptor is shared
            one_rows = grp["row"].apply(
                lambda s: s.iloc[rng.integers(len(s))]).values
            cl_lib = grp["lib"].first().values
            cl_size = size.values
            targets_mean = {n: grp[n].mean().values for n in PROGRAMMES}
            targets_one = {n: frame[n].values[one_rows] for n in PROGRAMMES}
            Pc = P_all[idx][mean_rows]

            results[arm][split] = {"n_clonotypes": int(len(cl_size))}
            print(f"\n  {split}  ({len(cl_size):,} clonotypes, "
                  f"{len(idx):,} cells)")
            for lo, hi, label in STRATA:
                sel = (cl_size >= lo) & (cl_size <= hi)
                if sel.sum() < min_clones:
                    continue
                row = {"n": int(sel.sum())}
                parts = []
                for name in PROGRAMMES:
                    a = probe(Pc[sel], targets_mean[name][sel], cl_lib[sel])
                    if a:
                        row[f"{name}_mean"] = a
                        parts.append(f"{name[:4]}={a['median_r']:+.3f}"
                                     f"[{a['iqr'][0]:+.2f},{a['iqr'][1]:+.2f}]")
                    if lo > 1:
                        b = probe(Pc[sel], targets_one[name][sel], cl_lib[sel])
                        if b:
                            row[f"{name}_single"] = b
                results[arm][split][label] = row
                print(f"    size {label:10s} n={sel.sum():>7,}  "
                      + "  ".join(parts))
                if lo > 1:
                    sp_ = [f"{n[:4]}={row[f'{n}_single']['median_r']:+.3f}"
                           for n in PROGRAMMES if f"{n}_single" in row]
                    if sp_:
                        print(f"      {'(single cell target)':20s} "
                              + "  ".join(sp_))

    with open(f"{RUNS}/clonotype_eval.json", "w") as fh:
        json.dump(results, fh, indent=1)
    vol.commit()

    print("\n" + "=" * 86)
    print("mean-target minus single-cell-target = the averaging benefit")
    print("single-cell-target vs singletons     = the biology")
    print("=" * 86)
    return results


@app.local_entrypoint()
def main():
    evaluate.remote()
