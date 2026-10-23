"""Embed CDR3s with ESM-2 and score them against the AS clonotype.

Why this exists. Edit distance says `CASSVGLASTDTQYF` is one step from
`CASSVGLYSTDTQYF`. Biologically they are not close at all: the first is carried
by 147 of 786 healthy people, the second by 15, and only the second tracks
HLA-B27 and AS. Any family built by string distance is therefore useless here --
within one edit, 74% of healthy subjects are "carriers".

The question is whether a protein language model's representation separates
those two, i.e. whether it encodes something beyond spelling. Mal-ID
(Zaslavsky et al., Science 2025) reached 0.986 AUROC on five diseases using
off-the-shelf ESM-2 150M embeddings of CDR3s, so this is not a long shot.

Only one number per sequence is kept -- its cosine similarity to the target --
so 4M sequences cost 16 MB instead of 10 GB of 640-dim vectors.

    modal run scripts/esm_embed.py::embed_pool
    modal run scripts/esm_embed.py::report
"""
import modal

app = modal.App("as-embed")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

MODEL = "facebook/esm2_t30_150M_UR50D"      # 640-dim; the size Mal-ID used
TARGET = "CASSVGLYSTDTQYF"
POOL = "/data/embed/embed_pool.parquet"
OUT = "/data/embed/esm2_cosine.parquet"

# Reference points for the sanity check, with their healthy carriage in the 786
# Emerson repertoires. If the model has learned spelling, the first three rank
# as the target's nearest neighbours despite being common public sequences.
PROBES = {
    "CASSVGLYSTDTQYF": "the target, 15/786 healthy, B27-linked",
    "CASSVGLASTDTQYF": "1 edit, 147/786 healthy",
    "CASSVGLSTDTQYF": "1 edit, 114/786 healthy",
    "CASSLGLSTDTQYF": "2 edits, 634/786 healthy",
    "CASSQETQYF": "unrelated CDR3",
    "CASSIRSSYEQYF": "unrelated CDR3",
}

image = (
    modal.Image.debian_slim(python_version="3.12")
    # Pinned: an unpinned sklearn once killed a finished run at its last step.
    .pip_install("torch==2.5.1", "transformers==4.46.3", "pandas==2.2.3",
                 "pyarrow==17.0.0", "numpy==2.1.3")
    # Bake the weights into the image so the GPU container does not spend
    # billable minutes downloading 600 MB.
    .run_commands(
        f"python -c \"from transformers import AutoModel, AutoTokenizer; "
        f"AutoTokenizer.from_pretrained('{MODEL}'); "
        f"AutoModel.from_pretrained('{MODEL}')\""
    )
)


def _load(device):
    import torch
    from transformers import AutoModel, AutoTokenizer

    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModel.from_pretrained(MODEL).to(device).eval()
    if device == "cuda":
        model = model.half()
    return tok, model, torch


def _embed(seqs, tok, model, torch, device):
    """Mean-pool over residues only: special tokens would add a constant that
    shrinks the angular spread between short sequences."""
    enc = tok(list(seqs), return_tensors="pt", padding=True, truncation=True,
              max_length=48)
    enc = {k: v.to(device) for k, v in enc.items()}
    with torch.no_grad():
        h = model(**enc).last_hidden_state
    mask = enc["attention_mask"].clone()
    mask[:, 0] = 0                                   # <cls>
    lengths = mask.sum(1, keepdim=True)
    idx = lengths.squeeze(1)                          # last real token is <eos>
    mask[torch.arange(mask.shape[0], device=device), idx] = 0
    m = mask.unsqueeze(-1).to(h.dtype)
    return (h * m).sum(1) / m.sum(1).clamp(min=1)


@app.function(image=image, volumes={"/data": vol}, gpu="A10G", timeout=5400)
def embed_pool(batch_size: int = 2048):
    import os
    import time

    import numpy as np
    import pandas as pd
    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"
    tok, model, _ = _load(device)
    print(f"{MODEL} on {device}", flush=True)

    tgt = _embed([TARGET], tok, model, torch, device)
    tgt = torch.nn.functional.normalize(tgt.float(), dim=1)

    print("probe cosines (before the sweep, so they cannot be cherry-picked):")
    pv = _embed(list(PROBES), tok, model, torch, device)
    pv = torch.nn.functional.normalize(pv.float(), dim=1)
    for s, c in zip(PROBES, (pv @ tgt.T).squeeze(1).tolist()):
        print(f"  {c:+.4f}  {s:<16} {PROBES[s]}", flush=True)

    seqs = pd.read_parquet(POOL)["cdr3"].astype(str).to_numpy()
    print(f"\nembedding {len(seqs):,} CDR3s", flush=True)
    cos = np.empty(len(seqs), dtype=np.float32)
    t0 = time.time()
    for i in range(0, len(seqs), batch_size):
        chunk = seqs[i:i + batch_size]
        v = _embed(chunk, tok, model, torch, device)
        v = torch.nn.functional.normalize(v.float(), dim=1)
        cos[i:i + len(chunk)] = (v @ tgt.T).squeeze(1).cpu().numpy()
        if (i // batch_size) % 200 == 0:
            done = i + len(chunk)
            rate = done / max(time.time() - t0, 1e-9)
            print(f"  {done:,}/{len(seqs):,}  {rate:,.0f} seq/s  "
                  f"eta {(len(seqs) - done) / max(rate, 1e-9) / 60:.1f} min",
                  flush=True)

    os.makedirs("/data/embed", exist_ok=True)
    pd.DataFrame({"cdr3": seqs, "cos": cos}).to_parquet(OUT, index=False)
    vol.commit()
    print(f"\nwrote {OUT} in {(time.time() - t0) / 60:.1f} min")
    return {"n": int(len(seqs)), "minutes": round((time.time() - t0) / 60, 1)}


@app.function(image=image, volumes={"/data": vol}, timeout=1800)
def report(k: int = 50):
    """Nearest neighbours of the target in embedding space."""
    import pandas as pd

    df = pd.read_parquet(OUT).sort_values("cos", ascending=False)
    print(f"{len(df):,} sequences scored\n")
    print(f"top {k} by cosine to {TARGET}:")
    print(df.head(k).to_string(index=False))
    print("\nwhere the string neighbours rank:")
    rank = {c: i for i, c in enumerate(df["cdr3"], 1)}
    for s, note in PROBES.items():
        r = rank.get(s)
        c = df.loc[df["cdr3"] == s, "cos"]
        print(f"  rank {r if r else '-':>9}  cos {float(c.iloc[0]) if len(c) else float('nan'):+.4f}"
              f"  {s:<16} {note}")
    return df.head(k)["cdr3"].tolist()
