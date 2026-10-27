"""Score arbitrary CDR3-beta sequences against a peptide set with the PRP joint model.

Why this exists. `prp_divergence.py` establishes, on Wang/Yeh/Garcia's own 21
cached profiles, that CDR3-beta edit distance carries no information about what
a TCR recognises (motif JS vs edit distance, Spearman +0.113, p = 0.10, flat
from 0 to 7 edits). The functional relatives of our clonotype `CASSVGLYSTDTQYF`
sit *three* edits away, and its one-edit neighbours are not in its top six. That
explains why both of this project's widening attempts collapsed at radius 1.

This asks the next question: among the clonotypes our own AS screens surfaced,
which are functionally near the target -- and are there any we should have tried
and did not?

The model takes a (CDR3b, peptide) pair and returns a binding logit, and
`inference_proteome.py` stamps any CDR3b onto every panel row, so novel betas
can be scored. The catch is that each *single-TCR* checkpoint was fine-tuned
with its beta held constant, so it never saw beta vary and its response to a
novel beta is unconstrained. The `19.2_allmutants` joint checkpoint is the right
tool instead: it was trained on 19.2 plus five engineered CDR3b mutants --
six distinct betas, one of which is our exact target -- and the paper reports
that joint modelling inside a recognition-defined neighbourhood improves
transfer.

Two measurements per candidate beta, on one fixed peptide set shared by every
beta so the numbers are comparable:

- **target-set AUC**: can this beta's scores separate the *target's* recognised
  peptides from background? High means it recognises what the target does. This
  is the cheap, direct statistic -- it needs only the target's top peptides plus
  a matched background, not the whole 293k proteome.
- **motif JS** to the target, as the paper defines it (per-position
  Jensen-Shannon distance, P2/P8/P9 excluded), for comparability with Fig. 6c.

The positive control runs first and the discovery run does not start until it
passes: score the 15 panel betas the joint model never trained on, and check
that their target-set AUC tracks the motif JS computed independently from the
cached single-TCR profiles. If a model cannot recover the ordering of TCRs whose
answer we already know, its ranking of our candidates means nothing -- a screen
that cannot recover a known signal is uninformative, not negative.

Leakage note. The 167-sequence panel was selected using AS/healthy labels, so
ranking it by functional distance is **descriptive**, not a test. No p-value
about AS association can come out of this list; any such claim needs an
evaluation set that selection never touched.

    modal run scripts/prp_score.py::smoke            # measure throughput first
    modal run scripts/prp_score.py::control
    modal run scripts/prp_score.py::score_candidates
"""
import os

import modal

app = modal.App("as-prp")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

CKPT_URL = ("https://huggingface.co/hyeh/PRP-TCR-Specificity/resolve/main/"
            "lightning_logs/19.2_allmutants_finetune_esm.zip")
CKPT_MEMBER = "19.2_allmutants_finetune_esm/checkpoints/best-checkpoint.ckpt"
CKPT_LOCAL = "/data/prp/19.2_allmutants.ckpt"
ESM_PT = "/opt/esm/esm2_t33_650M_UR50D.pt"
TARGET = "CASSVGLYSTDTQYF"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "curl")
    # Pinned: an unpinned sklearn once killed a finished run at its last step.
    .pip_install(
        "torch==2.3.1", "fair-esm==2.0.0", "pytorch-lightning==2.3.3",
        "einops==0.8.0", "easydict==1.13", "numpy==1.26.4",
        "pandas==2.2.2", "pyyaml==6.0.1", "tqdm==4.66.4",
        "scikit-learn==1.5.1", "scipy==1.14.1",
    )
    # Bake the base ESM weights and the source tree so GPU containers do not
    # spend billable minutes fetching 2.5 GB. fair-esm looks for the
    # contact-regression sibling next to the model file, so take both.
    .run_commands(
        "mkdir -p /opt/esm",
        f"curl -sL -o {ESM_PT} https://dl.fbaipublicfiles.com/fair-esm/models/"
        "esm2_t33_650M_UR50D.pt",
        "curl -sL -o /opt/esm/esm2_t33_650M_UR50D-contact-regression.pt "
        "https://dl.fbaipublicfiles.com/fair-esm/regression/"
        "esm2_t33_650M_UR50D-contact-regression.pt",
        "git clone --depth 1 https://github.com/akds/PRP-TCR-Specificity.git /opt/prp",
    )
)


# --------------------------------------------------------------------------
# checkpoint: stream one member out of a 3.2 GB remote zip


@app.function(image=image, volumes={"/data": vol}, timeout=3600)
def fetch_ckpt(force: bool = False):
    import io
    import os
    import urllib.request
    import zipfile

    if os.path.exists(CKPT_LOCAL) and not force:
        return {"ckpt": CKPT_LOCAL, "bytes": os.path.getsize(CKPT_LOCAL),
                "cached": True}

    class HttpFile(io.RawIOBase):
        def __init__(self, url):
            self.url, self._pos = url, 0
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=60) as r:
                self.size = int(r.headers["Content-Length"])

        def seek(self, off, whence=io.SEEK_SET):
            self._pos = (off if whence == io.SEEK_SET
                         else self._pos + off if whence == io.SEEK_CUR
                         else self.size + off)
            return self._pos

        def tell(self):
            return self._pos

        def seekable(self):
            return True

        def readable(self):
            return True

        def read(self, n=-1):
            if n is None or n < 0:
                n = self.size - self._pos
            if n == 0 or self._pos >= self.size:
                return b""
            end = min(self._pos + n, self.size) - 1
            req = urllib.request.Request(
                self.url, headers={"Range": f"bytes={self._pos}-{end}"})
            with urllib.request.urlopen(req, timeout=600) as r:
                data = r.read()
            self._pos += len(data)
            return data

    os.makedirs("/data/prp", exist_ok=True)
    with zipfile.ZipFile(HttpFile(CKPT_URL)) as z:
        with z.open(CKPT_MEMBER) as src, open(CKPT_LOCAL, "wb") as out:
            while chunk := src.read(1 << 24):
                out.write(chunk)
    vol.commit()
    return {"ckpt": CKPT_LOCAL, "bytes": os.path.getsize(CKPT_LOCAL),
            "cached": False}


# --------------------------------------------------------------------------
# model + scoring


def _load_model(device):
    import sys

    import torch
    from easydict import EasyDict

    sys.path.insert(0, "/opt/prp/source")
    from TCRPNet import TCRPNet
    from trainer_esm import DiscriminativePL

    args = EasyDict({
        "esm_weights": ESM_PT,
        "freeze_encoder": False,
        "use_pos_weight": False,
        "seed": 1,
        "device": device,
    })
    torch.manual_seed(args.seed)
    # The checkpoint was trained with use_pos_weight: true, so BCEWithLogitsLoss
    # registered a `pos_weight` buffer. Constructing with pos_weight=None omits
    # that buffer and strict loading then rejects the key. The value is
    # irrelevant at inference -- only the buffer's existence matters.
    model = DiscriminativePL(args=args, model=TCRPNet(),
                             pos_weight=torch.tensor(1.0))
    ckpt = torch.load(CKPT_LOCAL, map_location="cpu")
    model.load_state_dict(ckpt["state_dict"], strict=True)
    return model.to(device).eval(), args


def _score(model, args, cdr3s, peptides, batch_size, device, log_every=20):
    """Returns a (len(cdr3s), len(peptides)) float32 array of logits."""
    import sys
    import time

    import numpy as np
    import pandas as pd
    import torch
    from torch.utils.data import DataLoader

    sys.path.insert(0, "/opt/prp/source")
    from data import BatchConverterCollater, DeepSequencingRawSequenceDataset

    collate = BatchConverterCollater(alphabet=model.alphabet)
    out = np.empty((len(cdr3s), len(peptides)), dtype=np.float32)
    t0 = time.time()
    done = 0
    for i, cdr in enumerate(cdr3s):
        df = pd.DataFrame({"Epitope": peptides})
        df["CDR3_b"] = cdr
        df["Score"] = 1
        df["Split"] = "test"
        ds = DeepSequencingRawSequenceDataset(None, "test", data=df)
        dl = DataLoader(ds, batch_size=batch_size, shuffle=False,
                        collate_fn=collate)
        preds = []
        with torch.no_grad():
            for b in dl:
                y = model(
                    TCRb=b["cdr"].to(device),
                    epitope=b["peptide"].to(device),
                    TCRb_ohe=b["cdr_ohe"].to(device),
                    epitope_ohe=b["peptide_ohe"].to(device),
                    TCRb_mask=b["cdr_mask"].to(device),
                    epitope_mask=b["peptide_mask"].to(device),
                )
                preds.append(y.detach().float().cpu().numpy().ravel())
                done += len(preds[-1])
        out[i] = np.concatenate(preds)[:len(peptides)]
        if (i + 1) % log_every == 0 or i == len(cdr3s) - 1:
            rate = done / max(time.time() - t0, 1e-9)
            print(f"  {i+1}/{len(cdr3s)} betas, {done:,} pairs, "
                  f"{rate:,.0f} pairs/s", flush=True)
    return out, done / max(time.time() - t0, 1e-9)


@app.function(image=image, volumes={"/data": vol}, gpu="A10G", timeout=1800)
def smoke(n_peptides: int = 512, batch_size: int = 8):
    """Measure throughput before committing to a real run.

    One epoch of everything: loads the checkpoint, scores the target against a
    few hundred peptides, and reports pairs/s so the full run can be sized from
    a measured number rather than a guess.
    """
    import time

    import pandas as pd
    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"
    t0 = time.time()
    model, args = _load_model(device)
    load_s = time.time() - t0
    print(f"model loaded in {load_s:.0f}s on {device}", flush=True)

    pan = pd.read_csv("/opt/prp/data/netmhc_WBSB_9mers_clean.csv")
    peps = pan["Epitope"].astype(str).tolist()[:n_peptides]
    logits, rate = _score(model, args, [TARGET], peps, batch_size, device,
                          log_every=1)
    full = 293_467
    return {
        "device": device,
        "load_seconds": round(load_s),
        "batch_size": batch_size,
        "pairs_per_s": round(rate, 1),
        "logit_mean": float(logits.mean()),
        "logit_std": float(logits.std()),
        "hours_for_168_betas_full_panel": round(168 * full / rate / 3600, 1),
        "minutes_for_168_betas_6k_peptides": round(168 * 6000 / rate / 60, 1),
    }


# --------------------------------------------------------------------------
# batch scoring, resumable


@app.function(image=image, volumes={"/data": vol}, gpu="A10G", timeout=21600)
def score_betas(betas: list, peptides: list, tag: str, batch_size: int = 32):
    """Score every beta against one shared peptide set, one .npy per beta.

    Resumable by construction: a beta whose file already exists is skipped and
    the volume is committed as it goes, so a timeout or preemption costs only
    the betas still outstanding. That is the rule this repo learned by throwing
    away twelve finished epochs to a failure in an optional diagnostic.
    """
    import numpy as np
    import torch

    outdir = f"/data/prp/{tag}"
    os.makedirs(outdir, exist_ok=True)
    todo = [b for b in betas if not os.path.exists(f"{outdir}/{b}.npy")]
    print(f"{len(betas)} betas, {len(todo)} outstanding, "
          f"{len(peptides):,} peptides", flush=True)
    if not todo:
        return {"tag": tag, "scored": 0, "skipped": len(betas)}

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model, args = _load_model(device)
    np.save(f"{outdir}/_peptides.npy", np.array(peptides))
    vol.commit()

    rates = []
    for n, beta in enumerate(todo, 1):
        logits, rate = _score(model, args, [beta], peptides, batch_size,
                              device, log_every=10 ** 9)
        np.save(f"{outdir}/{beta}.npy", logits[0])
        rates.append(rate)
        if n % 10 == 0 or n == len(todo):
            vol.commit()
            eta = (len(todo) - n) * len(peptides) / max(float(np.mean(rates)), 1e-9)
            print(f"  {n}/{len(todo)}  {np.mean(rates):,.0f} pairs/s  "
                  f"eta {eta / 60:.0f} min", flush=True)
    vol.commit()
    return {"tag": tag, "scored": len(todo),
            "pairs_per_s": round(float(np.mean(rates)), 1)}


def _peptide_set(preds: str, repo: str, n_target: int, n_background: int,
                 seed: int = 0):
    """Target-recognised peptides plus a matched random background.

    The target set is TCR 4.1's top predictions. 4.1 is the patient-derived TCR
    carrying our exact clonotype, and its profile comes from a single-TCR model
    this run never uses, so the peptide set is defined independently of the
    joint model being tested.
    """
    import csv as _csv

    import numpy as np

    with open(os.path.join(repo, "data", "netmhc_WBSB_9mers_clean.csv")) as fh:
        panel = np.array([r["Epitope"] for r in _csv.DictReader(fh)])
    prof = np.load(os.path.join(preds, "4.1__y_pred_proteome.npy")).ravel()
    order = np.argsort(-prof)
    pos = order[:n_target]
    rng = np.random.default_rng(seed)
    bg = rng.choice(order[n_target * 10:], size=n_background, replace=False)
    peps = np.concatenate([panel[pos], panel[bg]])
    label = np.concatenate([np.ones(len(pos), int), np.zeros(len(bg), int)])
    keep = np.array(["X" not in p and len(p) == 9 for p in peps])
    return peps[keep].tolist(), label[keep]


@app.local_entrypoint()
def control(preds: str, repo: str, n_target: int = 3000,
            n_background: int = 3000, batch_size: int = 32):
    """Positive control: can the joint model rank TCRs we already have answers for?

    Scores the 15 panel betas the joint model never trained on. Each one's
    target-set AUC is then compared against its motif JS to the target, computed
    independently from the cached single-TCR profiles. If the joint model cannot
    reproduce that ordering, its ranking of our own candidates means nothing.
    """
    import csv as _csv

    import numpy as np

    peps, label = _peptide_set(preds, repo, n_target, n_background)
    with open(os.path.join(repo, "data", "tcr_cdr3b.csv")) as fh:
        rows = list(_csv.DictReader(fh))
    trained = {"19.2", "19.2C1", "19.2C2", "19.2C3", "19.2C4", "19.2C5"}
    betas = sorted({r["CDR3_b"] for r in rows if r["tcr_id"] not in trained}
                   | {TARGET})
    print(f"control: {len(betas)} betas unseen by the joint model, "
          f"{len(peps):,} peptides ({label.sum():,} target-recognised)")
    np.save("prp_control_labels.npy", label)
    print(score_betas.remote(betas, peps, "control", batch_size))


@app.local_entrypoint()
def candidates(preds: str, repo: str, panel: str, n_target: int = 3000,
               n_background: int = 3000, batch_size: int = 32):
    """Score our own AS candidate clonotypes against the target's peptide set.

    Descriptive only. The panel was selected using AS/healthy labels, so no
    p-value about AS association can come out of this ranking.
    """
    import csv as _csv

    import numpy as np

    peps, label = _peptide_set(preds, repo, n_target, n_background)
    with open(panel) as fh:
        betas = sorted({r["cdr3"] for r in _csv.DictReader(fh, delimiter="\t")}
                       | {TARGET})
    print(f"candidates: {len(betas)} betas, {len(peps):,} peptides "
          f"({label.sum():,} target-recognised)")
    np.save("prp_candidate_labels.npy", label)
    print(score_betas.remote(betas, peps, "candidates", batch_size))
