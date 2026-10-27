"""Does functional similarity between TCRs track sequence similarity? No.

Why this exists. This project's two attempts to widen `CASSVGLYSTDTQYF` into a
motif family both failed, and both failed the same way: the metric was spelling.
Edit distance makes 74% of healthy subjects carriers at radius 1. ESM-2 cosine
returns anagrams, because mean-pooling discards residue order. Neither has any
notion of what a receptor actually recognises, so neither could say which
neighbours are relatives and which are common public sequence.

Wang, Yeh & Garcia (*Nat Biotechnol*, 2026) supply the missing ground truth.
They ran yeast-display selections with a ~1e9-variant 9-mer library against 21
HLA-B*27:05-restricted TCRs from ankylosing spondylitis and acute anterior
uveitis patients, then fine-tuned a TCR-by-peptide interaction network on each
TCR's screen. Their panel is, by luck, a characterised sample of exactly our
motif family: it contains our target twice (`4.1`, `19.2C5`) and four of its
one-edit neighbours.

A *peptide recognition profile* (PRP) is a TCR's vector of predicted binding
scores across a peptide panel -- here 293,467 human-proteome 9-mers filtered to
NetMHC-predicted B27 binders. Two TCRs are functionally similar if they
recognise the same peptides. That is a derived quantity: you cannot read it off
sequence, which is the whole reason our sequence-based attempts could not get at
it.

So this asks, on their numbers and none of ours: within one motif family, does
functional similarity track sequence similarity? If it does not, the failure of
our widening attempts was a property of the biology rather than a bad choice of
metric.

Three things make this nearly free. The profiles are cached model outputs, so
there is no model to run and no GPU. They are ~1.2 MB inside 3.2 GB archives
(each archive carries a fine-tuned ESM2-650M), so `hfzip.py` range-fetches only
the member wanted -- 32 MB instead of 70.6 GB. And the 6 pairs of TCRs that
share an identical CDR3-beta fall out as the zero-edit row of the same matrix,
which is the test of whether beta alone determines specificity.

Pre-registered before looking at any number:

- These are **model predictions**, not measurements. Each profile is one model's
  interpolation of its own TCR's screen. `--activation` therefore cross-checks
  profile similarity against the 16-peptide wet-lab T-cell activation table,
  which is the only pure functional readout in the deposit. If the two disagree,
  the profiles are measuring the model and not the biology.
- All 21 TCRs sit in one motif family (edit distances 0-7), so this tests
  sequence-versus-function *within* a family. That is our question, but it is
  narrower than the paper's general claim.
- Each model was trained on its own data volume and seed, so a zero-edit pair
  can differ for reasons that have nothing to do with the alpha chain. Without
  replicate models per TCR that noise floor cannot be measured, so the zero-edit
  row is read *relative to* the one- and two-edit rows rather than absolutely.

    python scripts/prp_divergence.py --preds <dir> --repo <clone>
"""
from __future__ import annotations

import argparse
import csv
import itertools
import os
import sys

import numpy as np
from scipy.stats import spearmanr

TARGET = "CASSVGLYSTDTQYF"
TOP_FRAC = 0.01          # a TCR's "recognised set": its top 1% of peptides


def levenshtein(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def jaccard(top_a: np.ndarray, top_b: np.ndarray) -> float:
    """Overlap of two recognised sets, as |A and B| / |A or B|."""
    inter = np.intersect1d(top_a, top_b, assume_unique=True).size
    return inter / (top_a.size + top_b.size - inter)


AA = "ACDEFGHIKLMNPQRSTVWY"
# The paper excludes P2, P8 and P9 (1-indexed): P2 is the canonical B27 anchor
# (fixed Arg), and P8/P9 were held fixed during library construction, so those
# positions carry almost no amino-acid variability among binders. Including them
# would add three near-constant columns that dilute every pairwise distance.
EXCLUDED_POS = (2, 8, 9)


def motif(peptides: list[str], length: int = 9) -> np.ndarray:
    """Position-weight matrix over the included positions, shape (|I|, 20)."""
    keep = [i for i in range(1, length + 1) if i not in EXCLUDED_POS]
    idx = {a: j for j, a in enumerate(AA)}
    m = np.zeros((len(keep), len(AA)))
    for p in peptides:
        if len(p) != length:
            continue
        for row, pos in enumerate(keep):
            a = p[pos - 1]
            if a in idx:
                m[row, idx[a]] += 1
    return m / np.maximum(m.sum(axis=1, keepdims=True), 1e-12)


def js_distance(p: np.ndarray, q: np.ndarray) -> float:
    """Mean per-position Jensen-Shannon distance, as the paper defines it.

        d_JS(p,q) = (1/|I|) * sum_i sqrt( (KL(p_i||m_i) + KL(q_i||m_i)) / 2 )

    with m_i the pointwise mean of the two position distributions. 0 means
    identical binding motifs; larger means more functionally divergent.
    """
    m = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_p = np.where(p > 0, p * np.log(p / m), 0.0).sum(axis=1)
        kl_q = np.where(q > 0, q * np.log(q / m), 0.0).sum(axis=1)
    return float(np.mean(np.sqrt(np.maximum(0.5 * (kl_p + kl_q), 0.0))))


def load(preds: str, repo: str):
    """Returns (tcr ids, CDR3b per id, profile matrix)."""
    with open(os.path.join(repo, "data", "tcr_cdr3b.csv")) as fh:
        cdr3 = {r["tcr_id"]: r["CDR3_b"] for r in csv.DictReader(fh)}

    ids, rows = [], []
    for tcr in sorted(cdr3):
        path = os.path.join(preds, f"{tcr}__y_pred_proteome.npy")
        if not os.path.exists(path):
            print(f"  missing profile for {tcr}, skipping", file=sys.stderr)
            continue
        ids.append(tcr)
        rows.append(np.load(path).ravel())
    return ids, cdr3, np.vstack(rows)


def activation_similarity(repo: str, ids: list[str]):
    """TCR-by-TCR Jaccard over the 16 peptides measured in cells.

    This is the only readout here that is not a model output, so it is the one
    check on whether profile similarity reflects function or fitting.
    """
    path = os.path.join(repo, "data", "activation_binary.csv")
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    sets = {t: {r["Peptide"] for r in rows if r.get(t) == "1"} for t in ids}
    return sets


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds", required=True, help="dir of *__y_pred_proteome.npy")
    ap.add_argument("--repo", required=True, help="PRP-TCR-Specificity clone")
    ap.add_argument("--top-frac", type=float, default=TOP_FRAC)
    args = ap.parse_args()

    ids, cdr3, prof = load(args.preds, args.repo)
    n_pep = prof.shape[1]
    k = int(round(args.top_frac * n_pep))
    print(f"{len(ids)} TCRs x {n_pep:,} peptides; "
          f"recognised set = top {k:,} ({args.top_frac:.0%})\n")

    # Rank once per profile: Spearman between profiles is then Pearson on ranks.
    ranks = np.vstack([np.argsort(np.argsort(p)) for p in prof]).astype(np.float64)
    ranks -= ranks.mean(axis=1, keepdims=True)
    ranks /= np.linalg.norm(ranks, axis=1, keepdims=True)
    rho = ranks @ ranks.T

    tops = [np.sort(np.argpartition(-p, k)[:k]) for p in prof]
    act = activation_similarity(args.repo, ids)

    # The paper's metric operates on binding *motifs*, not on the raw score
    # vector. That matters: a Spearman over all 293,467 predictions is mostly
    # about how each model ranks the ~290k non-binders, which is per-model
    # idiosyncrasy rather than biology. A motif is a 6x20 summary of what the
    # TCR actually binds, and it is what reproduces the published result.
    with open(os.path.join(args.repo, "data",
                           "netmhc_WBSB_9mers_clean.csv")) as fh:
        panel = np.array([r["Epitope"] for r in csv.DictReader(fh)])
    mots = [motif(list(panel[t])) for t in tops]

    recs = []
    for i, j in itertools.combinations(range(len(ids)), 2):
        a, b = ids[i], ids[j]
        sa, sb = act[a], act[b]
        union = len(sa | sb)
        recs.append({
            "a": a, "b": b,
            "edit": levenshtein(cdr3[a], cdr3[b]),
            "js": js_distance(mots[i], mots[j]),
            "rho": float(rho[i, j]),
            "jac": jaccard(tops[i], tops[j]),
            "act_jac": (len(sa & sb) / union) if union else float("nan"),
        })

    edits = np.array([r["edit"] for r in recs])
    rhos = np.array([r["rho"] for r in recs])
    jacs = np.array([r["jac"] for r in recs])
    jss = np.array([r["js"] for r in recs])

    # Positive control, first and non-negotiable. 19.2C1-C5 are 19.2 engineered
    # with 1-3 CDR3b substitutions on its original alpha chain, and the paper
    # (Fig. 5a) reports they form a tight functional neighbourhood. If this
    # reimplementation cannot recover that, its numbers mean nothing.
    fam = {"19.2", "19.2C1", "19.2C2", "19.2C3", "19.2C4", "19.2C5"}
    within = np.array([r["js"] for r in recs if r["a"] in fam and r["b"] in fam])
    print("=== POSITIVE CONTROL: the engineered 19.2 neighbourhood ===")
    print(f"  median JS within the 19.2 family : {np.median(within):.4f} "
          f"({len(within)} pairs)")
    print(f"  median JS over all pairs         : {np.median(jss):.4f}")
    print("  recovered" if np.median(within) < np.median(jss) else
          "  NOT RECOVERED -- the metric is broken; stop here")

    print("\n=== does sequence distance predict functional similarity? ===")
    for name, vals in (("motif JS (the paper's)", jss),
                       ("profile Spearman", rhos),
                       ("top-1% Jaccard", jacs)):
        s = spearmanr(edits, vals)
        print(f"  {name:24s} vs edit distance: rho = {s.correlation:+.3f}, "
              f"p = {s.pvalue:.3g}  (n = {len(recs)} pairs)")

    print("\n=== functional similarity by edit distance ===")
    print(f"  {'edits':>5} {'pairs':>6} {'median JS':>10} "
          f"{'median rho':>11} {'median Jaccard':>15}")
    for e in sorted(set(edits.tolist())):
        m = edits == e
        print(f"  {e:>5} {m.sum():>6} {np.median(jss[m]):>10.4f} "
              f"{np.median(rhos[m]):>11.3f} {np.median(jacs[m]):>15.3f}")

    print("\n=== zero-edit pairs: identical CDR3-beta, different TCR ===")
    print("  (read against the one- and two-edit rows above, not absolutely)")
    for r in sorted((r for r in recs if r["edit"] == 0), key=lambda r: r["js"]):
        print(f"  {r['a']:>7} / {r['b']:<8} {cdr3[r['a']]:<16} "
              f"JS {r['js']:.4f}  rho {r['rho']:+.3f}")

    print(f"\n=== our target {TARGET}: everything ranked by profile similarity ===")
    tgt = [t for t in ids if cdr3[t] == TARGET]
    print(f"  carried by {tgt}")
    for t in tgt:
        print(f"\n  against {t}:")
        rows = [r for r in recs if t in (r["a"], r["b"])]
        for r in sorted(rows, key=lambda r: r["js"]):
            other = r["b"] if r["a"] == t else r["a"]
            print(f"    JS {r['js']:.4f}  rho {r['rho']:+.3f}  "
                  f"{r['edit']} edits  {other:>7}  {cdr3[other]}")

    print("\n=== control: do model profiles agree with wet-lab activation? ===")
    ok = np.array([not np.isnan(r["act_jac"]) for r in recs])
    if ok.sum() > 2:
        aj = np.array([r["act_jac"] for r in recs])[ok]
        s = spearmanr(rhos[ok], aj)
        print(f"  profile Spearman vs measured-activation Jaccard: "
              f"rho = {s.correlation:+.3f}, p = {s.pvalue:.3g} "
              f"(n = {ok.sum()} pairs, 16 peptides)")
        s2 = spearmanr(edits[ok], aj)
        print(f"  edit distance    vs measured-activation Jaccard: "
              f"rho = {s2.correlation:+.3f}, p = {s2.pvalue:.3g}")
        print("  If the first beats the second, PRP similarity carries real")
        print("  functional information that spelling does not.")

    out = os.path.join(args.preds, "prp_pairs.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(recs)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
