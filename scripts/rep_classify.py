"""Positive control: predict HLA-A*02 carriage from a TCR-beta repertoire.

Why this target. HLA molecules present peptides, and a thymocyte only survives
selection if its receptor binds self-HLA with the right affinity. So a person's
HLA type is a direct constraint on which receptors can exist in them -- the
label is receptor-encoded by construction, not by hypothesis. Emerson et al.
(Nat Genet 2017) report recovering HLA type from the repertoire in this same
cohort. The labels ship with the data: `subject.mhc` in the ADC metadata.

What a positive control is for. Every result in this project so far has been a
negative measured against a shuffled baseline, which cannot tell "no signal"
from "broken pipeline". If this script cannot find HLA-A*02 -- a signal that is
definitely there -- then nothing it says about a disease label means anything.
This runs first, and the disease work does not start until it passes.

The leakage trap. The feature family that carries the signal is "which specific
public CDR3s does this subject have", and choosing which CDR3s to look at
requires the labels. Selecting them once over all subjects and then
cross-validating would report a near-perfect AUC on pure noise. So selection
happens inside each training fold only, and the whole procedure -- selection
included -- is re-run under permuted labels to show what it produces when
there is nothing to find.

Reading the numbers. An AUC needs a floor and a ceiling. The floor is the
permutation null (what this exact pipeline returns on scrambled labels; if
selection leaks, this sits well above 0.5). The ceiling is the published
result for this cohort. The honest summary is where the real number sits
between them.

Usage:
    python scripts/rep_classify.py
    python scripts/rep_classify.py --allele HLA-B*07 --permutations 20
    python scripts/rep_classify.py --blocks summary       # no public clones
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import hypergeom
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
# Two kinds of feature get confused with each other, and separating them
# changes what a baseline means.
#
# TECHNICAL is how much material was sequenced and how much came back. It has
# no biological content: if it predicts a label, the label is tracking sample
# handling and the result is an artifact.
#
# CLONAL_STRUCTURE is the shape of the clone-size distribution. These are
# depth-*influenced* -- a deeper sample sees more rare clones -- but they are
# not depth, they are biology. CMV is the clean illustration: it drives large
# oligoclonal CD8 expansions ("memory inflation"), so CMV+ repertoires really
# are more clonal. Grouping these with TECHNICAL once made a genuine
# biological effect look like a sequencing artifact.
TECHNICAL = ("n_clones", "n_templates", "n_unique_cdr3",
             "log_n_clones", "log_n_templates")
CLONAL_STRUCTURE = ("shannon", "shannon_norm", "clonality", "gini",
                    "top1_share", "top10_share", "top100_share",
                    "n_expanded_gt10", "singleton_frac")
DEPTH_SENSITIVE = TECHNICAL  # what --drop-depth removes
SEED = 0


# --------------------------------------------------------------------------
# labels


def build_labels(meta: pd.DataFrame, allele: str,
                 require_two: bool = False) -> pd.DataFrame:
    """Carriage of one HLA allele, from documented typing only.

    `subject.mhc` mixes typed alleles with Emerson's own model-inferred ones
    ("Inferred HLA-A*01"). Using an inferred allele as a label would be
    scoring a model against another model's output, so only documented entries
    count -- `adc_fetch.py` keeps the two apart.

    A subject counts as negative only if they have documented typing at that
    *gene*. Someone with no HLA-A call at all is not an A*02 negative, they are
    unknown, and is dropped. Calling them negative would put real A*02
    carriers in the control group and depress the result for no reason.
    """
    gene = allele.split("*")[0]              # 'HLA-A*02' -> 'HLA-A'
    doc = meta["hla_documented"].fillna("")
    alleles = doc.str.split(";")
    n_at_gene = alleles.map(lambda xs: sum(a.startswith(gene + "*") for a in xs))
    carries = alleles.map(lambda xs: allele in xs)
    out = meta.assign(n_typed=n_at_gene, y=carries.astype(int))
    # In this cohort 82 of 626 subjects have only one documented HLA-A allele.
    # That is either homozygosity or incomplete typing, and the two are not
    # distinguishable from the metadata. If it is incomplete typing, a real
    # A*02 carrier can be sitting in the control group. `require_two` restricts
    # to the 544 fully typed subjects as a sensitivity check; the effect should
    # be a slightly *higher* AUC on cleaner labels, not a lower one.
    keep = out["n_typed"] >= (2 if require_two else 1)
    return out[keep].copy()


def build_cohort_labels(meta: pd.DataFrame) -> pd.DataFrame:
    """Which sub-cohort is this subject from? A label with no biology in it.

    Emerson's 786 subjects are split into Cohort 01 (666) and Cohort 02 (120).
    Both are healthy donors, same lab, same Adaptive protocol, same tissue,
    same assay. There is no biological difference between the groups, so a
    classifier that separates them is reading batch -- collection site, reagent
    lot, sequencing run, processing date -- and nothing else.

    This is a *negative* control, and the inverse of the HLA run: there the
    answer had to be findable, here it has to be absent. It matters because the
    project's intended disease comparisons draw cases from one study and
    controls from another. If cohort membership is predictable at 0.9 within a
    single lab, then any cross-study case/control AUC is measuring the study.
    Here the confound can be isolated completely, because there is no biology
    for it to hide behind.
    """
    y = (meta["study_group"].astype(str).str.strip() == "Cohort 02").astype(int)
    return meta.assign(y=y, n_typed=2)


EXTERNAL = Path("data/external/airr_hip/metadata.txt")


def load_external(root: Path) -> pd.DataFrame:
    """Per-subject CMV serostatus, age, sex, from a mirror of this cohort.

    The ADC copy of Emerson carries HLA typing but not CMV serostatus, so the
    labels come from a redistribution of the same cohort that does. An external
    label table is exactly the kind of join that fails silently, so it is
    verified rather than trusted: this file has its own `hla` column, and
    HLA-A*02 derived from it agrees with HLA-A*02 derived from the ADC metadata
    for all 626 subjects typed in both, with zero disagreements. Same cohort,
    same subjects, correct join.
    """
    ext = pd.read_csv(root / EXTERNAL, sep="\t")
    ext = ext.rename(columns={"sample_id": "subject_id"})
    ext["age"] = pd.to_numeric(ext["age"], errors="coerce")
    ext["is_male"] = (ext["sex"] == "male").astype(float)
    return ext[["subject_id", "age", "is_male", "sex", "cmv", "hla"]]


def build_cmv_labels(meta: pd.DataFrame, ext: pd.DataFrame) -> pd.DataFrame:
    """CMV seropositive vs seronegative. 25 subjects are untested and dropped.

    Unlike HLA, this is an *acquired* label: a record of whether this person
    has met a particular virus. That makes it the right analogue of a disease
    label, where HLA -- a germline constraint on which receptors can exist at
    all -- is not.

    It is also confounded with age. Seropositive subjects here are 4.7 years
    older on average (p = 4.5e-5), and age reshapes the repertoire through
    thymic involution, so some of any CMV signal is really age. `main` reports
    an age-and-sex-only baseline for this reason.
    """
    m = meta.merge(ext, on="subject_id", how="left")
    m = m[m["cmv"].isin(["+", "-"])].copy()
    m["y"] = (m["cmv"] == "+").astype(int)
    m["n_typed"] = 2
    return m


# --------------------------------------------------------------------------
# public-clone burden, selected inside the fold


def select_clones(X: sp.csr_matrix, y: np.ndarray, p_cut: float) -> np.ndarray:
    """Clones over-represented in label-positive training subjects.

    One-sided Fisher exact test per clone, done as a hypergeometric survival
    function so all candidates are tested at once. For a clone carried by
    `carriers` of `N` training subjects, of whom `K` are positive and `k` of the
    carriers are positive, the p-value is P(>= k | no association).
    """
    N = X.shape[0]
    K = int(y.sum())
    carriers = np.asarray(X.sum(axis=0)).ravel()
    k = np.asarray(X[y == 1].sum(axis=0)).ravel()
    # sf(k-1) == P(X >= k)
    p = hypergeom.sf(k - 1, N, K, carriers)
    return np.flatnonzero(p < p_cut)


def burden_features(X: sp.csr_matrix, idx: np.ndarray, depth: np.ndarray) -> np.ndarray:
    """How many of the selected clones a subject carries.

    Three columns, because the raw count confounds two things: a subject with a
    deeper sample finds more of everything. `burden` is the raw hit count,
    `burden_frac` divides by how many clones were selected (a per-subject rate),
    and `burden_per_depth` divides by the subject's own clone count, which is
    the comparison Emerson's classifier makes.
    """
    hits = np.asarray(X[:, idx].sum(axis=1)).ravel().astype(float)
    n_sel = max(len(idx), 1)
    return np.column_stack([hits, hits / n_sel, hits / depth])


# --------------------------------------------------------------------------
# one evaluation


def fit_fold(Xtr_pub, ytr, Xte_pub, Str, Ste, depth_tr, depth_te,
             p_cut: float, use_public: bool, use_summary: bool):
    """Train on one fold and return test scores. Nothing here sees test labels."""
    tr_blocks, te_blocks = [], []
    n_sel = 0
    if use_public:
        idx = select_clones(Xtr_pub, ytr, p_cut)
        n_sel = len(idx)
        if n_sel == 0:
            # No clone passed. Honest outcome: contribute a constant, which is
            # chance, rather than silently falling back to the other block.
            tr_blocks.append(np.zeros((Xtr_pub.shape[0], 1)))
            te_blocks.append(np.zeros((Xte_pub.shape[0], 1)))
        else:
            tr_blocks.append(burden_features(Xtr_pub, idx, depth_tr))
            te_blocks.append(burden_features(Xte_pub, idx, depth_te))
    if use_summary:
        tr_blocks.append(Str)
        te_blocks.append(Ste)

    Xtr = np.hstack(tr_blocks)
    Xte = np.hstack(te_blocks)
    model = make_pipeline(
        StandardScaler(),
        # L2 logistic regression -- L2 is the default, and passing
        # penalty="l2" explicitly is deprecated as of scikit-learn 1.8. The
        # same API drift (`multi_class=`) once cost this project a finished
        # 12-epoch training run, so arguments here stay at their defaults
        # unless there is a reason.
        #
        # Nothing larger is warranted: the cell-level work established that
        # plain features tie a 1.5M-parameter transformer on receptor data.
        LogisticRegression(C=1.0, max_iter=5000, solver="lbfgs"),
    )
    model.fit(Xtr, ytr)
    return model.predict_proba(Xte)[:, 1], n_sel


def cross_validate(Xpub, S, y, depth, folds: int, p_cut: float,
                   use_public: bool, use_summary: bool, seed: int = SEED):
    """Out-of-fold scores for every subject, then one AUC over all of them."""
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    sel_counts = []
    for tr, te in skf.split(np.zeros(len(y)), y):
        scores, n_sel = fit_fold(
            Xpub[tr] if use_public else None, y[tr],
            Xpub[te] if use_public else None,
            S[tr] if use_summary else None, S[te] if use_summary else None,
            depth[tr], depth[te], p_cut, use_public, use_summary)
        oof[te] = scores
        sel_counts.append(n_sel)
    return roc_auc_score(y, oof), oof, sel_counts


def boot_ci(y: np.ndarray, s: np.ndarray, n: int = 2000, seed: int = SEED):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) < 2:
            continue
        out.append(roc_auc_score(y[i], s[i]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


# --------------------------------------------------------------------------


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="emerson-2017-natgen")
    ap.add_argument("--target", default="allele", choices=["allele", "cohort", "cmv"],
                    help="'cohort' predicts which sub-cohort a subject is from, "
                         "a label with no biological content (see build_cohort_labels)")
    ap.add_argument("--heldout-cohort", action="store_true",
                    help="train on Cohort 01 and test on Cohort 02 instead of CV "
                         "-- the split Emerson's own paper used")
    ap.add_argument("--allele", default="HLA-A*02")
    ap.add_argument("--public", default=None,
                    help="incidence npz stem (default: the only one present)")
    ap.add_argument("--folds", type=int, default=10)
    ap.add_argument("--p-cut", type=float, default=1e-4,
                    help="Fisher p threshold for selecting a clone (Emerson used 1e-4)")
    ap.add_argument("--blocks", default="both", choices=["public", "summary", "both"])
    ap.add_argument("--permutations", type=int, default=20)
    ap.add_argument("--drop-depth", action="store_true",
                    help="exclude depth-sensitive summary features")
    ap.add_argument("--require-two-alleles", action="store_true",
                    help="only subjects with both alleles typed at this gene")
    args = ap.parse_args()

    base = ROOT / "data" / "airr" / args.study
    meta = pd.read_csv(base / "metadata.tsv", sep="\t")
    feats = pd.read_parquet(base / "features_summary.parquet")

    if args.target == "cohort":
        label_name = "cohort membership (negative control: no biology)"
        lab = build_cohort_labels(meta)
    elif args.target == "cmv":
        label_name = "CMV serostatus"
        lab = build_cmv_labels(meta, load_external(ROOT))
    else:
        label_name = args.allele
        lab = build_labels(meta, args.allele, args.require_two_alleles)
    lab = lab[lab["subject_id"].isin(feats["subject_id"])]
    print(f"label: {label_name}")
    print(f"  {len(meta)} repertoires, {len(lab)} usable with features on disk")
    print(f"  positives {int(lab['y'].sum())}, negatives {int((1 - lab['y']).sum())} "
          f"({lab['y'].mean():.1%} prevalence)")
    by = lab.groupby("study_group")["y"].agg(["size", "sum"])
    for grp, r in by.iterrows():
        print(f"    {grp}: {int(r['size'])} subjects, {int(r['sum'])} positive")
    if lab["y"].sum() < 20 or (1 - lab["y"]).sum() < 20:
        print("  too few in one class to be worth running")
        return

    # Align every matrix to one subject order.
    stems = sorted(p.stem for p in (base / "reps").glob("*.parquet"))
    pub_files = sorted(base.glob("public_*.npz"))
    use_public = args.blocks in ("public", "both")
    Xpub = None
    if use_public:
        if not pub_files:
            print("  no public-clone matrix found -- run rep_features.py --pass public")
            return
        stem = args.public or pub_files[-1].stem
        M = sp.load_npz(base / f"{stem}.npz").tocsr()
        pub_subj = pd.read_csv(base / f"{stem}.subjects.tsv", sep="\t")["subject_id"].tolist()
        print(f"  public clones: {stem}, {M.shape[1]:,} candidates")
        pos = pd.Index(pub_subj).get_indexer(lab["subject_id"])
        assert (pos >= 0).all(), "subject missing from incidence matrix"
        Xpub = M[pos]

    feats = feats.set_index("subject_id").loc[lab["subject_id"]]
    depth = feats["n_clones"].to_numpy(dtype=float)
    cols = [c for c in feats.columns
            if not (args.drop_depth and c in DEPTH_SENSITIVE)]
    S = feats[cols].to_numpy(dtype=float)
    y = lab["y"].to_numpy()
    use_summary = args.blocks in ("summary", "both")
    print(f"  summary features: {S.shape[1]}"
          f"{' (depth-sensitive dropped)' if args.drop_depth else ''}")

    results = {}

    # --- the transfer test ------------------------------------------------
    # Train on one cohort, test on the other. This is the measurement the
    # project most needs: cross-cohort transfer is where every previous result
    # died (CD4-vs-CD8 fell 0.658 -> 0.528 on a held-out study). The gap
    # between this and the CV number below is the cost of changing cohort.
    if args.heldout_cohort:
        grp = lab["study_group"].astype(str).str.strip().to_numpy()
        tr = np.flatnonzero(grp == "Cohort 01")
        te = np.flatnonzero(grp == "Cohort 02")
        if len(te) < 20 or len(np.unique(y[te])) < 2:
            print("\nheld-out cohort unusable: "
                  f"{len(te)} test subjects, {len(np.unique(y[te]))} classes")
        else:
            scores, n_sel = fit_fold(
                Xpub[tr] if use_public else None, y[tr],
                Xpub[te] if use_public else None,
                S[tr] if use_summary else None, S[te] if use_summary else None,
                depth[tr], depth[te], args.p_cut, use_public, use_summary)
            a = roc_auc_score(y[te], scores)
            lo, hi = boot_ci(y[te], scores)
            print(f"\ntrain Cohort 01 (n={len(tr)}, {int(y[tr].sum())} pos) -> "
                  f"test Cohort 02 (n={len(te)}, {int(y[te].sum())} pos)")
            print(f"  AUC {a:.3f}  [{lo:.3f}, {hi:.3f}]   ({n_sel:,} clones selected)")
            results["auc_heldout_cohort"] = float(a)
            results["ci_heldout_cohort"] = [lo, hi]

    # --- the real run -----------------------------------------------------
    t0 = time.time()
    auc, oof, sel = cross_validate(Xpub, S, y, depth, args.folds, args.p_cut,
                                   use_public, use_summary)
    lo, hi = boot_ci(y, oof)
    print(f"\n{args.folds}-fold CV, blocks={args.blocks}")
    print(f"  AUC {auc:.3f}  [{lo:.3f}, {hi:.3f}]   ({time.time() - t0:.0f}s)")
    if use_public:
        print(f"  clones selected per fold: {np.median(sel):,.0f} "
              f"(min {min(sel):,}, max {max(sel):,})")
    results["auc"] = auc
    results["ci"] = [lo, hi]
    results["clones_selected"] = [int(s) for s in sel]

    # --- the floor: same pipeline, scrambled labels -----------------------
    # This is the number that makes the result interpretable. Selection is
    # re-run inside every fold of every permutation, so if the selection step
    # leaks test information this sits well above 0.5 and the real AUC above
    # means nothing.
    if args.permutations:
        rng = np.random.default_rng(SEED)
        null = []
        t0 = time.time()
        for i in range(args.permutations):
            yp = rng.permutation(y)
            a, _, _ = cross_validate(Xpub, S, yp, depth, args.folds, args.p_cut,
                                     use_public, use_summary, seed=SEED + i)
            null.append(a)
            print(f"    permutation {i + 1}/{args.permutations}: {a:.3f}", flush=True)
        null = np.array(null)
        print(f"  permutation null: mean {null.mean():.3f}, "
              f"95th pct {np.percentile(null, 95):.3f}, max {null.max():.3f} "
              f"({time.time() - t0:.0f}s)")
        # Permutation p-value, with the +1 correction so it can never be 0.
        p = (1 + int((null >= auc).sum())) / (1 + len(null))
        print(f"  permutation p <= {p:.3f}")
        results["null_mean"] = float(null.mean())
        results["null_max"] = float(null.max())
        results["null_p"] = p

    # --- a depth-only baseline -------------------------------------------
    # Sequencing depth is the confounder that ruins repertoire analyses. HLA
    # carriage has no reason to track depth, so this should be at chance; if it
    # is not, something about the cohort is structured and must be chased down
    # before any disease label is trusted.
    for name, block, note in (
            ("technical-only", TECHNICAL, "should be ~0.5; artifact if not"),
            ("clonal-structure", CLONAL_STRUCTURE, "real biology, not an artifact")):
        cs = [c for c in feats.columns if c in block]
        a_b, _, _ = cross_validate(None, feats[cs].to_numpy(dtype=float), y, depth,
                                   args.folds, args.p_cut, False, True)
        print(f"  {name + ' baseline:':<26} AUC {a_b:.3f}  "
              f"({len(cs)} features; {note})")
        results[f"auc_{name.replace('-', '_')}"] = float(a_b)

    # --- a demographics baseline -----------------------------------------
    # Age and sex, no repertoire at all. For CMV this is not a nuisance check
    # but a competitor: seroprevalence rises with age, so a repertoire model
    # has to beat what you get from asking the patient their age. 87 subjects
    # have no recorded age and take the median, which if anything flatters this
    # baseline's coverage rather than its accuracy.
    if {"age", "is_male"}.issubset(lab.columns):
        dem = lab[["age", "is_male"]].copy()
        dem["age"] = dem["age"].fillna(dem["age"].median())
        a_dem, _, _ = cross_validate(None, dem.to_numpy(dtype=float), y, depth,
                                     args.folds, args.p_cut, False, True)
        print(f"  age+sex baseline:    AUC {a_dem:.3f}  (no repertoire data)")
        results["auc_demographics"] = float(a_dem)

    tag = {"cohort": "cohort", "cmv": "cmv"}.get(
        args.target, args.allele.replace("*", "").replace(":", ""))
    if args.heldout_cohort:
        tag += "_heldout"
    if args.drop_depth:
        tag += "_nodepth"
    if args.require_two_alleles:
        tag += "_twoalleles"
    out = base / f"classify_{tag}_{args.blocks}.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
