"""Does training a model across studies beat a hand-specified gene score?

The honest test is leave-one-STUDY-out: train on four cohorts, predict the
fifth. Splitting samples at random instead would let the model memorise each
study's platform and batch, which is why cross-validated AUCs in this
literature are often near 1.0 and mean nothing.

Three contestants, all scored on the same held-out studies:

  1. pre-specified NK effector score -- no training at all, 11 genes chosen
     from biology before seeing any performance number
  2. elastic-net logistic regression on the full shared transcriptome
  3. random forest on the same

Feature selection and all fitting happen strictly inside the training folds.

Usage:  python scripts/ml_loso.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bulk_meta import PREP, PRIMARY  # noqa: E402
from deconv_bulk import EFFECTOR  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "bulk" / "meta"
SEED = 0


def load_all():
    """Per-study z-scored expression on the genes shared by every study.

    z-scoring is per study and uses only expression, never labels -- it is the
    minimum needed to put a microarray intensity and a log-CPM on one scale.
    """
    expr, meta = {}, {}
    for gse in PRIMARY:
        e = pd.read_csv(PREP / f"{gse}_expr.tsv.gz", sep="\t", index_col=0)
        m = pd.read_csv(PREP / f"{gse}_samples.tsv", sep="\t", index_col=0)
        e = e[[c for c in e.columns if c in m.index]]
        expr[gse], meta[gse] = e, m.loc[e.columns]

    shared = set.intersection(*(set(e.index) for e in expr.values()))
    shared = sorted(shared)
    print(f"genes shared by all {len(PRIMARY)} studies: {len(shared)}")

    X, y, grp = [], [], []
    for gse in PRIMARY:
        e = expr[gse].loc[shared]
        z = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, np.nan), axis=0)
        z = z.fillna(0.0)
        X.append(z.T)
        y.append((meta[gse]["group"] == "AS").astype(int).to_numpy())
        grp.append(np.full(z.shape[1], gse))
    return (pd.concat(X), np.concatenate(y), np.concatenate(grp), shared)


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    X, y, grp, genes = load_all()
    print(f"samples: {X.shape[0]}  (AS {int(y.sum())}, HC {int((1 - y).sum())})\n")

    eff = [g for g in EFFECTOR if g in X.columns]
    print(f"pre-specified NK effector score uses {len(eff)} genes: {eff}\n")

    models = {
        "elastic net": Pipeline([
            ("sc", StandardScaler()),
            ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                       l1_ratio=0.5, C=0.05, max_iter=5000,
                                       random_state=SEED)),
        ]),
        "random forest": Pipeline([
            ("sc", StandardScaler()),
            ("clf", RandomForestClassifier(n_estimators=500, min_samples_leaf=3,
                                           random_state=SEED, n_jobs=-1)),
        ]),
    }

    rows = []
    for held in PRIMARY:
        tr, te = grp != held, grp == held
        ytr, yte = y[tr], y[te]
        if len(np.unique(yte)) < 2:
            continue
        base = max(yte.mean(), 1 - yte.mean())

        # Feature selection INSIDE the training fold only.
        f, _ = f_classif(X[tr], ytr)
        top = np.argsort(np.nan_to_num(f))[::-1][:2000]

        res = {"held_out": held, "n_test": int(te.sum()),
               "majority_baseline": base}

        s = X.loc[te.tolist() if False else te, eff].mean(axis=1)
        res["NK score (untrained)"] = roc_auc_score(yte, -s)

        for name, mdl in models.items():
            mdl.fit(X.iloc[tr, top], ytr)
            p = mdl.predict_proba(X.iloc[te, top])[:, 1]
            res[name] = roc_auc_score(yte, p)
        rows.append(res)
        print(f"  held out {held:10s} n={int(te.sum()):3d}  "
              f"NK score {res['NK score (untrained)']:.3f}   "
              f"elastic net {res['elastic net']:.3f}   "
              f"random forest {res['random forest']:.3f}   "
              f"(majority baseline accuracy {base:.2f})")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "ml_loso.tsv", sep="\t", index=False)

    print("\n" + "=" * 90)
    print("MEAN AUC ACROSS HELD-OUT STUDIES")
    print("=" * 90)
    for c in ("NK score (untrained)", "elastic net", "random forest"):
        v = df[c]
        print(f"  {c:24s} {v.mean():.3f}   (range {v.min():.3f}-{v.max():.3f})")

    print("\n  Same-study cross-validation, for contrast -- this is the number")
    print("  that gets published when the split ignores study membership:")
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    f, _ = f_classif(X, y)
    top = np.argsort(np.nan_to_num(f))[::-1][:2000]
    for name, mdl in models.items():
        s = cross_val_score(mdl, X.iloc[:, top], y, cv=cv, scoring="roc_auc")
        print(f"    {name:16s} AUC {s.mean():.3f}  <-- leaks: feature selection on all"
              f" data, and each study appears in train and test")


if __name__ == "__main__":
    main()
