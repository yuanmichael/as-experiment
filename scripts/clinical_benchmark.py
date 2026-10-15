"""If you ran bulk RNA-seq on one person's blood, what would it actually tell you?

An effect size answers "do these groups differ on average". A clinician needs
"given THIS person's result, what is the chance they have the disease" -- which
depends on sensitivity, specificity AND how common the disease is in the people
being tested. The last term is the one that decides whether a test is useful,
and it is the one effect sizes never show.

Benchmarked against CRP, the blood test actually used in axSpA.

Usage:  python scripts/clinical_benchmark.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bulk_meta import PREP, PRIMARY  # noqa: E402
from deconv_bulk import EFFECTOR, IDENTITY, zscore_rows  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "bulk" / "meta"

# Published CRP performance in axSpA, for comparison.
#   sens 0.35 / spec 0.90 at the standard 0.5 mg/dL cutoff (disease activity)
#   AUROC 0.71 for high disease activity
CRP = dict(name="CRP (published)", sens=0.35, spec=0.90, auc=0.71)
# HLA-B27: carried by ~90% of AS patients and ~8% of the (European) population.
# The existing test, and the bar any new blood test has to clear.
B27 = dict(name="HLA-B27 (published)", sens=0.90, spec=0.92, auc=0.91)

# Prevalence of axSpA in the populations a test would actually be used on.
SETTINGS = {
    "general population": 0.005,
    "chronic back pain clinic": 0.05,
    "back pain + inflammatory features": 0.14,
    "our studies (artificial)": 97 / 167,
}


def load(gse: str):
    expr = pd.read_csv(PREP / f"{gse}_expr.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(PREP / f"{gse}_samples.tsv", sep="\t", index_col=0)
    expr = expr[[c for c in expr.columns if c in meta.index]]
    return expr, meta.loc[expr.columns]


def ppv_npv(sens: float, spec: float, prev: float) -> tuple[float, float]:
    """Bayes. The whole point: the same test has a different meaning in a
    different population."""
    tp = sens * prev
    fn = (1 - sens) * prev
    fp = (1 - spec) * (1 - prev)
    tn = spec * (1 - prev)
    return tp / (tp + fp), tn / (tn + fn)


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    print("=" * 94)
    print("STEP 1 -- how well does the NK effector score separate AS from healthy?")
    print("=" * 94)

    rows = []
    for gse in PRIMARY:
        expr, meta = load(gse)
        z = zscore_rows(expr)
        eff = [g for g in EFFECTOR if g in z.index]
        ident = [g for g in IDENTITY if g in z.index]
        y = (meta["group"] == "AS").to_numpy().astype(int)

        for name, genes in (("effector", eff), ("identity", ident),
                            ("both", sorted(set(eff) | set(ident)))):
            if len(genes) < 4:
                continue
            s = z.loc[genes].mean(axis=0)
            auc = roc_auc_score(y, -s)
            fpr, tpr, thr = roc_curve(y, -s)
            j = int(np.argmax(tpr - fpr))  # Youden: best single cutoff
            rows.append(dict(gse=gse, score=name, n_genes=len(genes), auc=auc,
                             sens=tpr[j], spec=1 - fpr[j],
                             n_as=int(y.sum()), n_hc=int((1 - y).sum())))

    df = pd.DataFrame(rows)
    for name, g in df.groupby("score", sort=False):
        w = g["n_as"] + g["n_hc"]
        print(f"\n  {name} score")
        for r in g.itertuples():
            print(f"    {r.gse:10s} {r.n_genes:2d} genes  AUC {r.auc:.3f}  "
                  f"sens {r.sens:.2f}  spec {r.spec:.2f}")
        print(f"    {'WEIGHTED':10s}            AUC {np.average(g['auc'], weights=w):.3f}")

    df.to_csv(OUT / "clinical_benchmark.tsv", sep="\t", index=False)

    best = df[df["score"] == "effector"]
    w = best["n_as"] + best["n_hc"]
    auc = float(np.average(best["auc"], weights=w))
    sens = float(np.average(best["sens"], weights=w))
    spec = float(np.average(best["spec"], weights=w))

    print("\n" + "=" * 94)
    print("STEP 2 -- what that means for one patient, by setting")
    print("=" * 94)
    tests = [
        ("RNA (NK effector)", sens, spec, auc),
        ("CRP", CRP["sens"], CRP["spec"], CRP["auc"]),
        ("HLA-B27", B27["sens"], B27["spec"], B27["auc"]),
    ]
    for nm, se, sp, au in tests:
        print(f"  {nm:20s} AUC {au:.2f}  sens {se:.2f}  spec {sp:.2f}")

    print(f"\n  {'setting':34s} {'prev':>6s} | " +
          " | ".join(f"{nm:>17s}" for nm, *_ in tests))
    print("  " + "-" * 110)
    for name, prev in SETTINGS.items():
        cells = []
        for _, se, sp, _ in tests:
            p, n = ppv_npv(se, sp, prev)
            cells.append(f"PPV {p:5.1%} NPV {n:5.1%}")
        print(f"  {name:34s} {prev:6.1%} | " + " | ".join(cells))

    print("\n  PPV = of everyone the test calls positive, the share who truly have AS.")
    print("  NPV = of everyone it calls negative, the share who truly do not.")

    print("\n" + "=" * 94)
    print("STEP 3 -- how many people you would have to test to find one case")
    print("=" * 94)
    print(f"  {'setting':34s} " + " ".join(f"{nm:>20s}" for nm, *_ in tests))
    for name, prev in SETTINGS.items():
        cells = []
        for _, se, sp, _ in tests:
            p, _ = ppv_npv(se, sp, prev)
            cells.append(f"{(1 / p - 1) if p > 0 else np.inf:20.1f}")
        print(f"  {name:34s} " + " ".join(cells))
    print("\n  (false alarms per true case found)")


if __name__ == "__main__":
    main()
