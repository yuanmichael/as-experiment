"""Descriptive demographic reference on the composition pilot's exact patient folds.

Added after cohort metadata audit; not part of the frozen eight-model null family.
No AS-specific or incremental significance claim is made by this comparison.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from repertoire_composition import nested_scores, write_json


def main():
    root = Path(__file__).resolve().parents[1] / "data/repertoire_experiment"
    if (root / "demographic_reference.json").exists():
        raise ValueError("Refusing to overwrite inspected demographic reference")
    patients = pd.read_csv(root / "patients.tsv", sep="\t")
    runs = pd.read_csv(root / "runs.tsv", sep="\t")
    p = patients.merge(runs[["run", "age", "sex"]], on="run", validate="one_to_one", sort=False)
    assert p.subject.tolist() == patients.subject.tolist()
    p["age"] = pd.to_numeric(p.age, errors="raise")
    assert p.sex.isin(["male", "female"]).all()
    p["male"] = (p.sex == "male").astype(float)
    rows = json.loads((root / "folds.json").read_text())
    folds = [[(np.array(f["train"]), np.array(f["test"]))
              for f in rows if f["repeat"] == r] for r in sorted({f["repeat"] for f in rows})]
    auc, pred, choices, _ = nested_scores(p[["age", "male"]].to_numpy(),
                                         p.y.to_numpy(), False, folds, record=True)
    out = p[["subject", "y", "age", "sex"]].copy()
    for r in range(len(folds)):
        out[f"prediction_repeat{r}"] = pred[r]
    out.to_csv(root / "demographic_predictions.tsv", sep="\t", index=False)
    result = dict(mean_fold_auc=auc, choices=choices,
                  interpretation="descriptive age/sex reference; no adjusted disease effect tested")
    write_json(root / "demographic_reference.json", result)
    print(f"Demographic reference mean-fold AUROC: {auc:.3f}")


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
