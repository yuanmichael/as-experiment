"""Patient-held-out composition/rank pilot. No remote compute or supervised clone filter.

prepare writes immutable-by-convention features/manifests; evaluate writes OOF
predictions, fold membership, nested choices, and a rerun-tuning permutation null.
This coarse pooled representation does NOT replace the planned attention model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from threadpoolctl import threadpool_limits

TARGET = "CASSVGLYSTDTQYF"
SEED = 20261027
ALPHAS = (0.1, 1.0, 10.0)
PRIMARY = ("uniform", "frequency", "rank", "combined", "combined_rbf")


def write_json(path, value):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    tmp.replace(path)


def rank_weights(counts):
    # Average percentile for ties; bounded weights avoid arbitrary tie breaking.
    r = rankdata(counts, method="average") / len(counts)
    return r / r.sum()


def pooled(mat, counts):
    return np.stack([
        np.asarray(mat.mean(axis=0)).ravel(),
        np.asarray(mat.T @ (counts / counts.sum())).ravel(),
        np.asarray(mat.T @ rank_weights(counts)).ravel(),
    ])


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "features.npz").exists():
        raise ValueError("Features already exist; use another output directory to change inputs.")
    runs = {}
    for path in sorted(Path(args.tables).rglob("*_cdr3.json")):
        m = json.loads(path.read_text())
        if m.get("status") != "ok" or m.get("cell_type") != "PBMC":
            continue
        if m.get("disease") not in ("as", "healthy"):
            continue
        table = path.with_name(path.name.replace(".json", ".tsv.gz"))
        if not table.exists():
            raise ValueError(f"Missing table for {path}")
        sha = digest(table)
        if m["run"] in runs:
            if runs[m["run"]]["sha256"] != sha:
                raise ValueError(f"Conflicting duplicate run {m['run']}")
            continue
        runs[m["run"]] = dict(m, source=str(table.resolve()), sha256=sha)
    manifest = pd.DataFrame(runs.values()).sort_values(["subject", "run"])
    if manifest.empty:
        raise ValueError("No eligible input runs")
    if manifest.groupby("subject").disease.nunique().max() != 1:
        raise ValueError("Conflicting subject labels")
    manifest.to_csv(out / "runs.tsv", sep="\t", index=False)
    # Multiple BioSamples are not silently pooled. Largest recovered run is a
    # label-independent choice; record exclusions at the fixed depth.
    chosen = manifest.sort_values(["junction_reads", "run"], ascending=[False, True])
    chosen = chosen.drop_duplicates("subject").sort_values("subject")
    vectorizer = HashingVectorizer(analyzer="char", ngram_range=(3, 3),
                                  lowercase=False, n_features=512,
                                  alternate_sign=False, norm="l1")
    raw_meta = pd.read_csv(args.run_metadata, sep="\t").set_index("run_accession")
    arrays, removed, shuffled, tech, known, patients, excluded = [], [], [], [], [], [], []
    for m in chosen.to_dict("records"):
        df = pd.read_csv(m["source"], sep="\t")
        valid = df.cdr3.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]{10,}")
        df = df[valid].groupby("cdr3", as_index=False).reads.sum().sort_values("cdr3")
        counts = df.reads.to_numpy(dtype=np.int64)
        if (counts <= 0).any():
            raise ValueError("Non-positive counts")
        depth = int(counts.sum())
        if depth < args.depth:
            excluded.append(dict(subject=m["subject"], disease=m["disease"],
                                 depth=depth, reason="single selected run below fixed depth"))
            continue
        subject_seed = int.from_bytes(hashlib.sha256(m["subject"].encode()).digest()[:4], "little")
        rng = np.random.default_rng(SEED + subject_seed)
        thinned = rng.multivariate_hypergeometric(counts, args.depth)
        keep = thinned > 0
        counts = thinned[keep]
        seq = df.cdr3.to_numpy()[keep]
        mat = vectorizer.transform(seq)
        arrays.append(pooled(mat, counts))
        without = seq != TARGET
        removed.append(pooled(mat[without], counts[without]))
        shuffled.append(pooled(mat, rng.permutation(counts)))
        target_count = int(counts[~without].sum())
        known.append(target_count / args.depth)
        run = raw_meta.loc[m["run"]]
        mean_mate_length = float(run.base_count / run.read_count / 2)
        tech.append([np.log1p(depth), mean_mate_length])
        patients.append(dict(subject=m["subject"], run=m["run"], disease=m["disease"],
                             y=int(m["disease"] == "as"), raw_depth=depth,
                             thinned_depth=args.depth, clones=len(seq),
                             target_count=target_count, mean_mate_length=mean_mate_length,
                             candidate_runs=int((manifest.subject == m["subject"]).sum()),
                             replicate_policy="largest single run; no cross-BioSample merge"))
        print(f"prepared {m['subject']}: {len(seq):,} clones", flush=True)
    pd.DataFrame(patients).to_csv(out / "patients.tsv", sep="\t", index=False)
    pd.DataFrame(excluded).to_csv(out / "excluded.tsv", sep="\t", index=False)
    np.savez_compressed(out / "features.npz", composition=np.array(arrays),
                        without_target=np.array(removed), coupling_shuffled=np.array(shuffled),
                        technical=np.array(tech), known=np.array(known))
    write_json(out / "feature_config.json", dict(seed=SEED, depth=args.depth,
        n_hash_features=512, kmer_length=3, n_patients=len(patients), excluded=excluded,
        code_sha256=digest(__file__), numpy=np.__version__, pandas=pd.__version__))


def kernel_pair(xtr, xte, rbf=False):
    mu, sd = xtr.mean(axis=0), xtr.std(axis=0)
    # Constant training features are excluded, including unseen test-only features.
    active = sd > 1e-12
    if not active.any():
        return np.ones((len(xtr), len(xtr))), np.ones((len(xte), len(xtr)))
    a = (xtr[:, active] - mu[active]) / sd[active]
    b = (xte[:, active] - mu[active]) / sd[active]
    p = a.shape[1]
    ka, kb = a @ a.T / p, b @ a.T / p
    if rbf:
        aa, bb = (a*a).sum(axis=1) / p, (b*b).sum(axis=1) / p
        ka = np.exp(-np.maximum(aa[:, None]+aa[None, :]-2*ka, 0)/2)
        kb = np.exp(-np.maximum(bb[:, None]+aa[None, :]-2*kb, 0)/2)
    return ka, kb


def solve(ka, kb, y, alpha):
    center = y.mean()
    return kb @ np.linalg.solve(ka + alpha * np.eye(len(y)), y - center) + center


def nested_scores(x, y, rbf, folds, record=False):
    oof = np.zeros((len(folds), len(y)))
    aucs, choices, membership = [], [], []
    for rep, splits in enumerate(folds):
        for fold, (tr, te) in enumerate(splits):
            inner = StratifiedKFold(3, shuffle=True, random_state=SEED + rep * 5 + fold)
            val = {alpha: [] for alpha in ALPHAS}
            for it, iv in inner.split(x[tr], y[tr]):
                ka, kb = kernel_pair(x[tr[it]], x[tr[iv]], rbf)
                for alpha in ALPHAS:
                    val[alpha].append(roc_auc_score(y[tr[iv]], solve(ka, kb, y[tr[it]], alpha)))
            # Strongest regularization wins ties; grid is fixed before results.
            alpha = max(ALPHAS, key=lambda a: (np.mean(val[a]), a))
            ka, kb = kernel_pair(x[tr], x[te], rbf)
            pred = solve(ka, kb, y[tr], alpha)
            oof[rep, te] = pred
            aucs.append(roc_auc_score(y[te], pred))
            if record:
                choices.append(dict(repeat=rep, fold=fold, alpha=alpha,
                                    auc=aucs[-1], inner_auc=float(np.mean(val[alpha]))))
                membership.append(dict(repeat=rep, fold=fold,
                                       train=tr.tolist(), test=te.tolist()))
    return float(np.mean(aucs)), oof, choices, membership


def evaluate(args):
    out = Path(args.out)
    result_path = out / "results.json"
    if result_path.exists():
        raise ValueError("Results already exist; use a separate run directory for changes")
    p = pd.read_csv(out / "patients.tsv", sep="\t")
    f = np.load(out / "features.npz")
    c = f["composition"]
    y = p.y.to_numpy()
    models = {"uniform": (c[:, 0], False), "frequency": (c[:, 1], False),
              "rank": (c[:, 2], False), "combined": (c.reshape(len(y), -1), False),
              "combined_rbf": (c.reshape(len(y), -1), True),
              "without_target": (f["without_target"].reshape(len(y), -1), False),
              "coupling_shuffled": (f["coupling_shuffled"].reshape(len(y), -1), False),
              "technical": (f["technical"], False)}
    def splits(labels):
        return [list(StratifiedKFold(5, shuffle=True, random_state=SEED+r).split(y, labels))
                for r in range(args.repeats)]
    folds = splits(y)
    observed, preds = {}, p[["subject", "y"]].copy()
    for name, (x, rbf) in models.items():
        auc, oof, choices, membership = nested_scores(x, y, rbf, folds, record=True)
        observed[name] = dict(mean_fold_auc=auc, pooled_oof_auc=float(roc_auc_score(y, oof.mean(0))),
                              choices=choices)
        for rep in range(args.repeats):
            preds[f"{name}_repeat{rep}"] = oof[rep]
        write_json(out / "folds.json", membership)
        print(f"{name}: mean-fold AUC {auc:.3f}", flush=True)
    known_auc = np.mean([roc_auc_score(y[te], f["known"][te]) for s in folds for _, te in s])
    observed["known_clone"] = dict(mean_fold_auc=float(known_auc),
                                  pooled_auc=float(roc_auc_score(y, f["known"])))
    preds["known_clone"] = f["known"]
    preds.to_csv(out / "predictions.tsv", sep="\t", index=False)
    write_json(out / "observed.json", observed)
    rng = np.random.default_rng(SEED + 1)
    nulls = []
    for b in range(args.permutations):
        yp = rng.permutation(y)
        perm_folds = splits(yp)
        row = {name: nested_scores(x, yp, rbf, perm_folds)[0]
               for name, (x, rbf) in models.items()}
        row["max_primary"] = max(row[n] for n in PRIMARY)
        row["max_all"] = max(row[n] for n in models)
        nulls.append(row)
        if (b+1) % 10 == 0 or b+1 == args.permutations:
            pd.DataFrame(nulls).to_csv(out / "permutations.tsv", sep="\t", index=False)
            print(f"permutations {b+1}/{args.permutations}", flush=True)
    null = pd.DataFrame(nulls)
    summary = {}
    for name in models:
        auc = observed[name]["mean_fold_auc"]
        summary[name] = {k:v for k,v in observed[name].items() if k != "choices"}
        if args.permutations:
            summary[name].update(p_unadjusted=float((1+(null[name]>=auc).sum())/(len(null)+1)),
                p_max_all=float((1+(null.max_all>=auc).sum())/(len(null)+1)))
            if name in PRIMARY:
                summary[name]["p_max_primary"] = float((1+(null.max_primary>=auc).sum())/(len(null)+1))
    summary["known_clone"] = observed["known_clone"]
    write_json(result_path, dict(n=len(y), n_as=int(y.sum()), n_healthy=int((y==0).sum()),
        repeats=args.repeats, permutations=args.permutations, models=summary,
        code_sha256=digest(__file__), feature_sha256=digest(out / "features.npz"),
        limitations=["exploratory within-study; no external validation performed",
          "CDR3 3-mer pooling only; no full embeddings, attention or V/J model",
          "fixed read thinning; reads are not independent biological cells",
          "global label null; not proof of AS-specific or conditional incremental information",
          "single count-shuffle seed; no patient bootstrap uncertainty yet"]))
    print(json.dumps(summary, indent=2), flush=True)


def self_test():
    from scipy.sparse import csr_matrix
    counts = np.array([1, 1, 8])
    assert rank_weights(counts)[0] == rank_weights(counts)[1]
    a = pooled(csr_matrix(np.eye(3)), counts)
    assert np.allclose(a.sum(1), 1)
    assert np.allclose(a[1], counts/counts.sum())
    # Row permutation leaves all three pooled views unchanged.
    order = [2, 0, 1]
    assert np.allclose(a, pooled(csr_matrix(np.eye(3)[order]), counts[order]))
    # Strong distributed signal with disjoint held-out patients is recoverable.
    rng = np.random.default_rng(123)
    y = np.tile([0, 1], 20)
    x = rng.normal(size=(40, 32)) + y[:, None] * 0.9
    folds = [list(StratifiedKFold(5, shuffle=True, random_state=7).split(x, y))]
    assert all(not set(tr) & set(te) for tr, te in folds[0])
    assert nested_scores(x, y, False, folds)[0] > 0.9
    # Test samples never affect training standardization.
    ka, _ = kernel_pair(x[:20], x[20:])
    kb, _ = kernel_pair(x[:20], x[20:] * 100)
    assert np.allclose(ka, kb)
    print("self-test passed", flush=True)


def external(args):
    """Score every prespecified model, without using external labels for fitting."""
    out = Path(args.out)
    if (out / "external_results.json").exists():
        raise ValueError("External results already exist; do not overwrite after inspection")
    p = pd.read_csv(out / "patients.tsv", sep="\t")
    f = np.load(out / "features.npz")
    y = p.y.to_numpy()
    df = pd.read_parquet(args.external_table)
    df = df[(df.tissue == "blood") & df.group.isin(["AS", "HC"])]
    if df.groupby("subject").group.nunique().max() != 1:
        raise ValueError("Conflicting external labels")
    vec = HashingVectorizer(analyzer="char", ngram_range=(3, 3), lowercase=False,
                           n_features=512, alternate_sign=False, norm="l1")
    rows, xs, known = [], [], []
    for subject, t in df.groupby("subject", sort=True):
        if t["sample"].nunique() != 1:
            raise ValueError("External person has multiple blood samples; resolve visits")
        label = int(t.group.iloc[0] == "AS")
        t = t[t.cdr3.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]{10,}")]
        t = t.groupby("cdr3", as_index=False).cells.sum()
        counts = t.cells.to_numpy()
        xs.append(pooled(vec.transform(t.cdr3), counts))
        known.append(float(t.loc[t.cdr3 == TARGET, "cells"].sum() / counts.sum()))
        rows.append(dict(subject=subject, y=label, counted_cells=int(counts.sum()),
                         distinct_cdr3=len(t)))
    ext = pd.DataFrame(rows)
    c, e = f["composition"], np.array(xs)
    matrices = {
        "uniform": (c[:, 0], e[:, 0], False),
        "frequency": (c[:, 1], e[:, 1], False),
        "rank": (c[:, 2], e[:, 2], False),
        "combined": (c.reshape(len(y), -1), e.reshape(len(ext), -1), False),
        "combined_rbf": (c.reshape(len(y), -1), e.reshape(len(ext), -1), True),
    }
    summaries = {}
    cv = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(y, y))
    for name, (x, xe, rbf) in matrices.items():
        scores = {a: [] for a in ALPHAS}
        for tr, te in cv:
            ka, kb = kernel_pair(x[tr], x[te], rbf)
            for alpha in ALPHAS:
                scores[alpha].append(roc_auc_score(y[te], solve(ka, kb, y[tr], alpha)))
        alpha = max(ALPHAS, key=lambda a: (np.mean(scores[a]), a))
        ka, kb = kernel_pair(x, xe, rbf)
        ext[name] = solve(ka, kb, y, alpha)
        summaries[name] = dict(alpha=alpha, auc=float(roc_auc_score(ext.y, ext[name])))
    ext["known_clone"] = known
    summaries["known_clone"] = dict(auc=float(roc_auc_score(ext.y, known)))
    ext.to_csv(out / "external_predictions.tsv", sep="\t", index=False)
    write_json(out / "external_results.json", dict(study="GSE288581", n=len(ext),
        models=summaries, input_sha256=digest(args.external_table),
        code_sha256=digest(__file__),
        limitations=["previously inspected cohort; assay and sorting shift",
                     "four AS and four healthy; descriptive transfer check",
                     "native cell counts versus thinned junction reads",
                     "all model choices/tuning use Komech only; no external calibration"]))
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["prepare", "evaluate", "external", "self-test"])
    ap.add_argument("--tables")
    ap.add_argument("--out", default="data/repertoire_experiment")
    ap.add_argument("--run-metadata", default="data/airr/komech_runs.tsv")
    ap.add_argument("--depth", type=int, default=800000)
    ap.add_argument("--permutations", type=int, default=199)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--external-table", default="data/as_vdj/clonotypes.parquet")
    args = ap.parse_args()
    with threadpool_limits(limits=1):
        if args.mode == "prepare":
            if not args.tables:
                ap.error("prepare requires --tables")
            prepare(args)
        elif args.mode == "evaluate":
            evaluate(args)
        elif args.mode == "external":
            external(args)
        else:
            self_test()
