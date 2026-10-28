"""Test Komech's untestable candidates in Hanson's control arm.

The original sweep found 1,075 clonotypes carried by 6-11 of 25 Komech AS
patients and none of its 10 healthy controls, and could not call any of them
significant: with 25 cases and 10 controls the smallest reachable p-value for a
marker in k cases and 0 controls does not clear 0.01 until k >= 12. The
candidates were not weak, the control arm was too small to measure them.

Hanson 2020 supplies 36 controls, which drops that floor to 8 of 48 AS carriers
(17%). This script takes the candidate list from Komech and evaluates it in
Hanson, so discovery and evaluation sit in different cohorts, on different
platforms, with different control donors. Nothing about Hanson's labels informs
which sequences are tested.

Three built-in controls:
  * the CD4 fraction, where a CD8-restricted signal must be absent;
  * the motif family, which must come out positive or the pipeline is broken;
  * a size-matched random draw from clonotypes seen in Komech AS at any
    carrier count, to show the candidate set beats an arbitrary one.

    python scripts/komech_candidates_in_hanson.py \
        --tables <komech cdr3 dir> --export <hanson export dir>
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clono_sweep import carrier_matrix, load_tables, thin  # noqa: E402
from hanson_scan import parse_name, scan  # noqa: E402
from motif_family_test import PRESPECIFIED  # noqa: E402

THIN_DEPTH = 800_000
SEED = 0


def komech_candidates(tables: str, min_car: int, max_car: int, case: str):
    rng = np.random.default_rng(SEED)
    df = load_tables(tables)
    df = df[df["group"].isin([case, "healthy"])]
    thinned = thin(df, THIN_DEPTH, rng)
    mat, subjects, cdr3s = carrier_matrix(thinned, "presence")
    grp = (thinned.drop_duplicates("subject").set_index("subject")["group"]
           .reindex(subjects).to_numpy())
    is_case = grp == case
    a = np.asarray(is_case.astype(np.int32) @ mat).ravel()
    b = np.asarray((~is_case).astype(np.int32) @ mat).ravel()
    sel = (a >= min_car) & (a <= max_car) & (b == 0)
    print(f"Komech: {int(is_case.sum())} {case}, {int((~is_case).sum())} healthy; "
          f"{mat.shape[1]:,} clonotypes after thinning")
    print(f"candidates in {min_car}-{max_car} AS patients and 0 controls: "
          f"{int(sel.sum()):,}")
    pool = cdr3s[(a >= 1) & (b == 0)]
    return cdr3s[sel], a[sel], pool


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--export", required=True)
    ap.add_argument("--min-carriers", type=int, default=6)
    ap.add_argument("--max-carriers", type=int, default=11)
    ap.add_argument("--case", default="as")
    ap.add_argument("--out", default="data/panel/komech_candidates_in_hanson.tsv")
    args = ap.parse_args()

    cands, komech_car, pool = komech_candidates(
        args.tables, args.min_carriers, args.max_carriers, args.case)
    rng = np.random.default_rng(SEED)
    decoys = rng.choice(pool[~np.isin(pool, cands)], size=len(cands), replace=False)
    targets = set(cands) | set(decoys) | set(PRESPECIFIED)
    print(f"scanning Hanson for {len(targets):,} sequences "
          f"({len(cands):,} candidates + {len(decoys):,} decoys + family)")

    rows = []
    files = sorted(f for f in os.listdir(args.export) if f.endswith(".tsv"))
    for k, fn in enumerate(files, 1):
        meta = parse_name(os.path.splitext(fn)[0])
        if meta is None or meta["prefix"] not in ("A", "K"):
            continue
        found, total = scan(os.path.join(args.export, fn), targets)
        rows.append((meta["prefix"], meta["fraction"], meta["subject"], found))
        if k % 50 == 0:
            print(f"  {k}/{len(files)}")

    def table(fraction: str):
        """carriers[seq] -> (AS carriers, control carriers), plus group sizes."""
        subs = {"A": set(), "K": set()}
        car = {"A": {}, "K": {}}
        for prefix, frac, subject, found in rows:
            if frac != fraction:
                continue
            subs[prefix].add(subject)
            for seq in found:
                car[prefix].setdefault(seq, set()).add(subject)
        return car, len(subs["A"]), len(subs["K"])

    out = []
    car8, nA, nK = table("CD8")
    car4, nA4, nK4 = table("CD4")
    print(f"\nHanson CD8: {nA} AS vs {nK} control;  CD4: {nA4} vs {nK4}")

    for seq in sorted(targets):
        a = len(car8.get("A", {}).get(seq, ()))
        b = len(car8.get("K", {}).get(seq, ()))
        kind = ("family" if seq in PRESPECIFIED else
                "candidate" if seq in set(cands) else "decoy")
        out.append(dict(
            cdr3=seq, kind=kind,
            komech_as=int(komech_car[list(cands).index(seq)]) if kind == "candidate" else 0,
            hanson_as=a, hanson_ctrl=b,
            cd4_as=len(car4.get("A", {}).get(seq, ())),
            p=fisher_exact([[a, nA - a], [b, nK - b]], alternative="greater")[1]))
    df = pd.DataFrame(out)

    print("\nREPLICATION IN HANSON CD8")
    print(f"{'set':<12}{'n':>7}{'median AS carriers':>22}{'>=8 carriers, 0 ctrl':>24}"
          f"{'p<=0.01':>10}{'FDR<=0.10':>11}")
    for kind in ("candidate", "decoy", "family"):
        s = df[df["kind"] == kind]
        if s.empty:
            continue
        q = false_discovery_control(s["p"].to_numpy(), method="bh")
        strong = int(((s["hanson_as"] >= 8) & (s["hanson_ctrl"] == 0)).sum())
        print(f"{kind:<12}{len(s):>7}{s['hanson_as'].median():>22.1f}{strong:>24}"
              f"{int((s['p'] <= 0.01).sum()):>10}{int((q <= 0.10).sum()):>11}")

    cand = df[df["kind"] == "candidate"].copy()
    cand["q"] = false_discovery_control(cand["p"].to_numpy(), method="bh")
    cand = cand.sort_values("p")
    print("\nTOP CANDIDATES BY HANSON p")
    print(f"  {'cdr3':<20}{'Komech AS/25':>14}{'Hanson AS/48':>14}"
          f"{'ctrl/36':>9}{'CD4':>6}{'p':>10}{'q':>9}")
    for _, r in cand.head(15).iterrows():
        print(f"  {r['cdr3']:<20}{r['komech_as']:>14}{r['hanson_as']:>14}"
              f"{r['hanson_ctrl']:>9}{r['cd4_as']:>6}{r['p']:>10.3g}{r['q']:>9.3g}")

    df.to_csv(args.out, sep="\t", index=False, float_format="%.6g")

    # Burden: panel defined in Komech, scored per Hanson subject. No Hanson
    # label informs membership, so this is the leak-free version of the
    # burden score that was retracted when it was built and scored on Komech.
    from scipy.stats import mannwhitneyu
    sets = {"candidate": set(cands), "decoy": set(decoys),
            "family": set(PRESPECIFIED)}
    burden = []
    for prefix, frac, subject, found in rows:
        row = dict(prefix=prefix, fraction=frac, subject=subject)
        for name, members in sets.items():
            row[name] = len(members & set(found))
        burden.append(row)
    bdf = (pd.DataFrame(burden).groupby(["prefix", "fraction", "subject"])
           .max().reset_index())
    print("\nBURDEN PER SUBJECT (panel from Komech, scored in Hanson)")
    print(f"{'set':<12}{'fraction':<10}{'AS median (IQR)':>22}"
          f"{'ctrl median (IQR)':>22}{'AUC':>8}{'p':>11}")
    for frac in ("CD8", "CD4"):
        for name in sets:
            s = bdf[bdf["fraction"] == frac]
            a = s.loc[s["prefix"] == "A", name].to_numpy()
            k = s.loc[s["prefix"] == "K", name].to_numpy()
            if not len(a) or not len(k):
                continue
            u, p = mannwhitneyu(a, k, alternative="greater")
            fmt = lambda v: (f"{np.median(v):.0f} "
                             f"({np.percentile(v, 25):.0f}-{np.percentile(v, 75):.0f})")
            print(f"{name:<12}{frac:<10}{fmt(a):>22}{fmt(k):>22}"
                  f"{u/(len(a)*len(k)):>8.3f}{p:>11.3g}")
    bout = args.out.replace(".tsv", "_burden.tsv")
    bdf.to_csv(bout, sep="\t", index=False)
    print(f"\nwrote {args.out} and {bout}")


if __name__ == "__main__":
    main()
