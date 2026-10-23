"""What did TRUST4 actually recover, and is the AS motif there?

Runs locally on the downloaded reports -- no compute cost.

The honest purpose of this is to measure RECOVERY, not to test the hypothesis.
With 5 AS patients and Faham's ~35% motif prevalence you would expect ~1.75
positives even with a perfect assay, so neither a hit nor a miss settles
anything. What it does settle: how many clonotypes bulk RNA-seq yields, and
how expanded the largest clones are -- which is what determines whether a
bigger run is worth paying for.

Usage:  python scripts/tcr_analyse.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TCR = ROOT / "data" / "tcr" / "GSE181364"

# run -> sample, from the ENA run table
SAMPLES = {
    "SRR15334651": "AS_patient_1", "SRR15334652": "AS_patient_2",
    "SRR15334653": "AS_patient_3", "SRR15334654": "AS_patient_4",
    "SRR15334655": "AS_patient_5", "SRR15334656": "Healthy_1",
    "SRR15334657": "Healthy_2",    "SRR15334658": "Healthy_3",
}

# Faham's AS-associated motif; Komech's full clonotype is CASSVGLYSTDTQYF
MOTIF = "SVGLYSTDTQ"
# Looser relatives, to show what a motif search over a small repertoire returns
PARTIALS = ["VGLYST", "GLYSTDTQ", "YSTDTQ"]


def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    return df


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    files = sorted(TCR.glob("*_report.tsv"))
    print(f"{len(files)} samples recovered of {len(SAMPLES)} attempted\n")

    rows, all_trb = [], []
    for f in files:
        acc = f.stem.replace("_report", "")
        name = SAMPLES.get(acc, acc)
        df = load(f)
        if "CDR3aa" not in df.columns:
            print(f"{name}: unexpected columns {list(df.columns)[:8]}")
            continue

        df = df[df["CDR3aa"].notna()].copy()
        df["CDR3aa"] = df["CDR3aa"].astype(str)
        df["V"] = df["V"].astype(str)
        # TRUST4 reports BCR and TCR together; keep TCR beta
        trb = df[df["V"].str.startswith("TRBV")].copy()
        # productive only: no stop codons, no frameshift placeholders
        trb = trb[~trb["CDR3aa"].str.contains(r"[\*\?_]", regex=True, na=False)]
        trb = trb[trb["CDR3aa"].str.len() >= 8]

        tot = trb["count"].sum()
        trb["freq"] = trb["count"] / tot if tot else 0.0
        trb["sample"] = name
        all_trb.append(trb)

        v9 = trb["V"].str.startswith("TRBV9")
        hit = trb["CDR3aa"].str.contains(MOTIF, na=False)

        rows.append(dict(
            sample=name,
            group="AS" if name.startswith("AS") else "HC",
            all_clonotypes=len(df),
            trb_productive=len(trb),
            trb_reads=int(tot),
            top_clone_pct=100 * trb["freq"].max() if len(trb) else np.nan,
            top10_pct=100 * trb["freq"].nlargest(10).sum() if len(trb) else np.nan,
            n_trbv9=int(v9.sum()),
            trbv9_pct=100 * trb.loc[v9, "freq"].sum() if len(trb) else 0.0,
            motif_hits=int(hit.sum()),
        ))
        for _, r in trb[hit].iterrows():
            print(f"  *** MOTIF {name}: {r['CDR3aa']}  {r['V']}  "
                  f"count={r['count']}  freq={100*r['freq']:.3f}%")

    out = pd.DataFrame(rows).sort_values(["group", "sample"])
    print("=" * 104)
    print("RECOVERY -- what bulk RNA-seq actually yields")
    print("=" * 104)
    with pd.option_context("display.width", 220):
        print(out.round(3).to_string(index=False))

    print("\n  'top_clone_pct' is the largest single clone as a share of the TCRb")
    print("  repertoire. This is the number that decides whether shallow data can")
    print("  see expanded clones at all.")

    trb = pd.concat(all_trb, ignore_index=True)

    print("\n" + "=" * 104)
    print("THE MOTIF")
    print("=" * 104)
    n_hit = int(trb["CDR3aa"].str.contains(MOTIF, na=False).sum())
    print(f"  exact '{MOTIF}' across all samples: {n_hit} clonotypes")
    for p in PARTIALS:
        sub = trb[trb["CDR3aa"].str.contains(p, na=False)]
        if len(sub):
            print(f"  partial '{p}': {len(sub)} clonotypes")
            for _, r in sub.iterrows():
                print(f"      {r['sample']:14s} {r['CDR3aa']:22s} {r['V']:12s} "
                      f"count={r['count']}")
        else:
            print(f"  partial '{p}': 0")

    print("\n" + "=" * 104)
    print("TRBV9 USAGE (a gene-level number -- NOT the finding, shown for contrast)")
    print("=" * 104)
    for grp, g in out.groupby("group"):
        print(f"  {grp}: n={len(g)}  TRBV9 share of repertoire "
              f"median {g['trbv9_pct'].median():.2f}%  "
              f"(range {g['trbv9_pct'].min():.2f}-{g['trbv9_pct'].max():.2f}%)")
    print("\n  Faham's result is a CDR3 SEQUENCE on TRBV9, not TRBV9 abundance.")
    print("  TRBV9 is a common gene used by everyone; its usage level is exactly")
    print("  the kind of graded, overlapping measurement that cannot discriminate.")

    print("\n" + "=" * 104)
    print("DETECTION LIMIT")
    print("=" * 104)
    med_reads = out["trb_reads"].median()
    print(f"  median TCRb reads assembled per sample: {med_reads:,.0f}")
    for f_ in (0.05, 0.01, 0.005, 0.001):
        print(f"    a clone at {100*f_:5.2f}% of the repertoire -> "
              f"~{med_reads * f_:6.1f} reads expected  "
              f"{'(detectable)' if med_reads * f_ >= 3 else '(missed)'}")

    out.to_csv(ROOT / "data" / "tcr" / "recovery_summary.tsv", sep="\t", index=False)
    trb.to_csv(ROOT / "data" / "tcr" / "all_trb_clonotypes.tsv.gz",
               sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
