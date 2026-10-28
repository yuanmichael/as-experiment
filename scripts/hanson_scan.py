"""Scan the Hanson 2020 immuneACCESS export for the AS motif family.

The export (193 Adaptive v1 sample files) carries no diagnosis column, so group
membership has to be inferred from the sample names. The naming splits into
`A####` and `K####` subject prefixes with `_4` / `_8` suffixes for the sorted
CD4 and CD8 fractions. The hypothesis under test here is A = AS patient,
K = healthy control.

That hypothesis is checked against a published number rather than assumed:
Hanson reports `CASSVGLFSTDTQYF` in **12 of 37 HLA-B27+ AS patients and absent
from every B27- patient and every control** in the CD8 compartment. If the
A/K split reproduces that count, the mapping is confirmed by an independent
statistic and the pipeline is validated at the same time. If it does not, the
labels are wrong and nothing downstream should be computed.

Note what this cannot do: once labels are inferred from these data, a p-value
computed on the same data is conditioned on the inference. The confirmation
above is what makes it usable, and the authors should still be asked for the
real key before anything is published.

    python scripts/hanson_scan.py --export "<dir of 193 .tsv>"
"""
from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict

import pandas as pd

from motif_family_test import EXPLORATORY, PRESPECIFIED

PANEL = list(PRESPECIFIED) + list(EXPLORATORY)
HANSON_CLONE = "CASSVGLFSTDTQYF"
NAME = re.compile(r"^(?P<sub>.+?)_(?P<tag>[NT]?)(?P<frac>[48])(?P<rest>.*)$")


def parse_name(stem: str) -> dict | None:
    m = NAME.match(stem)
    if not m:
        return None
    sub = m.group("sub")
    prefix = ("A" if re.match(r"^A\d", sub) else
              "K" if re.match(r"^K\d", sub) else
              "digits" if re.match(r"^\d", sub) else "other")
    return dict(sample=stem, subject=sub, prefix=prefix, tag=m.group("tag") or "-",
                fraction="CD4" if m.group("frac") == "4" else "CD8",
                replicate=bool(m.group("rest")))


def scan(path: str, targets: set[str]) -> tuple[dict[str, int], int]:
    """Templates per target sequence, plus the file's total productive templates."""
    found: dict[str, int] = {}
    total = 0
    cols = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_aa = header.index("aminoAcid")
            i_n = header.index("count (templates/reads)")
            i_st = header.index("sequenceStatus")
        except ValueError:
            return found, 0
        cols = (i_aa, i_n, i_st)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(cols):
                continue
            if parts[i_st] != "In":
                continue
            try:
                n = int(parts[i_n])
            except ValueError:
                continue
            total += n
            aa = parts[i_aa]
            if aa in targets:
                found[aa] = found.get(aa, 0) + n
    return found, total


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", required=True)
    ap.add_argument("--out", default="data/panel/hanson_family_scan.tsv")
    ap.add_argument("--min-templates", type=int, default=1,
                    help="templates needed to call a sequence present")
    args = ap.parse_args()

    targets = set(PANEL)
    rows = []
    files = sorted(f for f in os.listdir(args.export) if f.endswith(".tsv"))
    print(f"scanning {len(files)} files for {len(targets)} sequences")
    for k, fn in enumerate(files, 1):
        meta = parse_name(os.path.splitext(fn)[0])
        if meta is None:
            meta = dict(sample=os.path.splitext(fn)[0], subject=os.path.splitext(fn)[0],
                        prefix="other", tag="-", fraction="?", replicate=False)
        found, total = scan(os.path.join(args.export, fn), targets)
        row = dict(meta, total_templates=total)
        for seq in PANEL:
            row[seq] = found.get(seq, 0)
        rows.append(row)
        if k % 40 == 0:
            print(f"  {k}/{len(files)}")

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print("\nFILES BY GROUP HYPOTHESIS")
    print(df.groupby(["prefix", "fraction"]).size().unstack(fill_value=0))

    thr = args.min_templates
    print(f"\nCARRIERS BY PREFIX AND FRACTION (>= {thr} template)")
    for frac in ("CD8", "CD4"):
        sub = df[df["fraction"] == frac]
        if sub.empty:
            continue
        print(f"\n  --- {frac} ---")
        hdr = f"  {'cdr3':<22}" + "".join(f"{p:>14}" for p in sorted(sub['prefix'].unique()))
        print(hdr)
        for seq in PANEL:
            cells = ""
            for p in sorted(sub["prefix"].unique()):
                s = sub[sub["prefix"] == p]
                carriers = s.loc[s[seq] >= thr, "subject"].nunique()
                n = s["subject"].nunique()
                cells += f"{carriers:>6}/{n:<7}"
            print(f"  {seq:<22}{cells}")
        # family union, per subject
        cells = ""
        for p in sorted(sub["prefix"].unique()):
            s = sub[sub["prefix"] == p]
            carrier_subs = set()
            for seq in PRESPECIFIED:
                carrier_subs |= set(s.loc[s[seq] >= thr, "subject"])
            cells += f"{len(carrier_subs):>6}/{s['subject'].nunique():<7}"
        print(f"  {'UNION (7 prespecified)':<22}{cells}")

    hit = df[(df["fraction"] == "CD8") & (df[HANSON_CLONE] >= thr)]
    print(f"\nPUBLISHED-COUNT CHECK — Hanson reports {HANSON_CLONE} in 12 of 37 "
          f"B27+ AS CD8 samples, 0 controls")
    print(f"  observed: {hit['subject'].nunique()} subjects carry it in CD8")
    print("  by prefix: " + ", ".join(
        f"{p}={g['subject'].nunique()}" for p, g in hit.groupby("prefix")))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
