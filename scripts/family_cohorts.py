"""Motif-family carriage in PRJEB51315 (AS vs PsA) and GSE288581."""
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, r"C:\Coding\as-experiment\scripts")
from motif_family_test import PRESPECIFIED, EXPLORATORY

PANEL = list(PRESPECIFIED) + list(EXPLORATORY)
B = r"C:\Users\myuan\AppData\Local\Temp\claude\C--Coding-as-experiment\2f391c41-410f-452c-9331-4e63aba51115\scratchpad"

runs = pd.read_csv(r"C:\Coding\as-experiment\data\airr\komech2022_runs.tsv", sep="\t")
meta = {}
for _, r in runs.iterrows():
    p = str(r["library_name"]).split("_")
    meta[r["run_accession"]] = dict(dx=p[0], subject=p[1], tissue=p[2],
                                    frac="_".join(p[3:]))

rows = []
for d in ("cdr3", r"volaudit\cdr3sf"):
    for f in glob.glob(os.path.join(B, d, "*_cdr3.tsv.gz")):
        acc = os.path.basename(f).split("_")[0]
        if acc not in meta:
            continue
        t = pd.read_csv(f, sep="\t")
        tot = int(t["reads"].sum())
        hit = t[t["cdr3"].isin(PANEL)].groupby("cdr3")["reads"].sum().to_dict()
        rows.append(dict(run=acc, **meta[acc], total=tot,
                         **{q: hit.get(q, 0) for q in PANEL}))
df = pd.DataFrame(rows).drop_duplicates("run")
print(f"PRJEB51315: {len(df)} runs, {df['subject'].nunique()} subjects")
print(df.groupby(["dx", "tissue"])["subject"].nunique().to_string(), "\n")

MIN = 2  # reads, matching MIN_READS_PRESENT


def carriers(sub, seqs):
    s = set()
    for q in seqs:
        s |= set(sub.loc[sub[q] >= MIN, "subject"])
    return s


for tissue in ("PB", "SF"):
    cut = df[df["tissue"] == tissue]
    if cut.empty:
        continue
    print(f"=== PRJEB51315 {tissue} (>= {MIN} reads) ===")
    groups = sorted(cut["dx"].unique())
    print(f"  {'cdr3':<24}" + "".join(f"{g:>12}" for g in groups))
    for q in PANEL + ["UNION"]:
        seqs = list(PRESPECIFIED) if q == "UNION" else [q]
        cells = ""
        for g in groups:
            s = cut[cut["dx"] == g]
            cells += f"{len(carriers(s, seqs)):>5}/{s['subject'].nunique():<7}"
        tag = "UNION (7 prespec)" if q == "UNION" else (
            q if q in PRESPECIFIED else "* " + q)
        print(f"  {tag:<24}{cells}")
    print()

p = r"C:\Coding\as-experiment\data\as_vdj\clonotypes.parquet"
t = pd.read_parquet(p)
print("=== as_vdj/clonotypes.parquet ===")
print(t.groupby(["group", "tissue"])["subject"].nunique().to_string())
for (g, ti), s in t.groupby(["group", "tissue"]):
    hits = {}
    for q in PANEL:
        hits[q] = s.loc[s["cdr3"] == q, "subject"].nunique()
    u = set()
    for q in PRESPECIFIED:
        u |= set(s.loc[s["cdr3"] == q, "subject"])
    n = s["subject"].nunique()
    named = {k: v for k, v in hits.items() if v}
    print(f"  {g:<10} {ti:<10} n={n:<3} union={len(u)}/{n}  {named}")
