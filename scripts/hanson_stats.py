"""Fisher tests for the motif family in the Hanson export, A (AS) vs K (control), CD8."""
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu

import sys
sys.path.insert(0, r"C:\Coding\as-experiment\scripts")
from motif_family_test import PRESPECIFIED, EXPLORATORY

df = pd.read_csv(r"C:\Coding\as-experiment\data\panel\hanson_family_scan.tsv", sep="\t")
cd8 = df[df["fraction"] == "CD8"]
panel = list(PRESPECIFIED) + list(EXPLORATORY)

A = cd8[cd8["prefix"] == "A"]
K = cd8[cd8["prefix"] == "K"]
nA, nK = A["subject"].nunique(), K["subject"].nunique()
print(f"CD8 compartment: {nA} A-prefix (AS) vs {nK} K-prefix (control) subjects\n")


def carriers(frame, seqs):
    s = set()
    for q in seqs:
        s |= set(frame.loc[frame[q] >= 1, "subject"])
    return s


def test(label, seqs):
    a, k = len(carriers(A, seqs)), len(carriers(K, seqs))
    tab = [[a, nA - a], [k, nK - k]]
    p1 = fisher_exact(tab, alternative="greater")[1]
    p2 = fisher_exact(tab)[1]
    print(f"{label:<30} AS {a:>2}/{nA}  ({a/nA:5.1%})   ctrl {k:>2}/{nK} "
          f"  sens {a/nA:5.1%}  spec {(nK-k)/nK:5.1%}   p = {p1:.3g}  (2-sided {p2:.3g})")
    return a, k


print("PER SEQUENCE")
for q in panel:
    tag = "  " if q in PRESPECIFIED else "* "
    test(tag + q, [q])

print("\nPRIMARY")
test("UNION of 7 pre-specified", list(PRESPECIFIED))
test("CASSVGLYSTDTQYF alone", ["CASSVGLYSTDTQYF"])
test("Hanson's Y+F pair", ["CASSVGLYSTDTQYF", "CASSVGLFSTDTQYF"])
test("UNION + 2 exploratory", panel)

cnt = {}
for frame, lab in ((A, "AS"), (K, "ctrl")):
    per = frame.groupby("subject")[list(PRESPECIFIED)].max()
    cnt[lab] = (per >= 1).sum(axis=1).to_numpy()
u, p = mannwhitneyu(cnt["AS"], cnt["ctrl"], alternative="greater")
print(f"\nCARRIER COUNT (0-7): AS median {np.median(cnt['AS']):.1f} range "
      f"{cnt['AS'].min()}-{cnt['AS'].max()}; ctrl median {np.median(cnt['ctrl']):.1f} "
      f"range {cnt['ctrl'].min()}-{cnt['ctrl'].max()}")
print(f"  Mann-Whitney p = {p:.4g}, AUC = {u/(len(cnt['AS'])*len(cnt['ctrl'])):.3f}")

print("\nDEPTH CHECK (total productive templates, CD8)")
for frame, lab in ((A, "AS"), (K, "ctrl")):
    t = frame.groupby("subject")["total_templates"].sum()
    print(f"  {lab:<5} median {t.median():,.0f}  range {t.min():,.0f}-{t.max():,.0f}")
tA = A.groupby("subject")["total_templates"].sum()
tK = K.groupby("subject")["total_templates"].sum()
print(f"  Mann-Whitney on depth: p = {mannwhitneyu(tA, tK).pvalue:.3f}")
