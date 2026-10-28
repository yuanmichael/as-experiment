"""Discrimination statistics Paley 2024 did not report, from its own figure source data."""
import sys
sys.path.insert(0, r"C:/Users/myuan/AppData/Local/Temp/claude/C--Coding-as-experiment/6750aed4-d886-4a0d-a46d-a64d428f2328/scratchpad/pylibs")
import openpyxl
import numpy as np
from scipy.stats import mannwhitneyu

XLSX = r"C:\Users\myuan\Downloads\jciinsight-9-174776-s254.xlsx"
wb = openpyxl.load_workbook(XLSX, data_only=True)


def cols(sheet):
    rows = list(wb[sheet].iter_rows(values_only=True))
    hdr, out = rows[0], {}
    for j, h in enumerate(hdr):
        if h is None or not str(h).strip():
            continue
        vals = []
        for r in rows[1:]:
            if j < len(r) and r[j] is not None and str(r[j]).strip() != "":
                try:
                    vals.append(float(r[j]))
                except ValueError:
                    pass
        out[str(h).strip()] = np.array(vals)
    return out


def auc(pos, neg):
    return mannwhitneyu(pos, neg, alternative="two-sided").statistic / (len(pos) * len(neg))


def summarise(label, sheet):
    c = cols(sheet)
    hc_key = [k for k in c if "ealthy" in k][0]
    pt_key = [k for k in c if "ealthy" not in k][0]
    hc, pt = c[hc_key], c[pt_key]
    rng = np.random.default_rng(0)
    bs = [auc(rng.choice(pt, len(pt), True), rng.choice(hc, len(hc), True)) for _ in range(4000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    p = mannwhitneyu(pt, hc, alternative="greater").pvalue
    best = (-1, None, 0, 0)
    for t in np.unique(np.concatenate([pt, hc])):
        se, sp = (pt >= t).mean(), (hc < t).mean()
        if se + sp > best[0]:
            best = (se + sp, t, se, sp)
    print(f"{label}")
    print(f"   n = {len(pt)} patients vs {len(hc)} B27+ healthy controls")
    print(f"   AUC = {auc(pt, hc):.3f}   95% CI {lo:.3f}-{hi:.3f}   Mann-Whitney p = {p:.4f}")
    print(f"   patients  median {np.median(pt):.4g}  IQR {np.percentile(pt, 25):.4g}-{np.percentile(pt, 75):.4g}")
    print(f"   controls  median {np.median(hc):.4g}  IQR {np.percentile(hc, 25):.4g}-{np.percentile(hc, 75):.4g}")
    print(f"   best Youden cut {best[1]:.4g}: sensitivity {best[2]:.0%}, specificity {best[3]:.0%}\n")


for label, sheet in [
    ("YeiH-specific, % of CD8 T cells            (Fig 4B left)", "Fig 4B left"),
    ("YeiH-specific, per million PBMC            (Fig 4B right)", "Fig 4B right"),
    ("Influenza NP-specific, % of CD8  [control] (Fig 4C left)", "Fig 4C left"),
    ("Influenza NP-specific, per million PBMC    (Fig 4C right)", "Fig 4C right"),
]:
    summarise(label, sheet)

c = cols("Supp Fig 3A left")
hc = c["Healthy Control"]
for arm in ("B27 AAU", "B27 axSpa"):
    v = c[arm]
    print(f"YeiH % of CD8, {arm} vs healthy: n={len(v)} vs {len(hc)}  "
          f"AUC = {auc(v, hc):.3f}  p = {mannwhitneyu(v, hc, alternative='greater').pvalue:.4f}")
