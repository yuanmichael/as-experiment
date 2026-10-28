"""Pool the motif-family result across cohorts and derive diagnostic metrics."""
import numpy as np
from scipy.stats import fisher_exact, beta, hypergeom

# (cohort, compartment, AS carriers, n AS, ctrl carriers, n ctrl)
HEALTHY = [
    ("Komech 2018", "unsorted blood, RNA", 14, 25, 0, 10),
    ("Hanson 2020", "sorted CD8, gDNA", 22, 48, 2, 36),
    ("GSE288581", "sorted CD8 mem, 10x", 3, 4, 0, 4),
]
PSA = [
    ("PRJEB51315 blood", 5, 7, 5, 16),
    ("PRJEB51315 synovial", 5, 7, 6, 17),
]
KNOWN = [
    ("Komech 2018", 11, 25, 0, 10),
    ("Hanson 2020", 22 - 10, 48, 2, 36),  # placeholder, overwritten below
]


def ci(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


print("AS vs HEALTHY — union of the 7 pre-specified sequences\n")
print(f"{'cohort':<16}{'compartment':<22}{'AS':>10}{'sens':>8}{'ctrl':>10}{'spec':>8}{'p':>11}")
num_or, den_or, sens_k, sens_n, spec_k, spec_n = 0.0, 0.0, 0, 0, 0, 0
for name, comp, a, na, b, nb in HEALTHY:
    p = fisher_exact([[a, na - a], [b, nb - b]], alternative="greater")[1]
    print(f"{name:<16}{comp:<22}{f'{a}/{na}':>10}{a/na:>8.1%}{f'{b}/{nb}':>10}"
          f"{(nb-b)/nb:>8.1%}{p:>11.2g}")
    n = na + nb
    num_or += a * (nb - b) / n
    den_or += b * (na - a) / n
    sens_k += a; sens_n += na; spec_k += (nb - b); spec_n += nb

print(f"\nMantel-Haenszel pooled odds ratio: {num_or/den_or:.1f}")
lo, hi = ci(sens_k, sens_n)
print(f"pooled sensitivity {sens_k}/{sens_n} = {sens_k/sens_n:.1%}  (95% CI {lo:.1%}-{hi:.1%})")
lo2, hi2 = ci(spec_k, spec_n)
print(f"pooled specificity {spec_k}/{spec_n} = {spec_k/spec_n:.1%}  (95% CI {lo2:.1%}-{hi2:.1%})")

se, sp = sens_k / sens_n, spec_k / spec_n
print("\nPPV / NPV at four prevalences (healthy comparator — optimistic)")
print(f"{'prevalence':>12}{'PPV':>10}{'NPV':>10}")
for pv in (0.005, 0.05, 0.20, 0.50):
    ppv = se * pv / (se * pv + (1 - sp) * (1 - pv))
    npv = sp * (1 - pv) / (sp * (1 - pv) + (1 - se) * pv)
    print(f"{pv:>12.1%}{ppv:>10.1%}{npv:>10.1%}")

print("\n\nAS vs PSORIATIC ARTHRITIS — the clinical differential")
print(f"{'cohort':<24}{'AS':>10}{'PsA':>10}{'spec vs PsA':>14}{'p':>10}")
for name, a, na, b, nb in PSA:
    p = fisher_exact([[a, na - a], [b, nb - b]], alternative="greater")[1]
    print(f"{name:<24}{f'{a}/{na}':>10}{f'{b}/{nb}':>10}{(nb-b)/nb:>14.1%}{p:>10.2g}")
print("  (single clone CASSVGLYSTDTQYF, blood: 5/7 AS vs 4/16 PsA, "
      f"p = {fisher_exact([[5,2],[4,12]], alternative='greater')[1]:.2g})")

print("\n\nWHAT THE CONTROL ARM BUYS — smallest reachable p for a marker in")
print("k cases and 0 controls, by control-arm size")
print(f"{'n controls':>12}{'  cases needed for p<=0.01':>30}")
for n_case, n_ctrl, label in ((25, 10, "Komech"), (48, 36, "Hanson"),
                              (73, 46, "pooled")):
    need = None
    for k in range(1, n_case + 1):
        p = hypergeom.sf(k - 1, n_case + n_ctrl, k, n_case)
        if p <= 0.01:
            need = k
            break
    print(f"{n_ctrl:>12}{f'{need} of {n_case} ({need/n_case:.0%})':>30}   {label}")
