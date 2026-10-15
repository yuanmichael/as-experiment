# Independent reanalysis of GSE277117 [Song et al., *Nature Communications* 16:5531 (2025)](https://doi.org/10.1038/s41467-025-60445-6)

Amongst Ankylosing Spondylitis patients, are there discernable differences in their single-cell sequencing data between responders to TNF inhibitors (like Humira, etc.) vs. non-responders? In theory, if there is a meaningful signal, you could perform a blood draw before prescribing AS patients to an anti-TNF inhibitor. Because ~40% of patients don't respond well to such medication, and may require an IL-17 inhibitor 

A from-scratch reanalysis of the data, built independently of the authors' pipeline.

**Short version:** two of the paper's single-cell claims replicate cleanly. One is directionally consistent but underpowered. One could not be tested reliably. One is not testable from this dataset at all. Along the way, three technical traps each produced a confident wrong answer — those are documented here for learnings' sake. 

---

## The data

| | |
|---|---|
| Accession | [GSE277117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277117) (`GSE277117_All_cell.h5ad`, 5.03 GB) |
| Contents | 357,651 cells × 36,764 features, **raw integer counts** |
| Processing already done | Patient demultiplexing (souporcell) and multiplet flagging only |
| Processing *not* done | QC, normalisation, clustering, cell type annotation |
| Patients | 25 (KAS01–KAS25), PBMCs, baseline and 1 month |

### The cohort is the binding constraint

| Cohort | Definition | Responders | Non-responders |
|---|---|---|---|
| **A (strict)** | TNF-arm baselines only | 12 | **4** |
| **B (expanded)** | all baselines; IL-17-arm patients counted as TNFi non-responders | 12 | **13** |

Every IL-17-arm patient is a TNFi non-responder by definition — failing TNFi is the admission criterion for that arm. That makes `Drug` partly confounded with the outcome, so cohort B buys sample size at the cost of a batch/timing confound. Both were analysed throughout.

**Power.** At 80% power, α = 0.05, two-sided:

```
cohort B (12 vs 13):  detectable effect size  d = 1.17
cohort A (12 vs 4):   detectable effect size  d = 1.74
power to detect a moderate effect (d = 0.5):  22%  /  13%
```

A Cohen's *d* of 1.17 means two groups whose distributions barely overlap. That is the only size of effect this dataset could have found.

---

## Pipeline

```
raw h5ad (357,651 × 36,764)
  → drop 163 Antibody Capture features        → 36,601 genes
  → drop Multiplet cells                      → 326,327
  → QC: 200–6000 genes/cell, ≤15% mito        → 314,368
  → genes in ≥3 cells                         →  31,813
  → CPM(1e4) + log1p, 2000 HVGs (batch-aware)
  → PCA(50) → Harmony(Sample_num) → Leiden
  → CellTypist (Immune_All_Low), majority voting
  → 23 cell types; pseudobulk per patient × cell type
```

**Clustering resolution** was chosen empirically rather than by convention. Leiden at 0.2–3.0 scored against independent per-cell reference labels:

| resolution | clusters | ARI | NMI |
|---|---|---|---|
| 0.2 | 12 | 0.552 | 0.619 |
| **0.6** | **20** | **0.571** | **0.627** |
| 1.0 | 27 | 0.499 | 0.607 |
| 3.0 | 50 | 0.231 | 0.552 |

---

## Scorecard

| Paper's claim | Verdict | Evidence |
|---|---|---|
| **Fig 1c** — CD14 monocytes ↓, T cells ↑ in non-responders at baseline | ✅ **Replicates** | CD14 33.6% → 23.0%, p = 0.011; T cell 33.8% → 40.9%. Near-identical to their reported percentages, from an independent pipeline |
| **Fig 3b** — AIF1 marks the ISG-high monocyte state | ✅ **Replicates** | Spearman r = 0.221, p = 7e-141 in non-classical monocytes (n = 91,275 cells) |
| **Fig 2e** — ISG-high ncMo fraction ↑ in non-responders | 🟨 **Directionally consistent, not significant** | 12/12 threshold comparisons higher in NR; ratio rises monotonically 1.3× → 2.2× as threshold tightens; all p > 0.34 |
| **Fig 7b** — activated Th17 / Treg ↑ in non-responders | 🟨 **Not replicated** | Right direction (19.3% vs 16.0%, p = 0.18), but the method needed to test it is unreliable at this n — see trap 3 |
| **Fig 3d–h** — serum AIF-1 predicts response, AUC 0.822 | ⬜ **Not testable here** | Established in separate bulk-RNA (n = 18) and serum-ELISA (n = 117) cohorts |

**The pattern:** the two claims that replicate are the two that reduce to cell-level or large-effect measurements. Everything that reduces to comparing 25 people is underpowered.

### Detail: the interferon signal

Per-cell Hallmark `INTERFERON_ALPHA_RESPONSE` score, fraction of cells above threshold, cohort B:

| threshold | non-responders | responders | ratio | p |
|---|---|---|---|---|
| top 25% | 28.4% | 21.6% | 1.31× | 0.53 |
| top 10% | 12.2% | 7.8% | 1.56× | 0.53 |
| top 5% | 6.5% | 3.5% | 1.83× | 0.34 |
| top 1% | 1.30% | 0.61% | 2.13× | 0.37 |

Consistent direction and a dose-response with stringency — but driven by a handful of non-responders with very high ISG fractions, while medians barely differ. **The interferon signal looks like a property of a subgroup of non-responders, not of non-responders as a class.** At n = 13, a subgroup effect is undetectable.

The per-cell score distribution is **unimodal with a long tail** — there is no discrete ISG-high population. Any "ISG-high" cluster is a cut through a gradient, and its size depends on where you cut.

---

## Traps encountered

Each of these produced a confident, plausible, wrong answer before being caught.

### 1. Pseudo-replication

357,651 cells, 25 patients. The outcome label is a patient attribute stamped onto every one of their cells. Splitting cells randomly into train/test puts the same patient on both sides, and the model learns to recognise the patient, not the biology. **Effective n is 25.**

### 2. A hidden confounded feature

The matrix mixes 36,601 genes with **163 CITE-seq antibody-capture features**, and assay type tracks the outcome almost perfectly:

```
CITE-seq signal present:   29/29 non-responder samples  (100%)
                           10/24 responder samples      ( 42%)
                           21/21 IL-17 arm              (100%)
```

Antibody tags were ~7% of total counts in non-responder groups and ~1% in responders. Including them in the CPM denominator deflated **every gene** in non-responders by ~6%, producing a clean, consistent, entirely artefactual "lower in non-responders" result across 9 of 10 cell types — in the exact direction needed to look like a failed replication.

FDR correction did **not** catch this. Multiple-testing control filters noise; it cannot detect systematic bias. Only reading the feature list does.

Filter with `adata.var["feature_types"] == "Gene Expression"` **by index, not name** — several antibodies share a symbol with a gene (`CD14`, `CD19`, `CD4`, `CD74`), so `var_names_make_unique()` renames one copy to e.g. `CD44-1` and name matching is unreliable.

### 3. Subclustering within T-cells / monocytes

The paper's findings are about **proportions of sub-populations identified by re-clustering within a compartment** (monocytes alone; memory CD4 T cells alone). We initially averaged expression across whole cell types. Averaging a signature over all monocytes cannot see a 1–3% subpopulation. The first round of analysis tested the wrong quantity and reported the null as informative.

---

## Figures

| File | Contents |
|---|---|
| `figures/qc_overview.png` | Genes/cell, % mitochondrial, depth-vs-complexity; the dying-cell mode at 75–95% mito |
| `figures/umap_batch_check.png` | UMAP before/after Harmony, coloured by cluster, batch and patient |
| `figures/annotation_v2.png` | 23 CellTypist populations; IL-17-axis cells highlighted |
| `figures/composition_baseline.png` | **Main result** — baseline composition by response, both cohorts, laid out like the paper's Fig 1c |
| `figures/isg_direct.png` | Per-cell IFN score distribution (unimodal) and ISG-high fraction by response |
| `figures/cd4.png` | Memory CD4 sub-clusters without batch correction, and the patient-specific artefacts it produces |

---

## Reproducing

```bash
python -m venv .venv && .venv/Scripts/pip install -r requirements.txt
modal run scripts/modal_app.py::download      # 5 GB GEO -> Modal volume
modal run --detach scripts/process.py         # QC, normalise, cluster
modal run --detach scripts/annotate_clean.py  # CellTypist
modal run --detach scripts/pseudobulk.py      # patient-level table
modal run --detach scripts/stats.py           # tests + FDR
modal run --detach scripts/isg_direct.py      # threshold test, no clustering
modal run --detach scripts/cd4_nohrm.py       # CD4 compartment
```

| Script | Purpose |
|---|---|
| `modal_app.py` | Modal app, image, volume, data download |
| `qc.py`, `qc_followup.py` | QC metrics; whether quality tracks the outcome |
| `process.py` | Antibody filter, QC, normalise, HVG, PCA, Harmony, Leiden |
| `annotate_clean.py` | CellTypist reference mapping |
| `pseudobulk.py` | Composition, pseudobulk, signature scores |
| `stats.py` | Patient-level tests, FDR, genome-wide scan |
| `subcluster.py`, `cd4_nohrm.py`, `isg_direct.py` | Compartment-level analyses |
| `check_features.py`, `check_hvg.py` | Antibody contamination forensics |

---

## What this dataset cannot do

A genome-wide scan at baseline — **158,440 gene × cell-type tests** — produced **zero** results surviving FDR at q < 0.05. That is the honest characterisation of what 25 patients support, and it is not a criticism of the study: the paper's predictive claim rests on 225 patients across four cohorts, with the single-cell arm serving to generate the hypothesis rather than test it.

---

*Data: Song W, Lee EE, Park S, et al. Nat Commun 16:5531 (2025), GEO GSE277117.*
