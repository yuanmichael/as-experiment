# Is there a reproducible blood transcriptome difference between AS patients and healthy controls?

A random-effects meta-analysis of **five independent studies, 97 AS patients vs 70 healthy controls**, built so that every case-vs-control comparison happens *inside* one study and only standardised effect sizes cross study boundaries.

**Short version:** yes, and it is concentrated in one place. The dominant, reproducible signal is a **reduced NK / cytotoxic lymphocyte programme in AS blood** (pooled *g* = −0.83 across nine canonical markers, present in all five studies, surviving removal of any one of them). A weaker but specific **type-1 interferon elevation** runs in the opposite direction to the genome-wide background and is unlikely to be an artefact. **AIF1 — the marker from the single-cell project — does not separate AS from healthy at all**, which is not a contradiction: it was a TNFi-response marker, not a disease marker.

![meta-analysis](figures/bulk_meta.png)

---

## Finding the data

GEO was searched programmatically (`scripts/geo_search.py`, NCBI E-utilities, five query spellings including "axial spondyloarthritis" and "spondyloarthropathy"), then **every candidate's sample labels were read directly from its GSM records** rather than trusted from the series summary (`scripts/geo_inventory.py`, 944 samples classified). ArrayExpress/BioStudies was searched separately and **added nothing** — every AS hit there is a GEO mirror.

That search produced two datasets the original shortlist missed, and killed one that looked like the biggest find in the field.

| Accession | Platform | Tissue | AS | HC | Verdict |
|---|---|---|---|---|---|
| **GSE73754** | Illumina HT-12 v4 | whole blood | 52 | 20 | primary |
| **GSE25101** | Illumina HT-12 v3 | whole blood | 16 | 16 | primary |
| **GSE18781** | Affy U133 Plus 2.0 | peripheral blood | 18 | 25 | primary — **not on the shortlist** |
| **GSE299639** | ONT long-read RNA-seq | peripheral blood | 6 | 6 | primary — **not on the shortlist** |
| **GSE181364** | Illumina RNA-seq | whole blood | 5 | 3 | primary — **not on the shortlist** |
| GSE221786 | Illumina RNA-seq | IL-17-enriched, CytoSim-stimulated | 20 | 8 | secondary |
| GSE11886 | Affy U133 Plus 2.0 | GM-CSF macrophages | 8 | 9 | secondary, flagged |
| GSE243689 | Illumina RNA-seq | sorted CD8 subsets | 6 | 5 | not used (per-sample files only) |
| GSE252708 | miRNA-seq | PBMC | 76 | 20 | other modality |
| GSE223717 | miRNA | sorted monocyte subsets | 128 | 63 | other modality |
| GSE179571 | methylation array | PBMC | 12 | 12 | other modality |
| GSE220130 | plasma exosomes | plasma | 17 | 10 | other modality |

**Dropped, with reasons:**

- **GSE39428** — 480 samples, titled *"Screening genes associated with rheumatoid arthritis and ankylosing spondylitis"*, 51 AS / 163 HC / 266 RA. It is a **384-SNP VeraCode genotyping array, not expression**. It would have been by far the largest cohort in the analysis and contributes nothing to it.
- **GSE141646** — 44 PBMC samples, 22 AS patients pre/post TNF inhibitor. **No healthy controls.** Useless for this contrast; it is the obvious validation cohort for the *previous* project's TNFi-response question.
- GSE266295, GSE232131, GSE194315, GSE288581 — AS-only, or single-cell.
- GSE122624, GSE304063, GSE30023 — SpA vs rheumatoid arthritis, no healthy arm.
- GSE58667 — juvenile spondyloarthritis; different disease.

**Two traps caught at the inventory stage:**

1. **GSE25101 annotates itself as both `tissue: Whole blood` and `cell type: PBMC`.** These are different tissues. GEO's own overall-design text settles it — *"RNA was extracted from whole blood using PAXGene tubes"* — so the `PBMC` field is a mis-annotation.
2. **GSE18781's SpA arm was run as two sets at different times, normalised independently**, and set 1 shares its controls with the study's sarcoidosis arm. The series matrix carries `set: 1` / `set: 2` as an explicit sample characteristic, so this is modelled as a covariate rather than assumed away.

---

## Method

Nothing is merged across studies. Per study, per gene:

```
log2(expression) ~ AS + sex + age + batch          (least squares, one design per study)

g      = J · beta_AS / sqrt(MSE)                   standardised mean difference, SD units
var(g) = (X'X)^-1[AS,AS] + g² / (2·dof)            design term + SD-estimate term
J      = 1 − 3/(4·dof − 1)                         Hedges' small-sample correction
```

All genes in a study share one design matrix, so a single least-squares solve handles every gene. Effect sizes are then pooled with **DerSimonian–Laird random effects**, reporting pooled *g*, 95% CI, Benjamini–Hochberg *q*, and between-study heterogeneity *I*².

Covariates are used where reported, but **disease-defining variables are deliberately excluded** — adjusting for HLA-B27 or disease activity would adjust away part of the contrast. A covariate collinear with disease status (|r| > 0.95) is dropped and the drop is printed.

Genes are tested if measured in **≥3 of the 5** studies: 11,269 genes.

### Power

| Study | AS | HC | smallest *d* at 80% power |
|---|---|---|---|
| GSE73754 | 52 | 20 | 0.74 |
| GSE25101 | 16 | 16 | 0.99 |
| GSE18781 | 18 | 25 | 0.87 |
| GSE299639 | 6 | 6 | 1.62 |
| GSE181364 | 5 | 3 | 2.05 |
| **pooled** | **97** | **70** | **0.44** (upper bound — random effects with τ² > 0 has less) |

Every individual study could only ever have found a large effect. Pooling is what buys the ability to see a moderate one.

---

## Two normalisation bugs, caught by the diagnostics and not by FDR

The first run reported 652 genes at *q* < 0.05, nearly all in one direction, with *I*² = 0 and 5/5 study agreement. That is what a real signal looks like, and it is also exactly what a normalisation artefact looks like. `scripts/bulk_diagnostics.py` separates them.

**Check 1 — global shift.** `GSE181364` had mean *g* = **−1.3 × 10¹²**. Two compounding errors, both mine: its supplementary table is **raw read counts**, which I had log2'd without depth-normalising (library size varied 2.5×, and tracked disease status — overall signal *d* = −1.92, *p* = 0.03), and at dof = 6 genes with near-zero residual variance send *g* = β/*s* to infinity. Fixed with CPM normalisation, an expression filter, and a residual-SD floor. `GSE221786` had the same count problem (4.9× library-size spread).

**Check 2 — GSE11886's matrix is already batch-corrected.** Its series matrix is per-gene median-centred log ratios **after Distance Weighted Discrimination**, applied by the submitters "to account for any batch effects". A supervised batch correction can delete or manufacture exactly the contrast being tested. It stays in the secondary analysis with its effect size flagged as untrustworthy; doing it properly needs the raw CELs.

**Check 3 — permutation null.** Shuffling disease labels within study and re-running the entire pipeline 20 times produced **0 hits every time** against 850 observed. This proves the pipeline does not manufacture signal from noise. It proves *less than it appears to*: shuffling also destroys any confound that tracks disease status, so it cannot distinguish biology from a disease-linked artefact. It is a necessary check, not a sufficient one.

**Check 4 — the one that carries the weight.** A global artefact moves every lineage. The observed pattern does not.

---

## Result

**850 genes at *q* < 0.05** (660 down in AS, 190 up). But the genome-wide background is *itself* mildly displaced — only **36.4%** of genes have pooled *g* > 0, median *g* = −0.135 — so a 50/50 null is the wrong test for any gene set. Everything below is calibrated against that 36.4% background.

### The signal is NK and T cell, and nothing else

| lineage | markers | mean pooled *g* | up / expected | Mann–Whitney vs all other genes |
|---|---|---|---|---|
| **NK** | 9 | **−0.834** | 0 / 3.3 | **5.9 × 10⁻⁷** |
| **T cell** | 8 | **−0.536** | 0 / 2.9 | **2.2 × 10⁻⁴** |
| B cell | 7 | −0.156 | 2 / 2.5 | 0.80 |
| monocyte | 8 | −0.051 | 3 / 2.9 | 0.50 |
| neutrophil | 8 | +0.104 | 4 / 2.9 | 0.038 |
| erythroid | 5 | +0.091 | 4 / 1.8 | 0.14 |

No global shift can push NK down 0.83 SD while leaving erythroid and monocyte flat. Individually: *EOMES* −1.03, *NCR3* −0.95, *KLRD1* −0.94, *NKG7*, *GNLY*, *PRF1* all at *q* < 0.05, with *I*² = 0.

**Leave-one-study-out** (panel D): NK stays between **−0.78 and −0.94** across all seven subsets, including arrays-only and dropping the dominant 52-patient study. T cell stays between −0.49 and −0.65. No single study carries it.

### Type-1 interferon: weak, but real and specific

18 of 22 interferon-alpha response genes are up in AS, against **8.0 expected** from the background. Sign test *p* < 10⁻⁴; Mann–Whitney against all other genes *p* = 6.3 × 10⁻⁵. The panel runs *against* the genome-wide downward drift, which is what makes it credible.

Magnitude is small (mean pooled *g* = +0.14) and only **MX2** survives gene-level FDR (*g* = +0.61, 95% CI 0.26–0.95, *q* = 0.014, *I*² = 0). This is a coherent set-level shift, not a set of individually significant genes — and it is the same direction as the interferon signal in the single-cell project.

### AIF1 does not separate AS from healthy

| | pooled *g* | 95% CI | *p* | *I*² |
|---|---|---|---|---|
| **AIF1** | +0.28 | −0.18 to +0.74 | 0.23 | 36% |

A clean null at *n* = 167, with 80% power to detect *d* = 0.44. **This is not a failed replication.** AIF1 was established as a predictor of *TNF-inhibitor response among AS patients* — a within-disease contrast. Nothing in that paper claims it distinguishes AS from healthy, and this result does not bear on it. The right test for AIF1 is GSE141646 (22 AS patients, pre/post TNFi, no controls needed).

### Monocyte markers

The single-cell project's CD14-monocyte finding was also a TNFi-response contrast and likewise is not tested here. For completeness: the monocyte lineage as a whole is flat (*g* = −0.05, *p* = 0.50), though *CSF1R* (−0.57, *q* = 0.023) and *CDKN1C* (−0.53, *q* = 0.038) are individually down.

### Secondary studies

Sorted/stimulated cells (GSE221786 + GSE11886, *k* = 2) yield **0 genes at *q* < 0.05** — expected at that *k* and those sample sizes. Of the 850 primary hits, 763 are measured there and **56% agree in direction** (sign test *p* = 5 × 10⁻⁴). Statistically significant but only 6 points above chance: **weak convergence**, and one of the two studies has a pre-applied supervised batch correction. It supports the primary result mildly; it does not independently confirm it.

---

## What this does and does not establish

**Does:** AS blood has a reproducibly reduced NK/cytotoxic and T-cell transcriptional programme, consistent across five cohorts, two platform families, three countries and fifteen years of array technology, robust to dropping any one study.

**Does not:** distinguish *composition* from *regulation*. A coordinated drop across an entire lineage's marker set is what you see when there are **fewer of those cells in the sample** — not necessarily when the cells present are behaving differently. Whole blood is a mixture, and this analysis measures the mixture. Neutrophil markers are only mildly up (+0.10), so it is not simple neutrophil-dilution of everything, but these data cannot resolve it further.

**That is exactly what Project 2a is for.** Deconvolving these same five cohorts with the annotated 314k-cell GSE277117 reference would turn "the NK programme is down" into "there are fewer NK cells" or "NK cells are transcriptionally altered" — and it would test the CD14-monocyte finding at *n* ≈ 167 instead of *n* = 25.

**Caveats on the record:**

- **GSE299639 and GSE181364 share an identical file template** (`Gene / samples / Symbol / Chr / GeneType / Description`), implying the same processing vendor and possibly overlapping cohorts. They are treated as independent; if they are not, the two smallest studies are double-counted. Together they carry ~13% of the pooled weight, and the leave-one-out columns show the result is unchanged without either.
- Median *I*² is 30%, and 34% of genes exceed *I*² = 50%. There is real between-study heterogeneity; the random-effects model accounts for it in the CI, but it means a single pooled number hides genuine disagreement for a third of the genome.
- No ML classifier was trained. Per the plan, that comes only after the signal supports it, with leave-one-study-out validation and feature selection inside the training fold.

---

## Reproducing

```bash
python scripts/geo_search.py        # GEO E-utilities sweep -> data/geo_candidates.tsv
python scripts/geo_inventory.py     # verify every GSM label -> data/geo_sample_inventory.tsv
python scripts/bulk_download.py     # series matrices + supplementary tables
python scripts/bulk_prep.py         # per-study normalise, annotate, collapse to genes
python scripts/bulk_meta.py         # effect sizes, DL random effects, panel, sensitivity
python scripts/bulk_diagnostics.py  # global shift, library size, permutation, lineage
python scripts/bulk_figure.py       # figures/bulk_meta.png
```

Runs locally in a few minutes; no Modal needed — the bulk matrices are small.
