# Reduced NK-cell cytotoxic transcription in ankylosing spondylitis is not explained by NK-cell abundance

**Draft preprint — for review. Not submitted.**

Target: bioRxiv (Immunology). ~3,200 words.

**Author:** Michael Yuan
**Affiliation:** *[to complete]*
**Correspondence:** *[to complete]*

---

## Abstract

Blood transcriptomic studies of ankylosing spondylitis (AS) consistently report reduced expression of natural killer (NK) cell genes. Whether this reflects fewer circulating NK cells or reduced transcription per cell has not been established, and the two have different implications. A meta-analysis of 95 flow-cytometry studies (4,020 patients, 3,065 controls) reports no difference in NK cell percentage, which is difficult to reconcile with the transcriptomic literature.

We assembled five independent AS blood cohorts (97 patients, 70 healthy controls; three microarray, two RNA-seq) and computed AS-versus-control effect sizes within each study, combining them by random-effects meta-analysis so that no comparison crossed a study boundary. We then deconvolved all samples against a 312,325-cell single-cell PBMC reference to estimate cell-type proportions, and separated proportion from per-cell transcription using disjoint gene sets: NK proportion was estimated from identity genes with all effector genes removed from the signature matrix, and the NK effector programme was tested as the outcome.

NK marker genes showed the strongest lineage-level reduction in AS (mean pooled *g* = −0.83; Mann–Whitney versus all other genes *P* = 5.9 × 10⁻⁷), robust to leave-one-study-out (range −0.78 to −0.94). Estimated NK proportion was modestly reduced (*g* = −0.54, 95% CI −0.91 to −0.17, *P* = 0.004). The NK effector programme was reduced by *g* = −0.98 (*P* = 4 × 10⁻⁸), and **82% of this effect persisted after adjustment for NK proportion** (*g* = −0.80, *P* = 1 × 10⁻⁵; *I*² = 0 for both).

Reduced NK cytotoxic transcription in AS blood is therefore predominantly a per-cell phenomenon rather than a change in cell abundance. Because transcript-based deconvolution mistakes a transcriptionally quiescent cell for an absent one, our proportion estimate is an upper bound; combined with the flow-cytometry evidence, the composition component is likely near zero. The effect is not diagnostically useful (AUC 0.77 against healthy controls; positive predictive value 17.7% at back-pain-clinic prevalence, versus 37.2% for HLA-B27).

---

## Introduction

*[Draft — needs citation numbering and a fuller treatment of AS immunopathology.]*

Ankylosing spondylitis is a chronic immune-mediated arthritis in which the IL-23/IL-17 axis and HLA-B27 are central, but the contribution of innate cytotoxic lymphocytes is less settled. Several blood transcriptomic studies have reported reduced expression of NK-associated genes in AS, and functional genomics integrating AS GWAS signal with leukocyte expression has implicated NK cells in pathogenesis. Single-cell work has described reduced NK frequency with impaired cytotoxic gene expression and lower plasma granzymes.

Against this, a meta-analysis of peripheral lymphocyte subsets pooling 95 flow-cytometry studies found **no significant difference in NK cell percentage** between AS patients and controls. Flow cytometry enumerates cells by surface protein and is insensitive to transcriptional state; bulk transcriptomics measures an average over a cell mixture and cannot, on its own, distinguish a change in mixture from a change in per-cell output. The two literatures are therefore not necessarily in conflict — but no study has tested that directly.

Two further features of the existing literature motivated a reanalysis. First, published integrative analyses of the two largest AS blood datasets (GSE25101 and GSE73754) merged them into a single expression matrix before differential testing, which permits comparison of one study's cases to another study's controls and confounds disease with batch. Second, no analysis has attempted to decompose the NK signal into abundance and per-cell components.

We therefore (i) re-derived the AS blood transcriptome signal using a per-study random-effects design across five cohorts, and (ii) used single-cell-reference deconvolution with a disjoint-gene-set design to separate NK abundance from NK per-cell transcription.

---

## Methods

### Dataset identification

GEO was queried programmatically (NCBI E-utilities, `esearch`/`esummary` against the `gds` database) using five spellings: "ankylosing spondylitis", "axial spondyloarthritis", "spondyloarthritis", "spondyloarthropathy", "radiographic axial spondyloarthritis". Results were restricted to human expression or methylation profiling series. ArrayExpress/BioStudies was searched separately and returned no series not already mirrored in GEO.

Sample labels were verified directly from GSM records (944 samples classified) rather than from series summaries. This was consequential: GSE39428, which is titled as a study of genes associated with rheumatoid arthritis and ankylosing spondylitis and contains 51 AS, 266 RA and 163 control samples, is a 384-SNP VeraCode genotyping array, not expression data. GSE25101 carries contradictory `tissue: Whole blood` and `cell type: PBMC` annotations; its series design text specifies PAXgene whole blood, and neutrophil transcript abundance in the data (85th expression percentile) confirms whole blood. GSE181364, annotated as whole blood, has neutrophil markers at the 39th percentile and is neutrophil-depleted.

### Cohorts

Five cohorts of circulating, unstimulated blood with both AS cases and healthy controls were included.

| Accession | Platform | Tissue | AS | HC | Covariates modelled |
|---|---|---|---|---|---|
| GSE73754 | Illumina HT-12 v4 | whole blood | 52 | 20 | sex |
| GSE25101 | Illumina HT-12 v3 | whole blood | 16 | 16 | — (matched by design) |
| GSE18781 | Affymetrix U133 Plus 2.0 | peripheral blood | 18 | 25 | sex, age, normalisation set |
| GSE299639 | ONT long-read RNA-seq | peripheral blood | 6 | 6 | sex |
| GSE181364 | Illumina RNA-seq | whole blood (neutrophil-depleted) | 5 | 3 | — |

Studies of sorted or stimulated cells (GSE221786, GSE11886), other modalities (miRNA, methylation, exosome proteomics) and AS-only cohorts without healthy controls (e.g. GSE141646) were excluded from the primary analysis. GSE18781's spondyloarthropathy arm was processed as two independently normalised sets; set membership is recorded in the series matrix and was modelled as a covariate.

### Preprocessing

Microarray series were taken from GEO series matrices and probes mapped to gene symbols using the platform annotation, retaining per gene the probe with highest mean signal. RNA-seq count tables were converted to counts per million before log transformation; this step is not optional, as library size varied 2.5-fold in GSE181364 and correlated with disease status, which would otherwise impose a genome-wide shift on every gene. Low-expression features were removed (microarray: bottom quintile by mean; RNA-seq: expressed above 1 CPM in <50% of samples).

### Per-study effect sizes and meta-analysis

For each study and gene, log₂ expression was modelled as a function of disease status and available covariates by least squares, with all genes in a study sharing one design matrix. Disease-defining variables (HLA-B27 status, disease activity) were excluded as covariates. The disease coefficient was converted to a standardised mean difference:

*g* = *J* · β_AS / √MSE, with Var(*g*) = (X′X)⁻¹[AS,AS] + *g*²/(2·df) and *J* = 1 − 3/(4·df − 1).

Genes measured in ≥3 of 5 studies (*n* = 11,269) were combined by DerSimonian–Laird random effects, reporting pooled *g*, 95% CI, Benjamini–Hochberg *q* and *I*².

### Deconvolution

A cell-type signature matrix was built from the annotated GSE277117 PBMC object (312,325 cells, 25 AS patients, 23 CellTypist types, including 20,105 CD16⁺ NK cells). Non-gene features were removed by index on `feature_types`. Profiles were computed from raw counts, CPM-normalised per cell. Markers were selected per cell type by the largest log₂ margin between that type's mean and the next-highest type's (120 per type; 2,760 unique). Granulocyte and globin genes were excluded, as the reference is PBMC-derived while two cohorts are whole blood.

Patients were split in half: the signature was built from 12 and held-out pseudobulk (summed counts per patient, with known proportions) from the remaining 13. Mixtures were estimated by ν-support-vector regression on standardised matrices, with negative coefficients clipped and the remainder normalised to sum to one.

### Separating abundance from per-cell transcription

Estimating a cell-type proportion from expression data and then adjusting an expression signature for that proportion is circular, because both quantities derive from the same transcripts. We therefore partitioned NK genes into:

- **identity** (KLRD1, KLRF1, NCR1, NCR3, EOMES, TBX21, IL2RB, CD247, S1PR5, KLRC1, PLEKHF1, SH2D1B, PRSS23, CD160, KIR2DL3, FCGR3A, TXK, CLIC3)
- **effector** (GZMA, GZMB, GZMH, GZMK, GZMM, PRF1, GNLY, NKG7, CTSW, SERPINB9, FASLG)

All effector genes were deleted from the signature matrix before proportion estimation. The NK effector score (mean within-study *z* across effector genes) was then modelled with and without the identity-derived NK proportion as a covariate.

All code is in the repository below; analyses are reproducible from GEO accessions without manual steps.

---

## Results

### A reproducible NK-dominant signal

850 genes reached *q* < 0.05 (660 reduced, 190 increased in AS). The genome-wide distribution is itself mildly displaced (36.4% of genes with *g* > 0; median *g* = −0.135), so gene sets were calibrated against that background rather than against a 50/50 null.

| Lineage | Markers | Mean pooled *g* | Up / expected | Mann–Whitney vs other genes |
|---|---|---|---|---|
| **NK** | 9 | **−0.834** | 0 / 3.3 | **5.9 × 10⁻⁷** |
| **T cell** | 8 | **−0.536** | 0 / 2.9 | **2.2 × 10⁻⁴** |
| B cell | 7 | −0.156 | 2 / 2.5 | 0.80 |
| Monocyte | 8 | −0.051 | 3 / 2.9 | 0.50 |
| Neutrophil | 8 | +0.104 | 4 / 2.9 | 0.038 |
| Erythroid | 5 | +0.091 | 4 / 1.8 | 0.14 |

The NK effect was present in all five cohorts (per-study *g* for EOMES: −1.08, −0.69, −1.13, −1.28, −1.41; *I*² = 0) and remained between −0.78 and −0.94 across all leave-one-study-out subsets, including microarray-only and excluding the largest cohort. Permuting disease labels within study and re-running the full pipeline produced 0 genes at *q* < 0.05 in 20 of 20 permutations.

A subset of NK genes showed substantial heterogeneity (KLRF1, KLRB1, KLRC1, GZMA, GZMK; *I*² 57–75%), driven principally by GSE25101. These are predominantly killer-lectin receptors and granzymes; identity genes were homogeneous across studies. We note this but do not interpret it.

### Deconvolution validation

Recovery of known proportions from held-out pseudobulk:

| Lineage | True mean | *r* with truth | Interpreted |
|---|---|---|---|
| B | 16.0% | 0.944 | yes |
| Monocyte | 28.7% | 0.914 | yes |
| CD4 T | 14.9% | 0.899 | yes |
| NK | 15.6% | 0.833 | yes |
| DC | 1.4% | 0.759 | marginal |
| CD8 T | 22.6% | 0.567 | **no** |
| Platelet | 0.6% | 0.015 | **no** |

CD8 T cells were not recoverable, consistent with Tem/Temra cytotoxic T cells sharing GNLY, NKG7 and PRF1 with NK cells. CD8 T and platelet estimates are reported for completeness but not interpreted, including where nominally significant (CD8 T proportion *g* = −0.376, *P* = 0.029).

### NK abundance is modestly reduced

| Lineage | Pooled *g* | 95% CI | *P* | *I*² |
|---|---|---|---|---|
| **NK** | **−0.540** | −0.91 to −0.17 | **0.004** | 11% |
| Monocyte | +0.309 | −0.04 to +0.65 | 0.080 | 3% |
| CD4 T | +0.043 | −0.35 to +0.44 | 0.83 | 20% |
| B | −0.043 | −0.59 to +0.50 | 0.88 | 53% |

In real bulk samples 30–59% of signal was absorbed by rare reference types acting as a sink (versus 1.8% in pseudobulk), reflecting the unmodelled granulocyte and reticulocyte compartment. The sink was slightly larger in AS (*g* = +0.30 to +0.52 in four of five cohorts), which could inflate apparent decreases through renormalisation. Renormalising to the identified compartment moved the NK estimate only from −0.540 to −0.486 (*P* = 0.005, *I*² = 0).

### Most of the cytotoxic deficit is per-cell

| | Pooled *g* | 95% CI | *P* | *I*² |
|---|---|---|---|---|
| NK effector programme, unadjusted | −0.978 | −1.33 to −0.63 | 4 × 10⁻⁸ | 0% |
| **NK effector programme, adjusted for NK proportion** | **−0.804** | −1.17 to −0.44 | **1 × 10⁻⁵** | 0% |
| NK proportion (identity genes only) | −0.540 | −0.91 to −0.17 | 0.004 | 11% |

82% of the effector effect persists after adjustment. Individual effector genes: NKG7 −0.87 (*q* < 0.001), GZMK −0.94 (*q* = 0.038), GZMA −0.88, CTSW −0.77, GNLY −0.74, PRF1 −0.73 (*q* = 0.003), GZMB −0.61 (*q* = 0.013), GZMH −0.46; 0 of 8 increased against 2.9 expected (Mann–Whitney *P* = 1.0 × 10⁻⁸).

### Secondary observations

**Type-1 interferon.** 18 of 22 interferon-α response genes increased against 8.0 expected (binomial *P* < 10⁻⁴; Mann–Whitney *P* = 6.3 × 10⁻⁵), running against the genome-wide downward drift. The magnitude is small (mean *g* = +0.14) and only MX2 survives gene-level correction (*g* = +0.61, *q* = 0.014). Notably, the interferon score achieves AUC 0.496 — a statistically detectable mean shift with no classification ability.

**Platelet / integrin.** 10 of 11 canonical platelet-integrin genes increased against 4.0 expected (binomial *P* = 0.0003; Mann–Whitney *P* = 5.4 × 10⁻⁶; mean *g* = +0.44), with ITGA2B *g* = +0.83 (*q* = 0.001, *I*² = 0), TREML1 +0.65, THBS1 +0.64, ITGB3 +0.61, MYL9 +0.59. A recent deep plasma proteomic study of 104 samples independently identified integrin-pathway proteins (FERMT3, ILK, TLN1, with SAA1) as its leading AS biomarkers. The convergence is at pathway rather than gene level — those three specific transcripts are not individually significant here (mean *g* = +0.20, *P* = 0.11). We did not pursue this further and flag it as a hypothesis.

**AIF1.** Pooled *g* = +0.28 (95% CI −0.18 to +0.74, *P* = 0.23, *I*² = 36%), at 80% power for *g* = 0.44. AIF1 was established as a predictor of TNF-inhibitor response within AS patients, a contrast these data cannot address; this null does not bear on that claim.

### Diagnostic performance

| | AUC | Sensitivity | Specificity | PPV at 5% prevalence |
|---|---|---|---|---|
| NK effector score | 0.77 | 0.67 | 0.84 | 17.7% |
| CRP (published) | 0.71 | 0.35 | 0.90 | 15.6% |
| HLA-B27 (published) | 0.91 | 0.90 | 0.92 | 37.2% |

At the prevalence of a chronic back pain clinic, a positive NK score corresponds to 4.6 false positives per true case. Negative predictive value is 98.0%. These figures are obtained against healthy volunteers and would be expected to degrade substantially against the clinically relevant comparator of mechanical back pain.

---

## Discussion

The reduced NK transcriptional signature in AS blood is predominantly a reduction in per-cell cytotoxic transcription rather than a reduction in NK cell number. This reconciles the transcriptomic and flow-cytometric literatures: flow cytometry enumerates cells by surface protein and is blind to transcriptional state, so it correctly reports unchanged NK percentages while transcriptomics correctly reports a reduced NK programme.

Our deconvolution-derived abundance estimate (*g* = −0.54) should be read as an upper bound on the composition component. Transcript-based deconvolution cannot distinguish an absent cell from a transcriptionally quiescent one; if AS NK cells downregulate identity genes alongside effector genes, abundance is systematically underestimated. The flow-cytometry meta-analysis, which is not subject to this failure mode and is far larger, suggests the true abundance shift is smaller still. Both lines of evidence favour regulation.

This is consistent with single-cell work reporting impaired cytotoxic gene expression and reduced plasma granzymes in AS, with granzyme B inversely correlated with disease activity.

Mechanistically, reduced NK cytotoxic output is compatible with impaired clearance of infected or stressed cells, and with an NK regulatory role in restraining T-cell responses; we do not adjudicate between these.

### Limitations

**The signature reference derives from AS patients only.** GSE277117 contains no healthy controls, so cell-type profiles are AS-conditioned. If AS alters the NK expression profile, the signature is mis-specified for control samples. No public AS-and-control PBMC reference of comparable size exists.

**All controls are healthy volunteers.** The clinically relevant comparator is mechanical back pain. No public dataset addresses this, and it limits any diagnostic interpretation.

**No external validation of the per-cell effect.** All estimates derive from the same five cohorts. Targeted plasma granzyme/perforin measurement in an independent cohort would be the direct test. Publicly deposited AS plasma proteomics (ProteomeXchange/iProX IPX0011977000) may permit this without new sample collection.

**Marker coverage varies three-fold across cohorts** (469–1,378 of 2,739 signature genes), so deconvolution precision is uneven, worst in GSE299639.

**Tissue.** This analysis is confined to blood. AS pathology is in the entheses, sacroiliac joints and, subclinically, the gut; circulating cells may be a poor proxy.

**Replication, not discovery.** The direction of the NK finding has been reported previously. Our contribution is the per-study meta-analytic design and the abundance/per-cell decomposition, not the observation that NK genes are reduced.

---

## Data and code availability

All data are public. Accessions: GSE73754, GSE25101, GSE18781, GSE299639, GSE181364, GSE277117.

Code: https://github.com/yuanmichael/as-experiment — every analysis re-fetches its inputs from GEO; no manual preprocessing steps.

## Author contributions

*[to complete]*

## Competing interests

*[to complete — note patient status if you wish to declare it]*

## Acknowledgements

Analysis was carried out with assistance from Claude (Anthropic). *[Check target journal/preprint-server policy on AI-assistance disclosure; bioRxiv currently expects disclosure in Methods or Acknowledgements.]*

---

## Open items before submission

1. **External validation** of the per-cell effect — the single most important gap. Try IPX0011977000 for granzyme/perforin peptides first.
2. Citation numbering and a complete reference list; the Introduction needs a fuller treatment of AS immunopathology.
3. Decide whether to report the platelet/integrin observation as a named finding or leave it in secondary observations. It is currently the second-strongest signal in the data and under-explored.
4. Figures: currently two composite panels. A journal version would want the forest plot, the validation scatter and the decomposition as separate numbered figures with full legends.
5. Supplementary tables: per-study effect sizes (11,269 genes × 5 studies), validation metrics, sample inventory.
6. Confirm whether an institutional affiliation is needed — bioRxiv accepts unaffiliated submissions but requires a named corresponding author.
