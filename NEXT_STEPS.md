# Handoff prompt — next phase of the AS transcriptomics project

*Paste this into a fresh Claude Code session started in `C:\Coding\as-experiment`.*

---

I'm a patient with ankylosing spondylitis learning single-cell and bulk transcriptomics. I just finished an independent reanalysis of GSE277117 (Song et al., Nat Commun 16:5531, 2025 — type-1 IFN and AIF-1 in TNF-inhibitor response in AS). The full write-up is in `REPORT.md` in this repo — **read it first**, along with `scripts/` and the memory files, before proposing anything.

**How I'd like you to work with me:** act as a tutor. Do the low-level grunt work (data cleaning, pipeline code, file prep) yourself, but surface the high-level decisions and check my understanding with questions rather than handing me conclusions. Explain concepts from fundamentals — don't assume I know. Grade my answers honestly, including when I'm wrong.

## What already exists

- Modal app `as-experiment`, volume `as-data`, containing a cleaned, QC'd, annotated AnnData object: `/data/annotated_clean/as_annotated.h5ad` — **314,368 cells × 31,813 genes, 25 AS patients, 23 CellTypist cell types, raw counts preserved in `layers["counts"]`**.
- Local venv at `.venv`, `requirements.txt`, scripts in `scripts/`, figures in `figures/`.
- Key established findings: CD14 monocytes 33.6% → 23.0% in TNFi non-responders (p = 0.011, replicates the paper's Fig 1c); AIF1 expression correlates with the interferon program within monocytes (Spearman r = 0.22, p < 1e-140).
- **GSE277117 contains no healthy controls.** All 25 patients have AS. It cannot be used for any AS-vs-healthy contrast directly.

## Project 1 (do this first): bulk RNA-seq meta-analysis, AS vs healthy controls

**Goal:** is there a reproducible gene expression difference between AS patients and healthy controls in blood, across independent studies? Then, if the signal supports it, test whether a classifier generalises.

**Design — this is not a merge-everything-and-train exercise.** Compute the AS-vs-control effect size *within each study separately*, then combine across studies with a random-effects meta-analysis. The batch effect cancels because every comparison happens inside one study. Never compare study A's patients to study B's controls.

Candidate datasets (verify sample labels directly — GEO summaries have been unreliable):

| Accession | Platform | n | Tissue | Notes |
|---|---|---|---|---|
| GSE73754 | array | 72 | whole blood | largest; also splits by sex |
| GSE25101 | array | 32 | whole blood | 16 AS vs 16 control |
| GSE221786 | RNA-seq | 28 | IL-17-producing cells | AS vs control, sex-split |
| GSE11886 | array | 33 | monocyte-derived macrophages | 8 AS vs 9 HC |
| GSE179571 | methylation array | 24 | PBMC | 12 v 12; different modality |

Realistic total: ~90 AS, ~60 controls across four expression studies. Search GEO for more — this list is not exhaustive and newer RNA-seq cohorts may exist.

**Steps:**
1. Inventory: for each candidate, confirm it has both AS cases and healthy controls, confirm the tissue, and confirm cases and controls were processed in the same batches. Drop any study where they weren't.
2. Per study: standard preprocessing for its platform. Do not merge raw data across platforms.
3. Per study: differential expression, AS vs control, with covariates for sex and age where available.
4. Meta-analyse effect sizes (random-effects; `metafor` in R or `PyMARE`/manual DerSimonian-Laird in Python). Report pooled effect, CI, and between-study heterogeneity (I²).
5. Check whether the genes we care about move: `AIF1`, Hallmark `INTERFERON_ALPHA_RESPONSE` members, and monocyte markers (`CD14`, `LYZ`, `S100A8/9`, `VCAN`).
6. Only then, ML: **leave-one-study-out cross-validation** — train on all studies but one, test on the held-out study. Feature selection must happen inside the training fold. Report the per-study AUCs, not a pooled number.

**Hard requirements:**
- Whole blood and PBMC are different tissues (whole blood includes neutrophils). Do not pool them; treat tissue as a stratifying variable.
- Microarray and RNA-seq have different dynamic ranges. Meta-analyse effect sizes; do not concatenate expression matrices.
- Class imbalance: always report the majority-class baseline accuracy alongside any model accuracy.
- If a classifier reports AUC > 0.95, assume leakage until proven otherwise.

## Project 2 (after Project 1): single-cell meta-analysis and deconvolution

**2a — Deconvolution.** Use the annotated GSE277117 object as a cell-type reference (23 types, 314k cells) to deconvolve the bulk cohorts from Project 1 into estimated cell type proportions (CIBERSORTx, MuSiC, or BayesPrism). Then test: **does the CD14 monocyte decrease we found at n=25 reproduce at n~150?** This is the single most valuable follow-up, because it tests our one solid finding at a sample size that could settle it.

Caveats: deconvolution is reliable for abundant, distinct cell types and poor for rare or similar ones. Validate the reference by deconvolving held-out pseudobulk from the single-cell data itself and checking the estimates match the known proportions.

**2b — Single-cell meta-analysis.** Assemble public AS scRNA-seq studies that each contain cases and controls, and meta-analyse cell type composition and per-cell-type pseudobulk. Starting points: GSE194315 (10 AS vs 29 HC, CITE-seq PBMC — note its published classifier reports AUROC up to 1.00 on a test set of 5 AS patients and should be treated sceptically), GSE163314 (colon tissue AND blood; Crohn's, spondyloarthritis, co-morbid, and healthy — small but the only paired gut/blood design found), GSE288581 (AS synovial fluid with paired blood). Expect each study to be small; the point is consistency of direction across studies, not any single p-value.

## Checklist of traps — each of these produced a confident wrong answer in the last project

1. **Pseudo-replication.** Cells are not independent samples; patients are. Any label attached to a patient must be tested at patient level, with one row per patient.
2. **Hidden non-gene features.** GSE277117 mixed 163 CITE-seq antibody-capture features into the gene matrix, and assay type was ~perfectly confounded with the outcome. Filter on `var["feature_types"] == "Gene Expression"` **by index, not by name** — several antibodies share a symbol with a gene. Always print the feature-type breakdown of any new dataset before using it.
3. **Batch correction can delete the signal.** Never use a batch key that is close to 1:1 with your biological unit of interest. Harmony on `Sample_num` erased between-patient interferon variation; removing it produced sub-clusters that were up to 88% a single patient. Prefer clustering-free tests (score cells, threshold, compare patients) when the question is about between-patient differences.
4. **Test the claim that was actually made.** The paper's findings were about proportions of sub-populations identified by re-clustering *within a compartment*, not mean expression across whole cell types. Read the figures before designing the test.
5. **FDR does not catch systematic bias.** Multiple-testing correction filters noise. It cannot detect a normalisation artefact that shifts every gene in one group.
6. **Report power.** State the detectable effect size for the n you have, before interpreting a null.

## Environment

- Windows 11, PowerShell. Local venv `.venv`; Python 3.13.
- Modal for anything that loads a full single-cell object (laptop has 14.8 GB RAM; the objects need 32–64 GiB). Use `modal run --detach` so a network drop doesn't kill the job.
- Modal costs so far: a few dollars total. Volume storage ~$0.45/month.
