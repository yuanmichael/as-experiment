# Handoff prompt — applying the AS method to IBD

*Paste this into a fresh Claude Code session started in `C:\Coding\as-experiment`.*

---

I'm a patient with ankylosing spondylitis learning transcriptomics. I've just finished two projects in this repo and want to apply the same method to inflammatory bowel disease.

**Read first, in this order:** `README.md`, then `BULK_META.md` (the full methods and every caveat), then `scripts/geo_search.py`, `scripts/bulk_prep.py`, `scripts/bulk_meta.py`, `scripts/deconv_validate.py` and `scripts/deconv_bulk.py`. `PREPRINT.md` is the draft write-up of the AS work. `REPORT.md` is the earlier single-cell project.

**How I'd like you to work with me:** act as a tutor. Do the grunt work yourself — data cleaning, pipeline code, file prep — but surface the high-level decisions and check my understanding with questions rather than handing me conclusions. Explain from fundamentals. Grade my answers honestly, including when I'm wrong. Push back when I'm reaching.

## Why IBD, and why now

The AS project hit a wall that is about biology, not statistics. Everything we could find is **blood**, and AS pathology is in the entheses, sacroiliac joints and gut. Our best result — NK cells in AS blood make less cytotoxic cargo — reaches AUC 0.77 against healthy volunteers and a positive predictive value of 17.7% in a back-pain clinic. That is not a diagnostic, and no amount of better statistics on blood will make it one.

IBD does not have that problem:

- **Diseased tissue is routinely biopsied**, so public cohorts contain the actual inflamed tissue.
- **Inflamed and uninflamed biopsies from the same patient** exist. That is a within-patient contrast — the single strongest design available, and one AS can never offer.
- **The cohorts are an order of magnitude larger.** A GEO sweep (run with the generalised `scripts/geo_search.py --disease ibd`, output in `data/ibd_geo_candidates.tsv`) returned 650 human expression/methylation series, 113 with case and control hints at n ≥ 40.
- **There is a real gut–joint axis**, so this is a continuation of the AS question rather than a change of subject. Roughly half of AS patients have subclinical gut inflammation.

Candidates the sweep surfaced, unverified — **check every sample label directly, GEO summaries are unreliable**:

| Accession | n | What it looks like |
|---|---|---|
| GSE193677 | 2,490 | adult IBD biopsy cohort (MSCCR) |
| GSE100833 | 1,717 | ustekinumab Crohn's trial, blood **and** biopsy |
| GSE186507 | 1,030 | adult IBD blood cohort |
| GSE177044 | 1,035 | whole-blood RNA-seq, PSC vs UC |
| GSE66407 | 368 | CD, UC and healthy gut biopsies |
| GSE57945 | 322 | paediatric Crohn's ileal transcriptome (RISK) |
| GSE179285 | 254 | **active and inactive lesions, ileum and colon** — the paired design |
| GSE112366 | 388 | Crohn's, intestinal epithelium |
| GSE111889 | 251 | HMP2/iHMP longitudinal multi-omics |

## What to do

**Project 1 — decide the question before touching data.** With this much data the binding constraint is no longer sample size, so the AS framing ("is there a difference?") is the wrong one. Better candidates, in rough order of how much I'd learn:

- **Within-patient inflamed vs uninflamed.** Removes between-person variation entirely. What is different about tissue that is inflamed in a person whose other tissue is not?
- **Does the AS blood signature appear in IBD blood?** We have a defined NK effector signature and a platelet/integrin signature. Testing them in IBD blood is a direct test of the gut–joint axis, and it is cheap because the signatures are already specified. Pre-register them; do not go fishing.
- **Tissue vs blood in the same patients.** GSE100833 has both. How much of what happens in the gut is visible in blood at all? That question generalises straight back to why the AS project was limited.

Ask me which before building. I lean toward the third.

**Project 2 — deconvolution, if the question needs it.** The AS pipeline works, but **the single-cell reference must change.** Ours is PBMC and AS-derived. For gut tissue you need a gut reference — the Gut Cell Atlas or a public IBD scRNA-seq study with healthy controls. Do not reuse the AS reference on biopsies.

## Rules carried over — each of these produced a confident wrong answer

1. **Verify sample labels from GSM records, not series summaries.** In the AS project, GSE39428 looked like the largest cohort in the field (51 cases, 163 controls) and is a SNP genotyping array. GSE25101 annotates itself as both whole blood and PBMC.
2. **Never merge expression matrices across studies.** Compute the effect inside each study, meta-analyse the effect sizes. Published AS analyses merged two cohorts and compared one study's cases to another's controls.
3. **Depth-normalise count data.** Skipping CPM gave one AS cohort a mean effect size of −1.3 × 10¹², because library size tracked disease status and every gene inherited the shift.
4. **Check what the submitter already did to the matrix.** One AS cohort's GEO matrix had a supervised batch correction (DWD) already applied, which can delete or create the contrast being tested.
5. **Calibrate gene sets against the genome-wide background**, not against zero. The AS background was 36.4% of genes up, median *g* = −0.135. Against a 50/50 null every downward set looked significant.
6. **Validate deconvolution before interpreting it.** It always returns a confident number. Deconvolve held-out pseudobulk with known proportions first. In AS this killed the CD8 T result (*r* = 0.57) despite *P* = 0.029.
7. **Watch for circularity.** Adjusting an expression signature for a proportion estimated from the same transcripts is guaranteed to return "it's all composition". Use disjoint gene sets.
8. **An effect size is not a classifier.** The AS interferon panel was significant at *P* = 0.004 with AUC 0.496.
9. **Report power before interpreting a null.** State the detectable effect size for the n you have.
10. **Cells are not independent samples; patients are.** One row per patient.
11. **Translate to PPV at realistic prevalence** before calling anything diagnostic, and benchmark against what clinicians already use (for IBD: faecal calprotectin, CRP, endoscopy).

## Environment

- Windows 11, PowerShell. Local venv `.venv`, Python 3.13 — already has scanpy, anndata, statsmodels, scikit-learn, modal.
- Modal for anything loading a full single-cell object (laptop has 15 GB; objects want 32–64 GiB). `modal run --detach`. App `as-experiment`, volume `as-data`. Costs so far: a few dollars.
- Bulk analyses run locally in minutes.

## Repo conventions

Keep IBD work in `scripts/ibd_*.py` and `data/ibd_*`, so it does not tangle with the AS files the preprint depends on. The generic pipeline (`geo_search.py`, `geo_inventory.py`, `bulk_prep.py`, `bulk_meta.py`, `deconv_*.py`) is shared — parameterise it rather than copying it. `geo_search.py` already takes `--disease ibd`.

Commit as you go. The repo is public: https://github.com/yuanmichael/as-experiment
