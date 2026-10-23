# Ankylosing spondylitis: are NK cells missing, or just quiet?

Independent reanalysis of public AS transcriptomics, done from scratch. Written by a patient learning the methods, so the wrong turns are documented alongside the results.

AS blood has a weaker NK-cell gene signature than healthy blood. That has been published several times. What nobody had separated is **why**: are there fewer NK cells, or are the same NK cells making less?

Those mean different things. Fewer cells is a trafficking or survival problem. Quieter cells is a functional one.

## What we found

Five independent cohorts, **97 AS patients vs 70 healthy controls**, each study analysed separately and combined by random-effects meta-analysis. Then deconvolved against a 312,000-cell single-cell reference to split proportion from per-cell expression.

| | pooled effect (*g*) | 95% CI |
|---|---|---|
| NK cytotoxic cargo (granzymes, perforin, granulysin) | **−0.98** | −1.33 to −0.63 |
| ...after accounting for how many NK cells are present | **−0.80** | −1.17 to −0.44 |
| NK cell proportion | −0.54 | −0.91 to −0.17 |

**82% of the cytotoxic deficit survives adjusting for cell number.** NK cells in AS blood have turned down their weapons. There are also somewhat fewer of them, but that is the smaller effect.

![deconvolution](figures/deconvolution.png)

Independent support, from a method that cannot be fooled the same way: a [meta-analysis of 95 flow-cytometry studies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8291033/) covering **4,020 patients and 3,065 controls** found **no difference in NK cell percentage** at all. Flow cytometry counts cells by surface protein and is blind to how much a cell is transcribing. Our method infers cell counts from transcripts and so mistakes a quiet cell for an absent one.

The two have opposite blind spots and point the same way. If the flow data is right, the true composition shift is near zero and essentially all of it is regulation.

## How we avoided fooling ourselves

The obvious version of this analysis is circular. Deconvolution reads cell proportions out of the same transcripts as the expression signature, so "adjust the NK signature for the NK proportion" removes the signal by construction and would return "it's all composition" no matter what is true.

So the two numbers use **disjoint gene sets**. Proportion comes only from NK *identity* genes (receptors, lineage transcription factors), with every effector gene deleted from the signature matrix first. The outcome is the *effector* genes, untouched.

Count the soldiers by their uniforms, then ask how much ammunition they carry. Count them by their ammunition and you can never tell a smaller army from a worse-supplied one.

Other checks that changed the answer:

- **Validate before interpreting.** Deconvolution always returns a confident number. We tested it on held-out pseudobulk from 13 patients the reference never saw, where the truth is known. NK recovered at *r* = 0.83. **CD8 T failed at 0.57** and was discarded, despite a nominally significant result (*p* = 0.029), because cytotoxic T cells share GNLY, NKG7 and PRF1 with NK cells.
- **Check the background.** The genome-wide median effect is −0.13, not zero. Against a 50/50 null every downward gene set looks significant.
- **Two normalisation bugs found by diagnostics, not by FDR.** Raw counts used without depth normalisation in one cohort (library size tracked disease status, so every gene inherited the shift), and a supervised batch correction already baked into another cohort's GEO matrix.
- **An unexplained compartment.** 30–59% of each real bulk sample lands in a sink of rare reference types, and the sink is slightly larger in AS, which could manufacture the whole result. Renormalising moves NK from −0.54 to −0.49. Not the explanation.

## Also in there

Platelet and integrin genes are the strongest *upward* signal: 10 of 11 canonical markers up, mean *g* = +0.44 (ITGA2B +0.83, *q* = 0.001). A [2025 deep plasma proteomics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12957791/) of 104 samples independently landed on integrin-pathway proteins as its top AS biomarkers. Different molecules, different compartment, same pathway. We have not followed this up.

**AIF1 is a clean null** (*g* = +0.28, CI −0.18 to +0.74) at 80% power for *g* = 0.44. Not a failed replication: AIF1 was established as a TNF-inhibitor *response* marker, a within-disease contrast this data cannot address.

## What this is not

A diagnostic. The NK score reaches AUC 0.77 against healthy controls. In a back-pain clinic, a positive result means a **17.7% chance of AS** — 4.6 false alarms per true case. CRP manages 15.6% and HLA-B27 manages 37.2%. All three are rule-out tests, not rule-in.

And every control here is a healthy volunteer. The comparison that matters clinically is AS versus mechanical back pain, and no public dataset we could find has it.

## What we are building now

The result above is a graded shift in host-cell expression, and graded shifts make poor diagnostics — hence the 17.7%. The data type that could do better is the **T-cell receptor repertoire**, because a receptor is a discrete sequence rather than a level. Faham et al. (*Arthritis Rheumatol* 2017) reported a CDR3 motif in 30–39% of AS cases, and [Komech et al.](https://doi.org/10.1093/rheumatology/kex517) found an HLA-B27-associated clonal cluster in 15 of 16 synovial samples.

Both results need reading carefully. Faham's headline compares a mostly-B27-positive case group against mostly-B27-negative controls; the comparison that isolates the motif from the HLA type, B27+ AS versus B27+ healthy, gives *p* = 0.049. So the incremental value over a routine B27 test is far smaller than the headline implies.

A count matrix cannot hold a sequence, so none of our existing data can test this. TRUST4 on bulk RNA-seq recovered **146–498 TCRβ reads per sample, 0.0009% of reads** — roughly 200 clonotypes out of a repertoire of 10⁶. The assay, not the analysis, is the limit.

### The data exists, and our first inventory said it didn't

An AIRR Data Commons sweep returned 5,356 repertoires and **zero** for AS, IBD or RA, and we concluded a repertoire model could not be built from public data. **That was wrong.** The ADC indexes curated repertoire studies; the archives hold the rest, filed as ordinary sequencing. A sweep of GEO, ENA and ArrayExpress found 2,715 candidate studies, and for AS specifically six — including Komech's two cohorts (44.5 GB and 47.1 GB of dedicated 5′RACE TCRβ) and paired single-cell data in psoriatic arthritis.

The systematic blind spot: **10x 5′ libraries with V(D)J enrichment get indexed as "single-cell RNA-seq."** They are dedicated, deep TCR sequencing *with a paired transcriptome*, and no repertoire inventory will surface them. Of 450 series checked, 86 have both — 6,008 samples.

### Built: 1.86M cells with a receptor and a transcriptome

Paired data makes one thing possible that sequence-only data does not: **one cell is one training example**, its own receptor as input and its own expression as the label. That inverts the scarcity we hit everywhere else — 167 patients in the meta-analysis, 666 in the largest HLA-typed repertoire cohort, against millions of cells.

| | |
|---|---|
| libraries built | **291** across 7 studies |
| cells with a productive TRB **and** a transcriptome | **1,862,114** |
| distinct clonotypes | 1,421,133 |
| **singleton clonotypes** | **91.3%** |
| CITE antibody features caught and dropped | 8,831 |
| held out, never trains | GSE288581 — AS, 27,375 cells, 85% pairing |

The singleton fraction is the number that matters. Expanded clones are overwhelmingly effector or exhausted cells, so a model graded on them is being graded on "big clones are effectors," which is textbook. Cells whose receptor appears exactly once are the honest stratum, and they are the overwhelming majority rather than a rump.

Two studies were dropped after measurement: liver-transplant biopsies yielded 3,029 cells from 30 libraries at 1.1% median pairing, and a CD45-sorted cohort gave 28,589 at 4.6%. Not wrong, just mostly not T cells.

### Three pipeline bugs that each looked like a result

The pattern worth recording, because it recurred:

- GEO's `filelist.txt` is `Archive/File | Name | Time | Size`, so the name is the **second** column. Reading the first gave every file the name "File". Paired-series count: 2 → 11 → **86**.
- One study deposits `Bac17L_CD3pos` and `Bac17L_CD3neg`. Pairing the TCR file against the CD3-**negative** barcodes gave **0.8%** overlap, which reads exactly like poor TCR recovery. Correct pairing: 83.3%.
- A low pairing rate was then flagged as a bad join — but in a T-sorted library the *unpaired* cells are also T cells, so the marker-enrichment test returns ~1 for a perfect join. One number almost never identifies a failure; the rate and the enrichment are both recorded for every library so they can be audited rather than trusted.

### Modelled: receptor sequence → that cell's own transcriptome

Two stages. Stage 1 compresses 2,009 genes into a learned latent using expression alone — no receptor, no gene lists, no labels. Stage 2 predicts that frozen latent from the CDR3 plus V and J genes. Stage 1 is frozen deliberately: trained jointly, the shared latent would be dominated by the expression branch and the receptor branch would wither while still looking load-bearing, which is [plapt-ablation](https://github.com/yuanmichael/plapt-ablation)'s `prot_zero` result.

Every measurement is run twice — once with real receptors, once with receptors permuted between cells of the same subject. **The gap between the two is the signal**; the absolute value is not interpretable.

| | real | shuffled |
|---|---|---|
| exhaustion, AS test, **per cell** | +0.247 [0.22, 0.29] | −0.003 [−0.01, 0.01] |
| exhaustion, AS test, **per clonotype** | **+0.126** [0.08, 0.18] | +0.010 [−0.03, 0.06] |
| exhaustion, held-out **cohort**, per clonotype | +0.012 | +0.005 |
| variance of the latent explained (R²) | 0.009 val, 0.022 AS | — |

So: a receptor explains on the order of **1–2%** of its own cell's expression state, and on a cohort the model has not seen, approximately none.

Four things learned the hard way, all of which inflated earlier numbers:

- **Half the cell-level signal was pseudo-replication.** Cells in one expanded clone share a receptor *and* a history, so they are not independent observations. Counting each clone once took exhaustion from 0.247 to 0.126 — the same error family as the KLRB1 result that read *p* = 0.0014 per cell and *p* = 0.289 per patient.
- **The loss was never the bottleneck.** The first design predicted raw counts and showed no real-vs-shuffled gap in its loss, so we rebuilt around a learned latent on the theory that 2,009 mostly-zero genes drowned the signal. The rebuild scored 0.247 against the original 0.244. Identical. The objective was not the constraint; a targeted probe found the signal either way.
- **Our marker scores were measurably broken.** Each was a mean of ~9 normalised genes, which only measures one thing if those genes co-vary. CTLA4 sat in both the exhaustion and regulatory sets, making those two scores correlate at **0.451**; CD160 and BTLA loaded at 0.089 and 0.075 on the exhaustion axis, i.e. barely participating. Programme coherence (mean pairwise *r* between its genes) predicts recoverability: cytotoxicity 0.397 → recovered at 0.86, regulatory 0.150 → 0.70. Scores are now PC1 of each set rather than its mean.
- **Compression discards sparse programmes.** A round trip through the frozen stage-1 encoder and decoder — no probe in the path, no linearity assumption — recovers naive/memory at 0.81 and exhaustion at only 0.60. Whatever stage 1 loses, stage 2 can never recover.

### What is not established, and the gap that matters most

**There is no positive control.** Every result above is a negative or near-negative measured against a shuffled baseline. That establishes "not much signal." It cannot establish that the pipeline is *capable* of finding signal that is present — a bug in tokenisation, the splits, the probe or the latent would produce these same numbers, and nothing here would reveal it.

Until a target known to be receptor-encoded is recovered by this pipeline, every negative in this section is uninterpretable. Chasing better losses and better targets before running that control was the wrong order of work.

**Resolved.** The controls were run and the pipeline is sound:

| control | val | held-out cohort | AS test |
|---|---|---|---|
| V gene from the latent (plumbing — V gene is an *input*, so this is circular by design) | 0.9998 | 0.9999 | 0.9999 |
| **CD4 vs CD8**, model / plain features | **0.658** / 0.641 | 0.528 / 0.527 | — |
| **MAIT (SLC4A10+)**, model / plain features | 0.612 / 0.603 | 0.507 / 0.510 | **0.802** / 0.804 |

CD4-vs-CD8 at 0.658 is squarely where the literature puts beta-chain-only classifiers, and MAIT reaches 0.802 on the AS cohort. Known receptor→phenotype relationships are recovered at the magnitudes they should be. **So the 1–2% programme signal is a property of the biology, not a bug.**

Two things fell out of those controls that matter more than the controls did.

**Plain features match the transformer.** Logistic regression on V-gene one-hot, J-gene one-hot, CDR3 length and amino-acid composition — no learning, no attention, no sequence order — ties a 1.5M-parameter transformer to within 0.02 AUC on every control and beats it on one. So the receptor information in play is *gene usage and composition*, not sequence order. (This is [plapt-ablation](https://github.com/yuanmichael/plapt-ablation)'s finding for a second time: a 4.3× larger encoder moved that score by −0.0016.)

**Nothing transfers across cohorts — not even the known biology.** CD4-vs-CD8 falls from 0.658 to 0.528 on a held-out study; MAIT from 0.612 to 0.507. Both arms, both feature sets. A relationship as solid as T-cell lineage essentially vanishes when the cohort changes. That is a larger obstacle than the weak programme signal, and the weak programme signal is plausibly a symptom of it.

### Where this is going

The cell-level model answers *does a receptor shape its own cell's state*. That is interesting biology and **not the question this project set out to ask**, which was whether a blood draw can predict disease — patient level, the way Faham's clonotype result is patient level. The cell-level route was adopted to work around label scarcity (millions of cells against ~666 labelled patients) and became the work.

Three results from it bear on the real question, and all three simplify it: per-cell receptor→state coupling is 1–2%, so a patient-level predictor should **not** be routed through predicted cell states; plain features match deep learning, so the patient-level version needs neither GPUs nor a large model; and [Mal-ID](https://www.science.org/doi/10.1126/science.adp2407) already demonstrates the patient-level approach working on lupus, type 1 diabetes, COVID and HIV.

The next step is in [NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md): patient-level classification from already-processed public repertoires, validated against a published benchmark before anything is trusted. **That is now built — see [Patient-level repertoires](#patient-level-repertoires-the-first-thing-that-transferred) below**, where CMV transfers across cohorts at 0.851 and the AS clonotype replicates in an independent cohort.

The architecture is prior art — [TESSA](https://doi.org/10.1038/s41592-020-01020-3), [mvTCR](https://www.nature.com/articles/s41467-024-49806-9), [scNAT](https://link.springer.com/article/10.1186/s13059-023-03129-y) all build joint TCR–transcriptome embeddings. What appears not to have been done is transferring such a model to a disease it never trained on, and reading it out on sequence-only AS repertoires.

And a live warning: the closest prior work ([Wang & Ji](https://www.biorxiv.org/content/10.1101/2023.11.27.568912v2)) retitled from "TCR sequences **predict** gene expression" to "transcriptional landscapes are **shaped by** TCR sequence similarity" between versions, tripled its data to 2.2M cells, and remains unpublished. The prediction claim may not survive an honest split. Finding that out is the point.

## Patient-level repertoires: the first thing that transferred

The plan above is now built, and the headline is that **patient-level repertoire classification works where the cell-level route did not.**

### The data, and a second inventory error

[NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) recorded that Emerson's CMV cohort needed immuneACCESS registration. **Wrong again** — it is in the AIRR Data Commons as `emerson-2017-natgen`, served by VDJServer. The same blind spot as before, in the opposite direction: we assumed a dataset was gated because one index didn't surface it.

One trap worth recording: the sync ADC endpoint caps every page at 1,000 rows, which puts 786 repertoires at **29 hours**. The async endpoint returns a whole repertoire as one file in under 10 seconds.

| | |
|---|---|
| repertoires fetched | **786** (666 Cohort 01 + 120 Cohort 02) |
| clones | **178,045,055** |
| median clones per repertoire | 222,256 (min 1,354, max 670,400) |
| HLA-A/B typed subjects | **626**, all in Cohort 01 |
| CMV serostatus | 340 positive, 421 negative, 25 unknown |
| public CDR3s (seen in ≥7 subjects) | 2,335,568 |

CMV labels are not in the ADC copy; they come from a Hugging Face mirror of the original metadata. The join was verified rather than assumed: **786/786 subject IDs matched, and HLA-A\*02 agreed with the ADC's own typing in 626/626 subjects with zero disagreements.**

### Controls first, because the batch effect here is enormous

| target | both blocks | public clones only | summary only | technical only |
|---|---|---|---|---|
| **HLA-A\*02** (positive control) | 0.953 [0.935, 0.969] | **0.959** | 0.607 | 0.476 |
| **cohort membership** (negative control — no biology exists) | 0.994 | 0.699 | 0.997 | 0.801 |
| **CMV**, cross-validated | 0.898 [0.875, 0.919] | 0.872 | 0.781 | 0.623 |
| **CMV, held-out cohort** (train 01 → test 02) | **0.851** [0.778, 0.915] | 0.842 | 0.769 | — |

How to read these. HLA-A\*02 is carried by 47% of these subjects and is directly visible in which receptors a person has, so 0.95 is the ceiling this data supports — it says the pipeline works. Cohort membership is the floor: the two cohorts differ by **nothing biological**, so anything above 0.5 is pure batch. It reaches **0.994**.

That one number reframes everything. Summary features (gene usage, clonality) hit 0.997 on batch — they are a *sequencing-run detector*. Only public-clone identity has a usable ratio: 0.959 on real biology against 0.699 on batch.

**CMV transfers across cohorts at 0.851.** This is the first cross-cohort transfer this project has achieved; the cell-level model's known-biology controls collapsed from 0.658 to 0.528 and 0.612 to 0.507. Age and sex alone give 0.584, so it isn't demographics. The permutation null sits at mean 0.508 over 10 shuffles, and with 10 permutations the smallest reportable *p* is 1/11 = 0.091 — that is a floor imposed by the number of shuffles, not a weak result.

One correction made along the way: a "depth-only baseline of 0.772" was reported as a possible artefact, having mislabelled Shannon entropy as depth-sensitive. Diagnostics showed the true technical measures (`n_templates`, `n_clones`) reach only 0.526–0.529 single-feature AUC (*p* = 0.18–0.22), while clonal expansion carries it (`top10_share` 0.798, clonality 0.761, *p* = 3.3 × 10⁻³⁵). That is CMV memory inflation — real biology, not depth.

### The AS clonotype, replicated in an independent cohort

`CASSVGLYSTDTQYF` (TRBV9/TRBJ2-3) is Komech's AS clonotype and carries Faham's `SVGLYSTDTQ` motif.

In **GSE288581** (4 AS patients, sorted CD45RO+CD8+) it appears in 3 of 4 AS and 0 of 4 controls. Against the 786 Emerson healthy subjects it is present in **15/786 (1.9%)**. Two other candidate CDR3s from the same screen were refuted outright at 50.9% and 79.5% healthy carriage.

**It is HLA-B27-restricted in healthy people**: 6 of 44 B27-positive subjects (13.6%) versus 7 of 582 B27-negative (1.2%), *p* = 1.0 × 10⁻⁴. So the comparison that isolates disease from HLA type is against B27+ healthy controls, not all healthy controls — the same correction Faham's headline needed.

Then two independent Komech cohorts, 5′RACE bulk from a different lab with a different assay, as the test.

**`PRJNA393498` (Komech 2018) was written off as unlabelled and is not.** Our own handoff doc recorded that its disease labels survive only inside FASTQ filenames, so only ~5 AS patients looked identifiable. ENA's mirror does serve an empty `disease` field — but NCBI's BioSample XML carries submitter-supplied `disease`, `health_state`, `Patient No`, `cell_type`, `tissue`, `age` and `sex` for all 78 runs. Verified against NCBI directly: **52 AS runs across 26 patients, 26 healthy runs, every label documented and none inferred.** That took the cohort from 7 AS patients to 31.

Two of its subsets must be excluded, and the first is the most dangerous mistake available in this dataset:

- The **14 TRBV9-sorted runs** are sorted for the V gene our target uses, and they are *all healthy donors*. Including them inflates healthy carriage and reads as a refutation of the finding.
- The **24 CD4/CD8-sorted AS runs** sample a different slice of the repertoire than whole PBMC. They are partial views of patients we already have, not extra patients — and a CD4-only library cannot contain a CD8 clone, so counting one as a negative manufactures a false negative.

Unsorted PBMC only, aggregated per subject, never per run.

### The finding is clonal expansion, not presence

The carrier frequencies barely overlap:

| group | frequency in carriers |
|---|---|
| healthy B27+ (Emerson, 6 carriers) | 1.7 × 10⁻⁶ – **9.7 × 10⁻⁶** |
| PsA (4 carriers) | 5.4 × 10⁻⁶ – **4.3 × 10⁻⁵** |
| **AS (12 carriers)** | **1.6 × 10⁻⁴ – 4.5 × 10⁻²** |

Healthy B27+ people and PsA patients carry this clonotype as a trace of a few templates. AS patients carry it **expanded 10- to 1,000-fold**, up to 4.5% of the entire repertoire in one patient. There is a clean gap between 4.3 × 10⁻⁵ and 1.6 × 10⁻⁴ with nothing in it. So presence/absence was measuring the wrong thing: the clonotype is not absent from healthy B27+ blood, it is present and small.

Carriage at frequency ≥ 10⁻⁴, restricted to subjects sequenced deeply enough to detect that level:

| group | carriers | 95% CI |
|---|---|---|
| **AS (2018 + 2022)** | **9/25 = 36.0%** | 18–57% |
| PsA (2022) | **0/16** | 0–21% |
| healthy (2018 study) | 0/8 | 0–37% |
| healthy B27+ (Emerson) | **0/44** | 0–8% |
| healthy, all (Emerson) | 0/785 | 0–0.5% |

| comparison | *p* (Fisher, one-sided) |
|---|---|
| AS vs healthy B27+ | **3.6 × 10⁻⁵** |
| **AS vs PsA** — same lab, same assay | **5.8 × 10⁻³** |
| AS vs all healthy | 5.2 × 10⁻¹⁵ |

**AS vs PsA is the load-bearing comparison**, because both arms are 5′RACE from one lab on one protocol — the only comparison here immune to the 0.994 batch effect. On presence/absence it was *p* = 0.16; on expansion it is *p* = 0.0058. **The expansion is specific to AS, not to spondyloarthritis generally.** And 36.0% lands inside Faham's reported 30–39%, reached from different data by a different method.

32 AS patients were scanned; 25 are evaluable at 10⁻⁴ and 7 are too shallow to be counted either way. Presence at any level, which is depth-confounded in the opposite direction, is 17/32 = 53.1% (CI 35–71%).

**Why it is not an artefact:** the AS carriers have *different nucleotide sequences* encoding the same peptide (`GTCGGCCTC`, `GTGGGACTA`, `GTAGGTCTC`, `GTAGGCCTT`). Contamination or index hopping would give identical nucleotides. Independent rearrangements converging on one protein is the signature of antigen-driven selection.

### Three caveats, one serious

**The 10⁻⁴ threshold was not pre-specified.** It was chosen after seeing the depth structure, which is selection. Mitigating it: the result holds across the whole plausible range — *p* = 0.0036 at 3 × 10⁻⁵, 4 × 10⁻⁵ at 10⁻⁴, 2 × 10⁻⁴ at 3 × 10⁻⁴, 10⁻⁴ at 10⁻³ — so it is not knife-edge. Read the expansion *p*-values as exploratory; the pre-specified presence test (*p* = 0.0033) stands on its own.

**RNA versus DNA confounds the magnitude claim.** The patient data is 5′RACE on RNA; Emerson is ImmunoSEQ on genomic DNA. One cell gives one DNA template but a variable amount of mRNA, and activated cells transcribe more TCR — so an activated expanded clone looks larger in RNA than in DNA. This does not touch AS vs PsA, both RNA from one lab, which is the other reason that comparison carries the weight.

**The 2018 libraries are ~20× shallower** in CDR3-spanning reads (1.4k–216k vs 1–4M), so 7 of 32 AS subjects are not evaluable at 10⁻⁴. Their zeros are excluded rather than counted as negatives; a subject with 15,000 CDR3 reads can only rule out frequencies above ~7 × 10⁻⁵.

### As a diagnostic: two operating points

| | presence (any level) | expanded (≥ 10⁻⁴) |
|---|---|---|
| sensitivity | 53.1% (17/32) | 36.0% (9/25) |
| specificity vs all healthy | 98.1% | **100%** (0/785) |
| specificity vs healthy B27+ | 86.4% | **100%** (0/44) |
| specificity vs PsA | 75.0% | **100%** (0/16) |
| PPV, back-pain clinic (5%) | 60% | 81–100% |
| PPV, B27+ back-pain clinic (30%) | 63% | 67–100% |
| NPV, B27+ back-pain clinic | 82% | 77–79% |

PPV ranges for the expanded point are given best-case to worst-case because **zero false positives were observed**, so the point estimate of 100% is not credible; the lower figure uses the 95% upper bound on the false-positive rate (0.47% against all healthy, 8.0% against B27+). For comparison, in the same back-pain-clinic setting CRP gives 15.6% PPV, HLA-B27 37.2% and our NK expression score 17.7%.

The threshold trades sensitivity for specificity, and expansion is the better trade: it is a **strong rule-in test and a poor rule-out test**, missing roughly 60% of patients. HLA-B27 status for the AS patients is in no structured metadata anywhere, which is the main thing worth asking the authors for.

### The $25 detour: assembly was never necessary

Three attempts, two of them wasted, for one lesson.

TRUST4 is built for bulk RNA-seq, where ~1% of reads are TCR and its first step discards the rest before assembling. These are **5′RACE amplicon** libraries where nearly 100% of reads are TCR, so nothing was discarded and it tried to assemble ~100× its design load. Assembly cost grows faster than linearly: one library ran past six hours unfinished. The tell was in plain sight — **325,158 productive clonotypes from 5.4M reads**, where RNA-seq yields hundreds.

Both failures were deadlines, not crashes:

1. The **orchestrator** had `timeout=7200`. It expired while its 8 workers were mid-assembly and Modal cancelled all of them. 2 of 51 runs survived.
2. The orchestrator limit was raised, but each worker kept a 6-hour cap. A worker hit it, and because `.map` re-raises by default, that single timeout killed the whole batch. 2 of 12 survived.

**The claim that made the second failure possible:** after the first, the code was described as unable to take down a batch, because every worker wraps its work in `try/except`. That is only true for Python exceptions. A container killed for exceeding its time limit is terminated from outside the interpreter — there is nothing to catch. One half of the problem was fixed and the whole thing declared solved. The fix is `.map(..., return_exceptions=True)` plus a per-run cap set from *measured* throughput.

What worked: **stop assembling.** The question was presence/absence of one sequence, which needs no repertoire catalogue. A nucleotide regex built from the peptide's codon alternatives finds it in any frame at compiled speed.

| | assembly | read-level scan |
|---|---|---|
| time per library | >6 h (unfinished) | **1.4–3.6 min** |
| cost, 51 libraries | ~$25 for 4 reports | **~$0.50 for 51** |
| detection floor | needs a contig to survive filters | **1 read in 7 million** |

The scan is also strictly *more* sensitive: one AS patient at 1.4 × 10⁻⁶ would have been invisible to the assembler. A `CASS` count per library is carried alongside, so a zero distinguishes "absent" from "scan broken" — without it, zero means nothing.

The general lesson: a full repertoire catalogue was built when the hypothesis only needed Ctrl-F. The cost came from reaching for the standard pipeline instead of asking what the question required.

### Where this is going: the representation, tested on a target we can actually evaluate

[Mal-ID](https://www.science.org/doi/10.1126/science.adp2407) (*Science* 2025) reaches 0.986 multi-class AUROC over 593 people with a three-model ensemble, and one of the three is **off-the-shelf ESM-2 (150M), not fine-tuned**, embedding CDR3s for a patient-level classifier. Our own exact-match features give 0.898/0.851 on CMV, so the open question is whether a learned representation beats exact matching — and specifically whether it beats it *on held-out cohort transfer*, which is where everything here has previously died.

AS cannot answer that question at n=7. **CMV can**, with 761 labelled subjects and a built-in cross-cohort split. So the benchmark holds the task and splits fixed and swaps only the representation: exact-match public clones (0.898/0.851, done) against edit-distance clusters, ESM-2 embeddings, and a TCR-specific encoder. Whichever wins on cohort transfer earns the right to be pointed at AS.

Two notes on expectations. Mal-ID's TCR clustering model underperformed, which the authors attribute to **not modelling HLA genotype** — we have HLA-A/B for 626 subjects, which they did not. And their known-binder enrichment for TCRs reached only 0.56 AUROC, which is a reason to expect little from curated TCR–antigen databases.

### Widening the search destroys it — twice

Two attempts to generalise from one sequence to a motif family, both negative, both informative.

**Edit distance.** Within one amino acid of the target, **584/786 healthy subjects (74%)** are carriers; within two, 783/786. Restricting to TRBV9 only drops it to 57%, not enough to matter. In the patients the same saturation: within one edit, AS 6/7 and PsA 13/16. Against healthy B27+ controls the comparison goes from ***p* = 0.0033 at the exact sequence to *p* = 0.58 one edit away.** The signal does not degrade gradually; it is gone at radius 1.

The reason is generation probability. `CASSVGLASTDTQYF` is carried by 147 healthy people and `CASSVGLSTDTQYF` by 114 — one edit out and 10–20× more common than the target, because TRBV9 + TRBJ2-3 with minimal junctional insertion produces `CASSVGL·STDTQYF` shapes readily. The target sits inside a dense cloud of easily-made public sequence, and the discrimination lives in the exact residues. All 15 healthy carriers of the exact sequence use TRBV9/TRBJ2-3, 15 for 15.

**ESM-2 embeddings.** [Mal-ID](https://www.science.org/doi/10.1126/science.adp2407) reaches 0.986 AUROC using off-the-shelf ESM-2 150M on CDR3s, so this was not a long shot. 3,959,135 CDR3s embedded, scored by cosine to the target, with the probe pre-registered before the sweep. It fails: the nearest neighbour is `CASSVGLSTDTQYF` (114 healthy carriers) at 0.9979, and the whole pool spans only 0.969–0.998 — an unrelated CDR3 sits about as close as a one-edit neighbour. The pre-registered *k* = 50 family covers 54.6% of healthy subjects and performs **worse** than exact matching (AS vs healthy B27+ *p* = 0.15 against 0.0033).

The top hits say exactly what it learned: `CASSVGLSYTDTQYF`, `CASSVLGYSTDTQYF`, `CASSVGYLSTDTQYF` — **anagrams**. Mean-pooling over residues discards order, so cosine measures amino-acid composition. A composition-matched anagram is not a functional relative, because TCR recognition depends on which residue sits where.

What this does *not* show: that the embeddings lack useful information. Mal-ID never used them as a similarity metric — they fed per-sequence vectors to a *supervised* classifier, which can learn which of the 640 dimensions matter while cosine weights all of them equally. The fair test is embeddings as features for a patient-level classifier, measured on CMV where 761 labelled subjects and a held-out cohort can produce a number. A TCR-specific encoder (TCR-BERT, SCEPTR, TULIP) or TCRdist3's position-weighted distance is the better starting point than a general protein model, since the anagram failure is direct evidence that position is what matters.

An honest caveat on widening the search: *p* = 0.021 above holds up because the sequence was pre-specified by two prior independent studies. Any embedding-based sweep has an effectively unlimited hypothesis space, so discovery has to run **AS versus PsA within one study** (assay-matched, immune to the 0.994 batch effect) with permutation nulls, and replicate against the 786 Emerson repertoires before anything counts.

### The assay sets the ceiling: bulk RNA-seq cannot test this

GSE221786 (= PRJNA916206) looked like a way to add 20 AS patients and 8 controls. It is **bulk RNA-seq, not TCR-seq** — `library_strategy=RNA-Seq`, `library_selection=cDNA` — but with a twist worth measuring rather than dismissing: the PBMCs were stimulated with CytoStim for 4 hours and enriched for IL-17 producers. Activated T cells transcribe far more TCR mRNA than resting blood, so the TCR content might have been lifted into usable range.

It was not. Measured on the two deepest libraries, 4M reads each:

| | AS_M_5 | AS_M_2 |
|---|---|---|
| TRBC anchor reads per 4M | 193 | 173 |
| TCR content | 0.0048% | 0.0043% |
| projected CDR3-spanning reads at full depth | 2,817 | 2,639 |
| **detection floor** | **3.6 × 10⁻⁴** | **3.8 × 10⁻⁴** |

The floor is ~4× coarser than the 10⁻⁴ threshold the expansion result is called at, so a patient sitting at 2 × 10⁻⁴ — comfortably expanded — would show **zero reads**. A table of zeros here would be uninterpretable, not negative. Depth cannot rescue it either: reaching 10⁻⁴ needs ~10,000 CDR3-spanning reads, i.e. ~190M reads per sample against the 70M deposited, and pooling is not available because the test is per-patient.

The interesting part is that the hypothesis was right in direction and wrong in size. Resting whole blood runs near 1.4 × 10⁻⁵ TCR content; these activated, IL-17-enriched PBMCs reach ~4.5 × 10⁻⁵. Activation enrichment bought a real **~3×**, not the ~10× required. That generalises: activation protocols do not turn bulk RNA-seq into a clonotype-frequency assay.

Cost of learning this: about 5 cents of measurement, after ~$1–2 lost to scanning two 72M-read libraries in full and hitting a 90-minute cap. Which was the wrong design — the quantity wanted is a *fraction*, and a fraction does not need the whole file. Streaming the first 4M reads and projecting onto the library's real depth answers it in ~4 minutes without the multi-GB download ever completing.

### A denominator that nearly produced a confident wrong answer

The first version of that pilot reported **0.50% of reads as TCR** — a 500-fold lift over resting blood, and a green light to spend real money on all 28 samples. It was noise, and the way it was caught is the generalisable part.

`CASS` is the four-residue framework motif used throughout this project as the denominator: a count of CDR3-spanning reads, so that a zero for the target can be told apart from a broken scan. Expanded through codon alternatives it is 12 nucleotides matching one position in ~58,000. A 150bp read offers 139 positions and both strands are scanned, so **278 × 1.7 × 10⁻⁵ = 0.48% of reads match with no TCR present at all.** The two libraries came in at 0.5055% and 0.5036%: on the chance line, and agreeing to within 0.4% of each other, which biology does not do between two people.

Replacing it fixed the measurement. The anchor is now an exact 25-mer, `GGCCACCTTCTGGCAGAACCCCCGC`, taken from a **173-nt stretch identical in TRBC1 and TRBC2** (Ensembl ENST00000633705 / ENST00000466254), so one pattern counts both constant-region alleles. Chance expectation is 2.2 × 10⁻¹³ per read: every hit is a real TCRβ transcript. A constant-region count is then converted to CDR3-spanning reads by window arithmetic — a CDR3 call needs the read to span ~45nt, giving a (150−45+1) start window against the anchor's (150−25+1), so CDR3 reads ≈ 0.84 × anchor reads. RNA-seq coverage leans 3′ and the constant region sits 3′ of CDR3, so the assumption errs towards caution.

**The obvious worry is whether this invalidates the project's own numbers, since they all use `CASS` as the denominator. It does not, and the check is worth recording.** Auditing all 39 scanned Komech 2018 libraries:

| | CASS as % of reads | spread | CV |
|---|---|---|---|
| Komech 2018 amplicon (39 libraries) | 0.022% – 0.966% | **44-fold** | **0.82** |
| RNA-seq pilot, background-dominated (2 libraries) | 0.5036% – 0.5055% | 1.004-fold | 0.003 |

A background-dominated count is nearly *constant* across libraries, because it depends only on read count and length. A signal-dominated count varies with the biology. The amplicon libraries vary 44-fold; the background-dominated ones agree to four significant figures. That settles it.

Two further points fell out of the audit. Most amplicon libraries sit **below** the uniform-random prediction, some 7× below — because amplicon reads are not random sequence. They are near-copies of a few hundred nucleotides of TRB with a variable CDR3, and a fixed sequence either contains `CASS` or does not; the random-sequence rate simply does not apply. Bulk RNA-seq samples the whole transcriptome, which is diverse enough that the uniform model held to within 5%. So the same motif is safe in one assay and worthless in the other, and which one it is depends on the sequence diversity of the library, not on the motif.

The lesson is not "degenerate motifs are bad". It is that **a denominator needs its own null**, and that two libraries agreeing too precisely is evidence of an artefact rather than reassurance. The pilot now also reports `CASS` as a multiple of its chance rate, and carries a third verdict — INCONCLUSIVE, when a library yields too few anchor reads to project from. The original had only USABLE and NOT USABLE, which is why a noise measurement came back as a confident green light.

### Synovial fluid: the compartment where the clone was first found

Running now, scope fixed in advance. PRJEB51315 holds more than the blood arm already analysed:

| | patients | runs | reads |
|---|---|---|---|
| AS SF | 7 (`Bel Kos Lap Non Pon Shor Zar`) | 17 | 86M |
| PsA SF | 17 | 47 | 245M |

Three features decided the scope. Only **3 of the 7 AS SF patients are new** (`Lap`, `Non`, `Shor`) — the other four also have blood, so they give a **within-person blood-versus-joint comparison** rather than a between-group one. `Lap` and `Shor` carry CD8-sorted libraries only, which is the compartment the target clone occupies, so they are pre-enriched rather than impoverished. And **PsA synovial fluid exists in the same deposit from the same lab on the same 5′RACE assay**, which matters most: there is no healthy synovial fluid anywhere, because healthy joints have no effusion to sample, so the only possible control for AS synovium is synovium from another arthritis. Without it an AS carriage figure in joints would have no denominator to be read against.

Synovial carriage will **not** be pooled into the blood denominator. The headline figure stays 9/25 at ≥ 10⁻⁴ in blood; SF is reported as its own comparison, AS joints against PsA joints, with the four paired patients reported as paired observations rather than as extra subjects. Each patient will carry its own detection floor, and any patient whose floor sits above 10⁻⁴ is **unevaluable rather than negative** — the rule that already set aside 7 of 32 AS subjects in blood.

## Reports

| | |
|---|---|
| [BULK_META.md](BULK_META.md) | The meta-analysis and the deconvolution: datasets, methods, every caveat. |
| [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md) | Why compartment beats modality, and the repertoire data inventory in full. |
| [NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) | Handoff: the patient-level plan, every dataset located, and the traps that cost time. |
| [REPORT.md](REPORT.md) | Earlier single-cell work — reanalysis of GSE277117 (Song et al., *Nat Commun* 16:5531, 2025), TNF-inhibitor response in 25 patients. |
| [NEXT_STEPS.md](NEXT_STEPS.md) | Project plan and the running list of traps. |

## Running it

```bash
python -m venv .venv && .venv/Scripts/pip install -r requirements.txt
```

The bulk meta-analysis runs locally in minutes. The single-cell reference needs [Modal](https://modal.com) — the object wants 64 GiB. Every script re-fetches its own data from GEO; nothing raw is committed.

### Never let a last-step failure discard a finished run

This cost us three times before it was fixed, most expensively when a stage-1 training run completed all 12 epochs and then died computing an optional diagnostic, throwing away the weights and the latents.

**`--detach` does not help with this, and it is worth being clear about why.** It detaches the app from the local client, so the run survives closing the terminal or losing the network. It does not checkpoint, does not persist partial results, and does not catch exceptions. The run that threw away 12 epochs *was* launched with `--detach`. Interrupted-connection and crashed-code are different failure modes and `--detach` only covers the first.

What actually works, in order of how much it saves:

1. **Persist the expensive artefact before anything optional.** Weights and derived arrays are what cost money; metrics are cheap to recompute. Save and `vol.commit()` immediately after training, then compute diagnostics.
2. **Wrap every optional diagnostic in `try/except`** and record the exception into the results dict instead of raising. A failed metric should appear as `"library_probe_error": "TypeError: ..."` in the output, not as a dead run.
3. **Smoke-test the whole path with one epoch before the real run.** The bug that killed stage 1 was reachable in 25 seconds and we found it after 5 minutes of GPU time. One epoch exercises every code path including the save and the diagnostics; the only thing it does not test is convergence.
4. **Checkpoint inside the training loop**, not only at the end, so a timeout or a preemption leaves something usable.
5. **Pin library versions in the image.** The failure was `LogisticRegression(multi_class=...)`, removed in recent scikit-learn — an API that existed when the code was written and not when it ran.

`scripts/vdj_expr_vae.py` is the worked example of all five.
