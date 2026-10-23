# Next experiment: patient-level disease prediction from TCR repertoires

Handoff document. Written at the end of the session that built and then
abandoned a cell-level model. Read the "What went wrong" section before
proposing anything.

## The actual goal

**Can a blood draw predict disease likelihood across multiple conditions?**

The inspiration is Faham et al. (*Arthritis Rheumatol* 2017): a CDR3 motif
present in 30–39% of ankylosing spondylitis patients and 0 of 227 controls.
That is a **patient-level** claim — does this person's repertoire contain a
particular clonotype — and it suggests the repertoire carries diagnostic
information when compared across people.

The target is a minimally invasive multi-disease diagnostic. Not a biological
model of T cell differentiation.

**Read Faham carefully before using it as motivation.** The headline compares a
mostly-HLA-B27-positive case group against mostly-B27-negative controls. The
comparison that isolates the motif from the HLA type — B27+ AS versus B27+
healthy — gives *p* = 0.049. The incremental value over a routine B27 test is
far smaller than the headline implies. The broad premise (repertoires differ by
disease) is better supported by Mal-ID and the coeliac classifier papers below
than by Faham's number.

## What went wrong last session, so it is not repeated

The session spent roughly $9 and most of its effort on a **cell-level** model:
given one T cell's receptor, predict that same cell's transcriptome. It was
adopted to work around label scarcity — millions of paired cells against ~666
labelled patients — and it quietly became the project.

It answers a different question and was never on the critical path. Four
specific process failures, in the order they happened:

1. **A loss function was diagnosed as the problem before a positive control
   existed.** The first design showed no real-vs-shuffled gap in its loss, which
   was attributed to a blunt objective (2,009 genes, 81% zeros, dominated by
   sequencing depth). An entire two-stage redesign followed. It scored 0.247
   against the original 0.244 — identical. The objective was never the
   constraint.
2. **The measurement unit was wrong for a long time.** Cells within one
   expanded clone share a receptor *and* a history, so they are not independent.
   Counting each clone once took the headline from 0.247 to 0.126. Half the
   signal was pseudo-replication.
3. **The targets were broken and unexamined.** Each "programme score" was a
   mean of ~9 normalised genes, which only measures one thing if those genes
   co-vary. CTLA4 sat in both the exhaustion and regulatory sets, making those
   scores correlate at 0.451. CD160 and BTLA loaded at 0.089 and 0.075 on the
   exhaustion axis. Coherence predicts recoverability: cytotoxicity (mean
   pairwise *r* = 0.397) recovered at 0.86; regulatory (0.150) at 0.70.
4. **No positive control for five rounds of redesign.** Every result was a
   negative measured against a shuffled baseline, which cannot distinguish "no
   signal" from "broken pipeline." Running the control first would have saved
   most of the session.

**The rule that follows: establish a positive control before interpreting any
negative, and before redesigning anything.**

## What the cell-level work established (and why it simplifies the real question)

| finding | consequence for the patient-level question |
|---|---|
| A receptor explains **1–2%** of its own cell's expression state (R² 0.009 val, 0.022 AS) | Do **not** route a patient-level predictor through predicted cell states. The base signal is too weak. |
| Plain features tie a 1.5M-parameter transformer: CD4/CD8 0.641 vs 0.658; MAIT 0.804 vs 0.802 | No GPU, no deep model needed. The usable receptor information is **V/J gene usage and amino-acid composition, not sequence order**. |
| Known biology does not transfer across cohorts: CD4/CD8 0.658 → **0.528**; MAIT 0.612 → **0.507** | Cross-cohort generalisation is the central obstacle. Any patient-level result must be validated **leave-one-cohort-out** or it means nothing. |

That third row is the most important thing in this document. A relationship as
solid as CD4-vs-CD8 lineage essentially vanished when the cohort changed.
Candidate causes, unexplored: the held-out study (GSE329527) is gut tissue with
only 5% CD8 cells, an extreme composition shift; V-gene *assignment* may differ
between cohorts because different Cell Ranger versions and reference sets call
V genes differently; or the subject/split handling.

## The plan

| step | what | cost |
|---|---|---|
| 1 | Pull **already-processed** clonotype tables from the AIRR Data Commons API — MS (93 repertoires), coeliac (10), SLE (10), COVID (3,356), plus controls. No FASTQ processing. | free |
| 2 | Patient-level features: V/J usage frequencies, CDR3 length distributions, clonality (e.g. Gini, top-clone share), public-clone sharing, generation-probability-weighted statistics | free, local |
| 3 | Classify disease, **leave-one-cohort-out**, against a frequency-matched null | free, local |
| 4 | **Validate against a published benchmark before trusting anything** | free |
| 5 | Only then pay to process the AS FASTQs | ~$5–15 |

### Step 4 is not optional

Best target: **CMV serostatus from repertoire** on Emerson et al.'s 666-subject
cohort — an established, replicated result. If the pipeline can call CMV status,
it works. That cohort lives in Adaptive's **immuneACCESS**, which is outside the
AIRR Data Commons and needs a free registration the user would have to do.

Fallback: PRJNA650391, "Classification of intestinal T cell receptor
repertoires using machine learning methods can identify patients with coeliac
disease regardless of dietary status" — a published repertoire classifier with a
number to beat. Also PRJNA678347 ("Immunosequencing of T Cell Repertoire
Reveals Celiac Disease State") and PRJNA1330789 (duodenal TCR-seq, coeliac plus
healthy controls).

### Prior art to read first

- **[Mal-ID](https://www.science.org/doi/10.1126/science.adp2407)** (*Science*,
  Boyd lab) — the closest thing to this project's goal. Scores individual
  BCR/TCR sequences with a language model, then aggregates to a patient-level
  call. Validated on lupus, type 1 diabetes, COVID and HIV. **Start here.**
- **BertTCR** — donor-level embeddings via a BERT encoder plus aggregation.
- An arXiv paper doing multimodal multiple-instance learning for autoimmune
  repertoire classification (SLE, RA).
- **OLGA / SONIA / soNNia** — generation-probability models. Pgen is a null
  model: observed frequency ÷ probability of arising by chance from V(D)J
  recombination. A clone abundant *despite* being unlikely to generate is
  evidence of antigen selection, found without knowing the antigen. soNNia is
  HLA-conditioned. Trained on non-productive (never-selected) sequences, which
  is why it is trustworthy. Good candidate for step 2 features.
- **ALICE** (Pogorelyy/Minervina, *PLOS Biology* 2019) — finds responding TCRs
  from a *single* snapshot by spotting sequences with more one-mismatch
  neighbours than the generation model predicts. Reanalysed 4 AS synovial
  repertoires from Komech.

## Data already available

### On the Modal volume `as-data` (app `as-experiment`)

| path | what |
|---|---|
| `/data/vdj/_dataset/` | 1,862,114 cells × 2,009 genes, paired receptor + transcriptome. `X.npz`, `cells.parquet`, `genes.txt`, `vocab.json` |
| `/data/vdj/_harmonised/` | 31,908 common genes, 2,009 selected, per-gene stats |
| `/data/vdj/_expr_latent/` | stage-1 expression latents at z = 16, 50, 100 plus checkpoints |
| `/data/vdj/_stage2/` | stage-2 checkpoints and all run results |
| `/data/tcr/GSE181364/` | TRUST4 output from bulk RNA-seq — 87–332 clonotypes/sample |

**This cell-level dataset is probably not what the next step needs**, but it
cost real money to assemble and it is the asset if anyone wants to revisit the
cell-level question with the lessons above applied.

### Repertoire data located but not downloaded

The first inventory (`scripts/airr_inventory.py`) queried the AIRR Data Commons
and returned **zero** AS/IBD/RA repertoires, from which it was concluded that a
multi-disease repertoire model could not be built from public data.
**That conclusion was wrong.** The ADC indexes *curated* repertoire studies; the
archives hold the rest, filed as ordinary sequencing. The systematic blind spot:
**10x 5′ libraries with V(D)J enrichment are indexed as "single-cell
RNA-seq."**

`scripts/repertoire_inventory.py` swept GEO, ENA and ArrayExpress instead:
2,715 candidate studies, 1,002 mentioning an inflammatory or autoimmune disease.

**AS / spondyloarthritis — six studies:**

| accession | design | size |
|---|---|---|
| **PRJNA393498** | Komech 2018 *Rheumatology* (`kex517`). 5′RACE TCRβ, blood + synovial fluid, TRBV9-sorted subsets. **Paper reports 25 AS + 108 healthy donors** — but only 78 runs are deposited here, so the control arm is probably a separately deposited Chudakov-lab healthy cohort (likely Britanova's aging cohort; unresolved, paper is paywalled) | 78 runs, 44.5 GB, 378M reads |
| PRJEB51315 = E-MTAB-11498 | Komech 2022 *Front Immunol*. 29 AS + 86 PsA, 64 SF + 51 PB. **No healthy controls** | 115 runs, 47.1 GB |
| GSE288581 | AS, mass cytometry + scRNA-seq + VDJ. "Resistance of CD8+ T cells to immune exhaustion" — AS CTLs express PD-1/TIGIT/LAG-3 yet retain granzyme/perforin. Held out as the test set last session | 30 samples, 27,375 paired cells |
| E-MTAB-10948 | AS scRNA-seq, PB + SF, 2 HLA-B27+ patients, >55,000 cells. No VDJ stated | — |
| PRJNA853355 = GSE207037 | single-cell TCR, DC-stimulated CD4 T cells | 1 run, 3.9 GB |
| GSE290921 | SpA, CD4 TRM Th17 driving IL-17A joint pathology | 6 samples |

Adjacent and on-target: **E-MTAB-9492** — 10x gene expression *plus* V(D)J from
blood, synovial fluid and synovial tissue in psoriatic arthritis. Only 5
patients, so a test set rather than training data.

**Cross-disease** (dedicated repertoire studies / with paired single-cell):
IBD 33/8 · coeliac 11/0 · type 1 diabetes 8/5 · rheumatoid arthritis 6/3 ·
AS 5/1 · MS 5/0. Standouts: GSE256256 (RA, seven-chain AIRR, n=112, controls),
GSE272431 (T1D, >2×10⁸ TCRβ from sorted CD4 subsets, patients and healthy
donors), PRJEB27352 (IBD in monozygotic twins — genetics-controlled).

Caveats: the `dedicated` flag is inferred from title text and "AMPLICON"
matches 16S microbiome studies, so a share of the IBD and psoriasis hits are gut
microbiome papers. GSE121082, logged as MS, is actually HTLV-I/HAM-TSP.

## Constraints and preferences

- **Get an explicit cost estimate approved before any non-trivial Modal job.**
  The user cancelled a run once on cost and has said so directly. Flag the
  number *before* launching, not alongside.
- **Do not create markdown files or commit without being asked.** This document
  and the README updates were explicitly requested.
- **Explain, do not narrate.** The user has asked more than once for less
  running commentary and clearer fundamentals. Define jargon on first use.
  Explain how to read a number, not just what it is: a loss value means nothing
  without a floor (null model) and a ceiling (information-saturated model), and
  the honest summary is the model's position between them.
- **The user pushes back well and is usually right.** They correctly rejected
  weighting marker genes by hand (it defeats the point of letting the model
  learn), correctly suspected the target definition rather than the model when a
  probe scored badly, and correctly identified that the project had lost its
  original question. Argue the substance, concede when wrong, and do not
  re-litigate.
- **`--detach` does not protect against code that fails at the last step.** It
  survives a lost terminal, not an exception. Persist the expensive artefact
  before any optional work, wrap diagnostics in try/except, smoke-test the whole
  path with one epoch, checkpoint to the volume rather than `/tmp`, and pin
  library versions. A run lost all 12 epochs to
  `LogisticRegression(multi_class=...)`, removed in recent scikit-learn.

## Traps that cost time last session

- **GEO deposits have no schema.** Files live under the sample directory
  (GSM-prefixed, bundled in a series tar) *or* the series directory
  (GSE-prefixed, loose). `filelist.txt` columns are
  `Archive/File | Name | Time | Size | Type` — **the name is the second field**;
  reading the first gave every file the name "File" and the paired-series count
  went 2 → 11 → 86 once fixed.
- **Sorted libraries silently mispair.** One study deposits `Bac17L_CD3pos` and
  `Bac17L_CD3neg`. Pairing the TCR file against the CD3-**negative** barcodes
  gave 0.8% overlap, which reads exactly like poor TCR recovery. Correct
  pairing: 83.3%.
- **One number rarely identifies a failure.** A low pairing rate means "wrong
  barcodes" or "not T-cell-sorted"; low marker enrichment means "bad join" or
  "every cell is a T cell." Three flags in a row were wrong because they tested
  a condition with two possible causes.
- **Antibody-capture features hide in gene matrices.** 8,831 CITE features were
  caught and dropped across the training set; GSE277117 had 163 confounded with
  the outcome. Filter on `feature_type` for `.h5` and the third column of
  `features.tsv`.
- **Sparse genes break dispersion-based gene selection.** `var/mean` rewards a
  gene detected in 0.3% of cells that spikes in a few; the first HVG run
  returned 2,000 genes whose *median* detection rate was 0.4%. A 1% detection
  floor fixed it and also raised cross-study HVG agreement from 32% to 50% —
  most of the apparent batch disagreement was sparse-gene noise.
- **TCR reads in bulk RNA-seq are a dead end.** TRUST4 recovered 146–498 TCRβ
  reads per sample, 0.0009% of reads, ~200 clonotypes out of a repertoire of
  10⁶. Reaching 10,000 clonotypes would need ~1.7 billion reads per sample,
  which costs more than a dedicated TCR kit and gives worse data.

## Things established earlier in the project, worth not rediscovering

- **A count matrix cannot hold a sequence.** GSE277117 could not find TRBV9
  because 10x 3′ chemistry discards V segments: TRBC2 in 33.4% of cells, TRBV9
  in 0.11%, V:C ratio 0.024. An assay limit, not a power limit.
- **Compartment beats modality.** HMP2, same metabolomics and same subjects:
  stool AUC 0.93 against plasma 0.74. AS has no draining compartment — the
  enthesis is small, poorly vascularised, and sheds little.
- **Discrete or foreign analytes make diagnostics; graded host-cell shifts do
  not.** Fetal cfDNA, dd-cfDNA, ctDNA, calprotectin and anti-CCP all work for
  the same reason. This is the argument for the repertoire being the right data
  type: a receptor is a discrete sequence, not a level.
- **Models cannot create information.** A cascade (predict B from A, then C from
  B̂) is capped by what A told you about C. Predicting TCR sequences from
  expression fails twice over: the CDR3 is produced by random nucleotide
  insertion that expression cannot encode, and even a perfect version could not
  support a Faham-style statistic, because any "shared clonotype" discovered
  among predicted sequences would be the model's output for the disease
  expression pattern.
- **Synthetic data can multiply inputs but not labels**, absent a validated
  oracle for the label. AlphaFold is an oracle for structure; OLGA is one for
  recombination; nothing is one for "does this repertoire come from an AS
  patient."

## Scripts

| script | what |
|---|---|
| `repertoire_inventory.py` | GEO + ENA + ArrayExpress sweep for dedicated repertoire data |
| `vdj_paired_inventory.py` | finds series with both V(D)J and GEX files via `filelist.txt` |
| `vdj_join_check.py` | measures the barcode join before trusting it |
| `vdj_manifest.py` | resolves per-library download URLs and sizes |
| `vdj_build.py` | Modal: downloads, joins, writes per-library arrays |
| `vdj_harmonise.py` | common gene axis, cross-study HVG selection by voting |
| `vdj_dataset.py` | consolidates into one training array |
| `vdj_expr_vae.py` | stage 1: expression → latent, with effective-rank and batch probes |
| `vdj_stage2.py` | stage 2: receptor → frozen latent, real and shuffled arms |
| `vdj_latent_probe.py` | round-trip and nonlinear probe diagnostics |
| `programme_coherence.py` | whether a marker gene set is internally coherent |
| `clonotype_eval.py` | re-evaluation at clonotype level, stratified by clone size |
| `positive_controls.py` | V gene, CD4/CD8, MAIT — model versus plain features |
| `airr_inventory.py` | the original ADC sweep whose conclusion was wrong |

Earlier phases: `geo_search.py`, `bulk_*.py`, `deconv_*.py`, `ml_loso.py`,
`clinical_benchmark.py`, `tcr_trust4.py`, `trbv_check.py`.
