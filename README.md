# Ankylosing spondylitis: what a T-cell repertoire can and cannot tell you

Independent reanalysis of AS repertoire and transcriptome data, done from scratch by a patient learning the methods. Three arms: patient-level repertoire analysis (the useful one), bulk transcriptomics, and a cell-level receptor→expression model.

## What is established

- A **seven-member TCRβ sequence family**, specified in the literature before we looked, is carried by **50.6% of AS patients and 4.0% of healthy controls** across three cohorts, three assays and three labs (MH odds ratio 23.3).
- The family is **CD8-restricted** and built on **TRBV9/TRBJ2-3**, which is what HLA-B27 class-I restriction predicts. Crystal structures show its members are **one shape** — 0.5–0.9 Å apart.
- **Nothing else is there.** Two open-ended sweeps, each with a working positive control, find no second marker and no distributed repertoire signal. The family is the signal.
- It does **not** separate AS from psoriatic arthritis. Given overlapping first-line treatment, the useful target condition is inflammatory spondyloarthritis versus no inflammatory disease.
- The best blood discrimination anyone has measured for this disease is **AUC 0.84**, and it is an antigen-specific tetramer assay, not a sequence test.

---

## The receptor family

[Hanson 2020](https://eprints.qut.edu.au/210039/) defines the AS-associated motif as **`CASSVG(V/I/L)(Y/F)STDTQYF`** — six sequences from two substitutions. [Yang 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC10511244/) solved crystal structures for three receptors, which adds a seventh sequence (`CASSVATYSTDTQYF`, their AS4.3). Faham 2017 and Komech 2018 reported `CASSVGLYSTDTQYF`, which is one member.

Seven sequences, fixed by other people's papers. Testing their union carries no selection leakage.

### Carriage, by cohort

| cohort | compartment | AS | sensitivity | control | specificity | one-sided Fisher p |
|---|---|---|---|---|---|---|
| Komech 2018 | unsorted blood, RNA amplicon | 14/25 | 56.0% | 0/10 | 100% | 0.0019 |
| Hanson 2020 | sorted CD8, Adaptive gDNA | 22/48 | 45.8% | 2/36 | 94.4% | 2.9×10⁻⁵ |
| GSE288581 | sorted CD8 memory, 10x | 3/4 | 75.0% | 0/4 | 100% | 0.071 |
| **pooled** | | **39/77 = 50.6%** (39.0–62.2) | | **48/50 = 96.0%** (86.3–99.5) | MH OR **23.3** |

Cohorts are not merged — each contributes its own 2×2, combined by Mantel-Haenszel. The single clone alone gives 12/48 and p = 0.016 in Hanson, so **the family nearly doubles sensitivity at identical specificity**.

Per sequence, in the two cohorts with control arms:

| cdr3 | source | Komech AS /25 | Hanson AS /48 | Hanson ctrl /36 | Emerson healthy /786 |
|---|---|---|---|---|---|
| `CASSVGLYSTDTQYF` | Hanson motif; Yang AS4.1 | 11 | 12 | 2 | 15 (1.9%) |
| `CASSVGLFSTDTQYF` | Hanson's headline clone; Yang AS4.2 | 6 | 12 | 0 | 20 (2.6%) |
| `CASSVGVYSTDTQYF` | Hanson motif | 4 | 2 | 0 | 11 (1.4%) |
| `CASSVATYSTDTQYF` | Yang AS4.3 | 4 | 3 | 0 | 7 (0.9%) |
| `CASSVGIFSTDTQYF` | Hanson motif | 0 | 3 | 0 | — |
| `CASSVGVFSTDTQYF` | Hanson motif | 2 | 1 | 0 | — |
| `CASSVGIYSTDTQYF` | Hanson motif | 0 | 0 | 0 | — |
| `CASSLGLFSTDTQYF` | *exploratory — found by us* | 7 | 4 | 1 | 31 (3.9%) |
| `CASSLGLYSTDTQYF` | *exploratory — found by us* | 4 | 8 | 0 | 33 (4.2%) |

`CASSVGIYSTDTQYF` has never been observed in any cohort. The two exploratory sequences are one edit from members but carry Leu at residue 95, outside Hanson's motif; we found them, so they stay out of the primary test. They add no new carriers in Komech after thinning.

Komech counts are presence at ≥2 reads after thinning to 800K junction reads ([`build_panel.py`](scripts/build_panel.py) definition); Hanson is ≥1 template at native depth. Hanson CD8 depth does not differ by group (median 66K vs 62K templates, p = 0.55).

### Hanson's labels are inferred, and checked

The 193-sample immuneACCESS export carries no diagnosis column. Groups come from the `A####` / `K####` subject prefixes, with `_4` and `_8` marking sorted CD4 and CD8. [`hanson_scan.py`](scripts/hanson_scan.py) checks that inference against the paper's own published counts rather than assuming it:

- Published: `CASSVGLFSTDTQYF` in 12 of 37 B27+ AS and 0 controls, CD8. Observed: **12 A-prefix, 0 of 36 K-prefix**.
- Published: such clones in CD4 "at 1 count in 2 samples". Observed: **2 A-prefix in CD4**, against 22 in CD8.

Both reproduce exactly. The authors' key should still be requested. Hanson states 19 of the 38 controls are B27+, and the 2 control carriers match the 2-of-19 B27+ carriers the paper reports, which would put specificity against B27+ controls near 89%.

### As a diagnostic

| prevalence | PPV | NPV |
|---|---|---|
| 0.5% (population) | 6.0% | 99.7% |
| 5% | 40.0% | 97.4% |
| 20% (back-pain clinic) | 76.0% | 88.6% |
| 50% | 92.7% | 66.0% |

A rule-in test: a positive result is informative, a negative one is not. Sensitivity is depth-limited — half the Komech carriers sit between 10⁻⁷ and 10⁻⁵ — and CD8 sorting concentrates the signal, so both deeper sequencing and sorting should raise it.

Against psoriatic arthritis it does not separate (PRJEB51315, the only cohort with the differential):

| comparison | AS | PsA | specificity vs PsA | p |
|---|---|---|---|---|
| family union, blood | 5/7 | 5/16 | 68.8% | 0.092 |
| family union, synovial fluid | 5/7 | 6/17 | 64.7% | 0.12 |
| `CASSVGLYSTDTQYF` alone, blood | 5/7 | 4/16 | 75.0% | 0.052 |

The union costs specificity against PsA, because the extra members are carried by PsA patients too. How much that matters depends on the target condition. AS and PsA share their first-line biologics and either diagnosis routes a patient to a rheumatologist, while the failure that actually occurs is years of inflammatory back pain treated as mechanical — a 7–9 year typical diagnostic delay. So the useful target is **inflammatory spondyloarthritis versus no inflammatory disease**, where a positive result means "refer". Two places the AS/PsA distinction still costs something: IL-23 inhibitors work in PsA and failed in axial SpA, and methotrexate helps peripheral PsA but not axial disease. The comparison never yet possible is AS versus *mechanical* back pain, for which no public dataset exists.

### Why it is CD8, and TRBV9

Across all 64 PRJEB51315 synovial libraries, on a within-subset denominator:

| group | subset | target frequency | junction reads |
|---|---|---|---|
| **AS** | **CD8** | **1.00×10⁻³** | 4,378,333 |
| AS | unsorted | 2.69×10⁻⁴ | 3,569,466 |
| AS | CD4 | 2.88×10⁻⁵ | 2,533,848 |
| PsA | CD8 | 8.78×10⁻⁵ | 10,914,688 |
| PsA | CD4 | 1.56×10⁻⁶ | 6,397,279 |

AS joint CD8 sits **35× above AS joint CD4** and 11× above PsA joint CD8. HLA-B27 is class I, read by CD8 cells. In Hanson's sorted blood the same pattern holds: 22 AS carriers in CD8 against 2 in CD4.

[`vgene_pilot.py`](scripts/vgene_pilot.py) assigns V genes by comparing the 60 nt of V sequence 5′ of the conserved cysteine against every germline TRBV allele (OLGA's IMGT reference; no aligner needed, as V segments are not spliced here). **`CASSVGLYSTDTQYF` is TRBV9** in all three AS patients tested — 71 reads, none assigned elsewhere, 35 of 39 at distance 0 with a median margin of 3. Three practical points:

- **TRBV9 is not itself a marker.** It is 3.61% of V-assigned reads in an AS patient and 2.53% in a healthy control. The rare thing is the pair, not either part.
- **Assignment needs a margin rule, not a distance cutoff.** TRBV9 sits 3 mismatches from TRBV13 and TRBV7-4. Requiring a ≥2 margin over the runner-up refuses a third of reads and shifts the usage profile (TRBV20-1 moves 18.2% → 29.8%). Irrelevant for naming one clone; essential for comparing profiles.
- **Read length is a group confound.** 24 of 27 AS libraries are 2×100 bp against 8 of 12 healthy, so `--truncate 100` is mandatory for cross-group usage comparisons.

### Structurally, the family is one shape

Yang deposited seven TCR–peptide–HLA-B27 complexes; three are family members holding the *same* peptide (PRPF3). [`cdr3_structure_compare.py`](scripts/cdr3_structure_compare.py) superposes each on the HLA platform — not on the loops, which would flatter the comparison — and measures the CDR3β backbone:

| CDR3β RMSD (Å) | AS3.1 + PRPF3 | AS4.2 + PRPF3 | AS4.3 + PRPF3 | AS4.2 + YeiH |
|---|---|---|---|---|
| **AS3.1** `CASSVGLYSTDTQYF` | — | 0.51 | 0.94 | 1.94 |
| **AS4.2** `CASSVGLFSTDTQYF` | 0.51 | — | 0.62 | 1.52 |
| **AS4.3** `CASSVATYSTDTQYF` | 0.94 | 0.62 | — | 1.10 |

Two independent structures of the same protein typically differ by 0.3–0.5 Å, so the family members are the same shape. Swapping the *antigen* under one unchanged receptor (last column) moves the loop further than any sequence difference within the family does.

Recognition is β-centric: CDR3β residue 98 — the Y or F that defines the family — is the peptide contact, CDR1α from TRAV21 touches peptide position 4, and CDR3α barely participates. [`figures/as-receptor-anatomy.html`](figures/as-receptor-anatomy.html) draws this.

The TRAV21 α chain is half the functional unit and is absent from every bulk TCRβ dataset here. [Paley 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11343591/) adds that **TRBV9 is characteristic but not required** — 5 of 20 motif-bearing receptors used TRBV5-4, TRBV5-5 or TRBV7-3 — so a TRBV9-conditioned search would discard a quarter of the pathogenic repertoire.

---

## Is anything else there? No, and the sweeps are now informative

Two open-ended sweeps, each reporting a positive control.

### Hanson, 2.2M clonotypes, 48 AS vs 36 controls

[`hanson_sweep.py`](scripts/hanson_sweep.py) tests every clonotype in the sorted-CD8 compartment against 200 label permutations, with the CD4 fraction swept identically as a negative control.

| α | observed | null mean | null 95th | empirical p |
|---|---|---|---|---|
| 0.05 | 655 | 583 | 1,071 | 0.34 |
| 0.01 | 77 | 72 | 155 | 0.37 |
| 0.001 | 3 | 2.5 | 10 | 0.34 |

Null at every threshold. CD4 lands *below* its null (624 against 1,033). And the positive control is as clean as it gets: **`CASSVGLFSTDTQYF` is the top hit of the entire 2.2-million-clonotype sweep**, 12/48 at p = 0.0006, with `CASSVGLYSTDTQYF` also recovered at 12/48.

A thinned variant (40K templates, which costs 14 of 48 AS patients) shows 7 clonotypes at α = 0.001 against 1.9 expected, p = 0.015 — but CD4 shows a comparable excess at α = 0.01 (p = 0.010), so that tail is not CD8-specific and is not read as signal. Thinning is unnecessary here because group depths match; both settings are reported.

### Komech's untestable candidates, tested in Hanson

The Komech sweep (25 AS, 10 controls) was power-limited by construction: with 10 controls a clonotype needs ≥12 of 25 AS carriers and 0 controls before α = 0.01 is reachable. The funnel:

| stage | clonotypes |
|---|---|
| distinct clonotypes observed | 5,329,665 |
| after thinning to 800K reads/subject | 2,678,691 |
| in ≥1 AS patient, 0 controls | 1,685,436 |
| in ≥6 AS patients, 0 controls | 1,081 |
| **in ≥12 AS patients, 0 controls — the only testable zone** | **6** |

Hanson's 36 controls drop that threshold to 8 of 48 (17%), so the 1,075 clonotypes in the 6–11 band became testable. [`komech_candidates_in_hanson.py`](scripts/komech_candidates_in_hanson.py) scores them in Hanson — discovery and evaluation in different cohorts, on different platforms, with 1,075 random decoys and the CD4 fraction alongside.

| set | n | median AS carriers | p ≤ 0.01 | survive FDR ≤ 0.10 |
|---|---|---|---|---|
| Komech candidates | 1,073 | 2 | 1 | **0** |
| decoys | 1,075 | 0 | 0 | 0 |
| **motif family** | 7 | 3 | 1 | **2** |

As a burden score — panel from Komech, counted per Hanson subject, which is the leak-free form of the score retracted below:

| set | fraction | AS median (IQR) | control median (IQR) | AUC | p |
|---|---|---|---|---|---|
| Komech candidates | CD8 | 62 (30–85) | 62 (41–75) | **0.492** | 0.55 |
| decoys | CD8 | 1 (0–2) | 2 (0–2) | 0.457 | 0.76 |
| **motif family** | **CD8** | 0 (0–1) | 0 (0–0) | **0.707** | **2.4×10⁻⁵** |
| motif family | CD4 | 0 (0–0) | 0 (0–0) | 0.507 | 0.36 |

0.492 is chance. The family reaches 0.707 on the same subjects through the same code and goes null in CD4.

Per sequence the candidates look better than decoys (median 2 carriers against 0), but that is prevalence, not signal — candidates were selected for ≥6 Komech carriers, and they are equally common in Hanson's *controls*. The within-cohort burden comparison settles it.

**Conclusion.** There is no second marker set and no distributed burden signal. The family is the signal, which supports an antigen-driven clonal model over a diffuse repertoire shift.

One sequence stays on a watchlist rather than in a result: `CASSFRASYEQYF` — 7/25 Komech AS, 8/48 Hanson AS, 0/36 controls, 0 in CD4, p = 0.0087, which does not survive correction for 1,073 tests.

### Two control arms, not one

An external cohort kills public clonotypes; the study's own controls kill batch artefacts.

The 14 top within-study Komech hits are all public — present in 35–85% of Emerson's 786 healthy repertoires (`CASSLAGYGYTF` in 664, `CASSFGTGGYEQYF` in 559). Ten controls cannot exclude a clonotype two-thirds of humanity carries. Our target was the sole exception at 1.9%.

Screening the other way, five of the nine Komech candidates rare in Emerson are a deposit-specific ladder — `CRVSGL[WGCR]AKNIQYF` plus `CSWGTHDALCSNQPQHF` — absent from all 786 external repertoires but carried by 30–50% of Komech's *own* controls. Against Emerson alone they reach p = 10⁻¹¹; against Komech's controls, p = 0.71.

### The burden score, retracted

[`as_burden_panel_167.tsv`](data/panel/as_burden_panel_167.tsv) keeps a CDR3 appearing in ≥5 of 25 Komech AS patients and ≤5% of Emerson's 786. It reached AUC 0.962 within Komech and 0.991 against an independent healthy cohort. **Both are withdrawn**: the candidate pool was defined across all 25 AS patients including the held-out one, so only the second-stage filter was fold-restricted.

Rebuilt with the pool from 24 training patients only:

| panel rule | median size | leak-free LOO AUC |
|---|---|---|
| ≥4 carriers | 907 | 0.528 |
| **≥5 carriers** | 150 | **0.672** |
| ≥6 carriers | 37 | 0.620 |

Permutation-tested with the carrier grid inside the null: observed 0.672, null mean 0.561, 95th percentile 0.744, **p = 0.17**. The Hanson burden test above is the independent confirmation.

One candidate survives from this work with a Bonferroni-corrected p but a shared select-and-test carrier count: `CASSFGTGGINEQYF`, 7/25 Komech AS and 5/796 healthy (0.6%), reaching 2.2×10⁻⁴ in P15. Not a PCR error variant of its abundant 2-edit neighbour `CASSFGTGGYEQYF` — it outnumbers the supposed parent 10× and 89× in the two patients carrying both. It partially replicates in the joint: 2 of 7 AS synovial samples above 10⁻⁴ against 1 of 17 PsA joints, Mann-Whitney p = 0.061.

---

## Representations that failed, and why

| representation | result |
|---|---|
| exact CDR3 strings | 5.3M features for 36 patients; produced a leaked AUC of 0.96 |
| edit distance ≤1 | saturates — 74% of healthy repertoires contain something within 1 edit, taking p from 0.0033 to 0.58 |
| ESM-2 cosine similarity | returns anagrams as nearest neighbours; mean-pooling discards order |
| V gene alone | 59 features; everyone carries every V gene |
| functional distance (peptide recognition profiles) | real, measured, and does not transfer to disease — below |
| **V gene + similar CDR3 families** | the biologically motivated unit; **untested** |

The first two share a cause: both ignore the V gene. "Within 1 edit" sweeps in receptors on different V segments that engage different MHC molecules.

### Functional distance has wet-lab ground truth and still does not transfer

[Wang, Yeh & Garcia 2026](https://doi.org/10.1038/s41587-026-03128-x) profiled 21 HLA-B\*27:05-restricted TCRs against a ~10⁹-variant peptide library — the Faham `S··YSTDTQYF` family, measured. [`prp_divergence.py`](scripts/prp_divergence.py) reproduces their result that **motif Jensen-Shannon distance is uncorrelated with CDR3β edit distance** (ρ = +0.113, p = 0.10), with their engineered 19.2 neighbourhood recovered as a positive control. The clone's nearest functional relatives sit *three* edits away, so radius-1 widening searched the wrong neighbourhood rather than too wide a one.

Scanning all 16 panel sequences against 31 AS, 10 healthy and 16 PsA subjects ([`prp_candidates_scan.py`](scripts/prp_candidates_scan.py)) inverts the expectation: AS carriage tracks **edit** distance (ρ = −0.644, p = 0.007), not functional distance (ρ = −0.028, p = 0.92). The clearest case is `CASSLGLYSTMEQYF`, the one panel sequence on a different J region (TRBJ2-5 against TRBJ2-3) — absent from all 57 subjects. It may well bind the same peptides; it is a different rearrangement, so it is not generated. **Functional similarity and disease association are different questions**: a clone must also be generatable and selected.

Read the negative weakly — functional distance exists for 16 sequences, 9 with zero AS carriers. Two by-products: their released checkpoints **cannot** score novel CDR3βs (single-TCR models held β constant while fine-tuning; the joint model learned β as a near-pure intercept, r = 0.9996 between β logit vectors, so [`prp_score.py`](scripts/prp_score.py) records a failed control), and the cached profiles sit inside 3.2 GB archives that are 99.96% model weights, so range-fetching one zip member cost 32 MB instead of 70.6 GB.

---

## The best blood discrimination measured: AUC 0.84, not a sequence test

[Paley 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11343591/) stained blood with an HLA-B27 tetramer loaded with **YeiH** — the bacterial peptide Yang showed this family recognises — and counted the T cells that bind it. They publish group comparisons but no discrimination statistic; [`paley_tetramer_roc.py`](scripts/paley_tetramer_roc.py) computes one from their open figure source data:

| measurement | AUC (95% CI) | p |
|---|---|---|
| **YeiH-specific, % of CD8** (31 patients vs **10 B27+ healthy**) | **0.837 (0.673–0.961)** | 0.0008 |
| YeiH-specific, per 10⁶ PBMC | 0.829 (0.645–0.974) | 0.0010 |
| influenza NP, % of CD8 — *specificity control* | 0.642 (0.473–0.805) | 0.093 |
| influenza NP, per 10⁶ PBMC — *same control* | 0.719 (0.545–0.871) | 0.020 |

Holding separately: B27+ uveitis 0.841, B27+ axSpA 0.850. This beats everything here — the bulk NK score (0.77), the composition model (0.787), the family carrier count (0.780) — and it measures antigen specificity rather than sequence identity, so depth and exact matching stop mattering.

Two caveats. The influenza control is not clean (p = 0.020 on one scale), so part of the separation may be a general difference in these patients' CD8 compartment; the per-person YeiH/NP ratio would settle it. And 10 controls cannot fix an operating point — the best Youden cut gives 97% sensitivity at 60% specificity.

The mechanism behind the reagent: Yang's yeast-display selection found YeiH (`LRVMMLAPF`, conserved across *E. coli*, *Salmonella*, *Shigella*, *Klebsiella*) and the self-peptides GPER1, PRPF3 and RNASEH2B. GPER1 activates these T cells when presented by disease-associated B\*27:05 but **not** by protective B\*27:09, while YeiH is presented by both — so the subtype genetics track the autoantigen, not the microbe. Paley adds that the cells carry gut-homing markers (α4β7, CCR6, CD103, CD161), suggesting intestinal priming.

---

## The earlier arms

### Bulk transcriptomics: NK cells are quiet, not absent

Five cohorts, 97 AS patients vs 70 healthy controls, random-effects meta-analysis, then deconvolution against a 312,000-cell reference.

| | pooled *g* | 95% CI |
|---|---|---|
| NK cytotoxic cargo (granzymes, perforin, granulysin) | **−0.98** | −1.33 to −0.63 |
| …after adjusting for NK cell proportion | **−0.80** | −1.17 to −0.44 |
| NK cell proportion | −0.54 | −0.91 to −0.17 |

82% of the cytotoxic deficit survives adjustment for cell number. The two numbers use **disjoint gene sets** — proportion from NK identity genes with every effector gene deleted, outcome from the effector genes — or the analysis is circular. A [meta-analysis of 95 flow studies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8291033/) (4,020 patients) found no difference in NK percentage; flow is blind to transcription, deconvolution mistakes a quiet cell for an absent one.

Not a diagnostic: AUC 0.77, and in a back-pain clinic a positive result means a 17.7% chance of AS (CRP 15.6%, HLA-B27 37.2%). Deconvolution was validated on held-out pseudobulk from 13 unseen patients: NK recovered at r = 0.83, **CD8 T failed at 0.57 and was discarded** despite p = 0.029, because cytotoxic T cells share GNLY, NKG7 and PRF1 with NK cells.

Unfollowed lead: platelet and integrin genes are the strongest upward signal (10 of 11 markers up, mean *g* = +0.44), matching a [2025 plasma proteomics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12957791/) whose top AS biomarkers were integrin-pathway proteins.

### Cell-level model: a receptor explains 1–2% of its own cell's state

1,862,114 cells across 291 libraries with both a productive TRB and a transcriptome; 91.3% singleton clonotypes. Stage 1 compresses 2,009 genes to a latent from expression alone; stage 2 predicts that frozen latent from CDR3 + V/J. Every measurement runs twice, once with receptors permuted within subject — the gap is the signal.

| | real | shuffled |
|---|---|---|
| exhaustion, per cell | +0.247 | −0.003 |
| exhaustion, **per clonotype** | **+0.126** | +0.010 |
| exhaustion, held-out cohort | +0.012 | +0.005 |

Half the per-cell signal was pseudo-replication. Positive controls recover known links at literature magnitudes (CD4-vs-CD8 0.658, MAIT 0.802). Two findings outlast the model:

- **Plain features match the transformer.** Logistic regression on V/J one-hot, CDR3 length and amino-acid composition ties a 1.5M-parameter transformer to within 0.02 AUC. The information is gene usage and composition, not sequence order.
- **Nothing transfers across cohorts, not even known biology.** CD4-vs-CD8 falls 0.658 → 0.528 on a held-out study; MAIT 0.612 → 0.507.

---

## Data-handling rules, earned the hard way

**Numerator and denominator in the same code path.** A 34× frequency inflation survived review because the motif was counted on both strands and its denominator on one: in these 5′RACE libraries read 2 carries the CDR3 reverse-complemented, so `CASS` was counted in 0.72% of reads where the truth is 24.3%. It cost the headline claim "AS 9/25 expanded ≥10⁻⁴ versus healthy 0/785" — on the honest denominator no AS patient in either cohort reaches 10⁻⁴. Presence numbers were never affected.

**Thin to the shallowest library, not a round number.** The first sweep thinned to 30,000 junction reads from libraries holding 5 million, leaving only clones above ~10⁻⁴ visible — it could not have found its own positive control. Thin only when depth differs by group; where it does not, thinning costs subjects to remove a confound that is not there.

**Separate selection from evaluation at every stage.** A fold-restricted final step makes a pipeline look cross-validated while an upstream filter has already seen the labels; that alone turned p = 0.17 into AUC 0.96. A permutation test catches it only if the leaked stage sits inside the shuffle.

**Two control arms.** External kills public clonotypes, internal kills batch artefacts, neither does the other's job.

**Every sweep reports a positive control.** A screen that cannot recover a known signal is uninformative, not negative.

**Permutation nulls are cheap, so they are non-negotiable.** With group sizes fixed, the one-sided Fisher p depends only on the carrier count among cases and the total carrier count, so one precomputed table serves every shuffle, and each shuffle is a sparse matvec. Milliseconds instead of hours.

**Persist the artefact.** The 167-sequence panel existed only in memory for several sessions, making every number from it unreproducible.

**Count distinct accessions before aggregating.** A recursive glob counted 153 runs as 306 because scan JSONs were duplicated across directories; the detection floors derived from them were 2× optimistic.

**Cross-deposit comparisons need calibration.** Cohort membership is predictable from repertoires at AUC 0.994 and HLA-A\*02 at 0.953, while CMV serostatus transfers across cohorts at 0.851 and age+sex at 0.584. Pooling subjects across deposits is therefore avoided in favour of within-study estimates combined afterwards.

Two assay limits worth recording: **bulk RNA-seq cannot test this** (GSE221786 gives a detection floor of 3.6–3.8×10⁻⁴, about 4× too coarse, needing ~190M reads/sample), and an early pilot's "0.50% of reads are `CASS`" was chance — `CASS` is 12 nt, matching ~1 position in 58,000, and 139 positions × 2 strands predicts 0.476%.

---

## Data inventory

**The archives were never the limit; email was.** The largest AS repertoire dataset here — Hanson 2020, 193 samples, 47 AS and 38 B27-matched controls with CD4 and CD8 sorted — is in no public archive; it sits on immuneACCESS under the authors' account and was obtained by asking. Faham 2017's 191 B27+ AS patients are absent for a different reason: Sequenta's LymphoSIGHT platform never migrated after the 2015 Adaptive merger. For AS, assume investigator-held data exceeds public data and write to the corresponding author first.

Also systematically missed: **10x 5′ libraries with V(D)J enrichment are indexed as "single-cell RNA-seq"** and no repertoire inventory surfaces them. Of 450 series checked, 86 have both VDJ and GEX — 6,008 samples.

Public AS archives are otherwise exhausted: `PRJNA853355` is a single CD4 sample, `E-MTAB-11498` mirrors `PRJEB51315`, `GSE207037`/`GSE290921` are CD4-only, `E-MTAB-10948` has no VDJ. Worth pursuing, in order:

1. **Steel 2019** (immuneACCESS `/pub/steel-2019-ar`) — IL-17A+CD8+ cells in PsA, Adaptive platform, so assay-matched to Hanson. The PsA arm rests on 16 subjects from one deposit; this is the cheapest fix and needs only a download.
2. **HLA-B27-associated acute anterior uveitis** — `GSE274471` (88 samples, aqueous humour + blood scTCR), `GSE232802`, `PRJNA973988`. If the family is B27-driven rather than AS-driven it should appear in B27+ uveitis too.
3. **GSE213095** — HLA-B27+ arthritis with paired GEX and VDJ, 46 samples, inventoried but never opened.
4. **HLA-B27 status for the Komech patients**, absent from all metadata.
5. **Emerson per-subject repertoires** are no longer on disk, so the family's specificity against 786 unselected healthy donors is bounded (2.5–6.7%) rather than measured. Re-exportable through immuneACCESS.

---

## Reports

| | |
|---|---|
| [BULK_META.md](BULK_META.md) | The meta-analysis and deconvolution: datasets, methods, caveats. |
| [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md) | Why compartment beats modality, and the repertoire inventory in full. |
| [NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) | Handoff: the patient-level plan and the traps that cost time. |
| [REPORT.md](REPORT.md) | Earlier single-cell work — GSE277117, TNF-inhibitor response in 25 patients. |
| [NEXT_STEPS.md](NEXT_STEPS.md) | Project plan and the running list of traps. |
| [important-papers.md](important-papers.md) | Faham, Hanson, Yang, Paley, Mal-ID and the PRP paper, with what each means here. |

## Running it

```bash
python -m venv .venv && .venv/Scripts/pip install -r requirements.txt
```

The bulk meta-analysis runs locally in minutes. Repertoire extraction and the single-cell reference need [Modal](https://modal.com). Every script re-fetches its own data; nothing raw is committed.

```bash
# extract CDR3s from raw reads (Modal, caches per run on the volume)
modal run --detach scripts/komech2018.py::cdr3_all
modal run --detach scripts/as_trust4.py::cdr3_all

# the family, in Komech and in Hanson
python scripts/motif_family_test.py --tables <komech dir>
python scripts/hanson_scan.py --export "<hanson export dir>"

# open-ended sweeps, each with a positive control
python scripts/hanson_sweep.py --export "<hanson export dir>" --min-carriers 3 --depth 0
python scripts/komech_candidates_in_hanson.py --tables <komech dir> --export "<hanson export dir>"

# panels, V genes (local, no Modal)
python scripts/build_panel.py --tables <dir>
python scripts/vgene_pilot.py --runs SRR5812652 SRR5812641 --truncate 100
```

| script | what it produced |
|---|---|
| [`motif_family_test.py`](scripts/motif_family_test.py) | the family union in Komech; **defines the seven pre-specified sequences, which the others import** |
| [`hanson_scan.py`](scripts/hanson_scan.py) | the A/K label inference and its check against Hanson's published counts |
| [`hanson_stats.py`](scripts/hanson_stats.py) | the Hanson Fisher tests — 22/48 vs 2/36 and that table |
| [`hanson_sweep.py`](scripts/hanson_sweep.py) | the 2.2M-clonotype sweep, its permutation null and the CD4 control |
| [`komech_candidates_in_hanson.py`](scripts/komech_candidates_in_hanson.py) | the 1,075-candidate replication and its burden score |
| [`family_cohorts.py`](scripts/family_cohorts.py) | family carriage in PRJEB51315 (AS vs PsA) and GSE288581 |
| [`family_meta.py`](scripts/family_meta.py) | Mantel-Haenszel pooling, the PPV table, the control-arm power floors |
| [`cdr3_structure_compare.py`](scripts/cdr3_structure_compare.py) | the CDR3β RMSD matrix from Yang's PDB entries |
| [`paley_tetramer_roc.py`](scripts/paley_tetramer_roc.py) | the AUC 0.837 Paley never published |
| [`bulk_score.py`](scripts/bulk_score.py) | the NK score's AUC 0.77 and the back-pain-clinic PPVs |
| [`rnaseq_pilot.py`](scripts/rnaseq_pilot.py) | the bulk RNA-seq ruling-out |
| [`esm_embed.py`](scripts/esm_embed.py) | the ESM-2 negative result |
| [`as_clonotypes.py`](scripts/as_clonotypes.py) | `data/as_vdj/clonotypes.parquet`, the GSE288581 numbers |

### Never let a last-step failure discard a finished run

A training run once completed all 12 epochs then died computing an optional diagnostic, throwing away the weights. `--detach` does not help — it survives a closed terminal, not a worker crash. What works, in order of value:

1. **Persist the expensive artefact before anything optional.** Save and `vol.commit()` immediately after training, then compute diagnostics.
2. **Wrap every optional diagnostic in try/except** and record the exception into the results dict.
3. **Cache per unit of work.** Each extraction function checks for its own output first; this turned a stalled 30-of-39 run into a 9-library job.
4. **Use `.map(..., return_exceptions=True)`.** A container timeout or OOM is killed outside the interpreter, so try/except inside the worker cannot catch it.
5. **Smoke-test the whole path with one epoch.** The bug that killed stage 1 was reachable in 25 seconds and found after 5 minutes of GPU time.
6. **Pin library versions in the image.**

[`vdj_expr_vae.py`](scripts/vdj_expr_vae.py) is the worked example. Modal also preempts containers, which costs wall time, not correctness.
