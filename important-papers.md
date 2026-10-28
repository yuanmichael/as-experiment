# Important papers

Papers this project leans on, with the findings worth citing and the caveats worth
remembering. Two paragraphs each: what they did, and what it means for us. Added
to as we go.

---

## Discovery of T Cell Receptor β Motifs Specific to HLA–B27–Positive Ankylosing Spondylitis by Deep Repertoire Sequence Analysis

*Faham et al., Arthritis & Rheumatology 2017; 69(4):774–84.*
[10.1002/art.40028](https://acrjournals.onlinelibrary.wiley.com/doi/10.1002/art.40028)

Sequenced TCRβ repertoires from **191 HLA-B27-positive AS patients, 43 B27-negative
AS patients and 227 controls**, yielding over 77 million clonotypes, and found CDR3
motifs enriched in B27-positive AS. The headline figure is a motif present in
**30–39% of AS patients**, and an independent cohort supported at least some motifs
tracking AS rather than B27 status alone. The `SVGLYSTDTQ` motif from this paper is
carried by `CASSVGLYSTDTQYF`, the clonotype this project's sequence work is built
around.

**Read the headline carefully.** It compares a mostly-B27-positive case group
against mostly-B27-negative controls, so it partly measures HLA type rather than
disease. The comparison that isolates the motif — B27+ AS versus B27+ healthy —
gives *p* = 0.049, so the incremental value over a routine B27 test is far smaller
than 30–39% suggests. Separately, this is **the dataset that would resolve our main
power limit**: 234 AS patients with 227 controls on a single assay, against our 25.
A comprehensive immuneACCESS search did not surface it, and asking the authors is
probably the highest-value outreach available. It is **not** on immuneACCESS because
Faham and most co-authors were at **Sequenta** (platform: LymphoSIGHT), which merged
with Adaptive in 2015; the data never moved to the public catalogue. Corresponding
author Robert Inman, `robert.inman@uhn.on.ca`.

---

## Altered Repertoire Diversity and Disease-Associated Clonal Expansions Revealed by T Cell Receptor Immunosequencing in Ankylosing Spondylitis Patients

*Hanson et al., Arthritis & Rheumatology 2020; 72(8):1289–1302.*
[10.1002/art.41252](https://doi.org/10.1002/art.41252) ·
[free full text at QUT](https://eprints.qut.edu.au/210039/)

Profiled the peripheral CD4 and CD8 TCRβ repertoires of **47 AS patients against 38
HLA-B27-*matched* healthy controls** — the only properly B27-matched control arm in
the AS repertoire literature — with CD4 and CD8 sorted separately, on Adaptive hsTCRB
at "Survey" resolution (~30,000 T-cell genomes per reaction). A subset of 19 patients
and 14 controls was additionally stimulated in vitro. At the aggregate level, AS
repertoires are **more** diverse than B27-matched controls by Normalized Shannon
Diversity Entropy (CD4 *p* = 7.8×10⁻⁶, CD8 *p* = 9.3×10⁻⁴), which the authors
attribute to AS patients holding *fewer* large clonal expansions; fewer patient T
cells produced IFN-γ (CD8 *p* = 0.03) or TNF (CD4 *p* = 0.01, CD8 *p* = 0.002) on
stimulation; and B27+ AS patients showed expanded EBV- (*p* = 0.03) and
CMV-specific (*p* = 0.02) clonotypes.

The clone-level results are the substance. Ten public CD8 CDR3s separate B27+ AS from
B27+ healthy controls at an sPLS-DA classification error rate of **0.29 ± 0.05**;
three of them (`CASSLGRAYEQYF`, `CASSSGTGGNQPQHF`, `CASSLSGGNTEAFF`) preferentially
use **TRBV27**, and only one (`CASSLGSTYEQYF`) was also associated with B27− disease.
Their single strongest clone is **`CASSVGLFSTDTQYF`**, detected in **12 of 37 B27+ AS
patients and absent from every B27− patient and every control** (*p* = 1.6×10⁻⁵).
They define the family as **`CASSVG(V/I/L)(Y/F)STDTQYF`** and report the union of its
two commonest members, `CASSVGLFSTDTQYF` + `CASSVGLYSTDTQYF`, in **18/37 (49%) of
B27+ AS versus 2/19 (11%) of B27+ healthy** (OR 8.05, *p* = 0.007). The motif's share
of the CD8 repertoire is 0.006% in B27+ AS vs 0.0003% in B27+ healthy (*p* = 0.002),
and the 24 clonotypes within one amino acid of the F variant sum to 0.008% vs 0.002%
(*p* = 0.003). Such clones were essentially absent from CD4 (1 count in 2 samples).
The association was independent of age, sex, TNF-antagonist use and BASDAI. Two
limits: Survey depth is ~14× shallower than Komech, so absence here is weak evidence;
and the diversity direction *contradicts* the 2018 Chinese pilot that reported
reduced diversity in AS, so repertoire-diversity features are unstable across
cohorts. Immunosequencing data are "available from the corresponding author upon
request" — Matthew Brown, `matt.brown@kcl.ac.uk` (KCL / Genomics England);
Supplementary Tables 7–9 hold the 10 public clones and the 24 homologous variants.

**Against this repo.** Their headline clone is not ours, and ours is the family's
second member rather than its only one. `CASSVGLFSTDTQYF`, `CASSVGVYSTDTQYF`,
`CASSLGLFSTDTQYF` and `CASSLGLYSTDTQYF` are all present in our panel files at 4–7
Komech AS carriers and 0 Komech controls (see the README's motif-family table), so
this is a replication we already held. Their B27-matched control arm is the thing we
cannot build from public data: our best substitute is 45 B27+ Emerson donors.

---

## Autoimmunity-associated T cell receptors recognize HLA-B*27-bound peptides

*Yang et al., Nature 2022; 612:771–777.*
[10.1038/s41586-022-05501-7](https://doi.org/10.1038/s41586-022-05501-7) ·
[free at PMC10511244](https://pmc.ncbi.nlm.nih.gov/articles/PMC10511244/)

Took the public AS/uveitis TCR motif into the wet lab and solved its antigens.
Starting from the TRBV9 + `(Y/F)STDTQ` + TRBJ2-3 CDR3β motif — clonally expanded
10–100× in AS joints and uveitis eyes, and **uniformly paired with TRAV21**
α-chains — they
displayed a stabilised single-chain HLA-B\*27:05 trimer on yeast with ~10⁹-member
9-mer and 10-mer peptide libraries, ran four rounds of TCR-driven selection, deep-
sequenced the survivors into a position-wise binding motif, then searched the human,
*Chlamydia*, *Salmonella*, *Shigella*, *Klebsiella* and *Yersinia* proteomes for real
peptides matching it. Hits were tested by CD69 upregulation on TCR-transduced SKW-3
cells against HLA-B\*27:05⁺ K562 cells, with affinities by SPR and **seven crystal
structures** (PDB 7N2N/O/P/Q/R/S, 8CX4). The microbial hit is **YeiH** (`LRVMMLAPF`),
a UPF0324 inner-membrane protein **conserved across *E. coli*, *Salmonella*,
*Shigella* and *Klebsiella*** (*Yersinia* differs at P1). The self hits are **GPER1**
(`GQMWLLAPR`), **PRPF3** (`TRLALIAPK`), **RNASEH2B**, JAK3, HELQ, IPO9, GLRB, MPP4
and SEC14L2 — GPER1 being expressed in cartilage, bone and the iris/ciliary body.
Recognition is β-centric: Tyr98/Phe98 of CDR3β contacts the peptide, CDR1α (Ala30,
Ile31, Tyr32) contacts P4, and CDR3α is minimally involved — which is why the J
region and the V scaffold matter as much as the CDR3 loop.

The named TCRs and their CDR3β cores: **AS2 / AS3.1 / AS4.1 = `VGLYSTDTQ`**,
**AS4.2 = `VGLFSTDTQ`**, **AS4.3 = `VATYSTDTQ`**. The mechanistic keystone, which the
[2026 review](https://pmc.ncbi.nlm.nih.gov/articles/PMC12963204/) draws out: the
**self**-peptide GPER1 activates these T cells when presented by disease-associated
HLA-B\*27:05 but **not** by protective B\*27:09, while **YeiH is presented by both**
— so the genetics of which B27 subtype causes disease tracks the autoantigen, not the
microbe. The reading is that a gut microbe primes the response and a self-peptide
sustains it. Caveats the authors state plainly: the peptides "are not known to be
linked to AS or AAU", several interactions are very weak, the AS4.1 tetramer failed
to stain, and some human peptides activated AS4.3 only at high concentration. **GEO
`GSE215018` from this paper is yeast-library sequencing (4 samples, *S. cerevisiae*,
PRJNA888081) — not patient repertoires.** Correspondence: `kcgarcia@stanford.edu`.

**Against this repo.** Three of their TCRs are sequences in our tables: AS4.1 is our
validated clone, AS4.2 is `CASSVGLFSTDTQYF` (Hanson's headline clone), and AS4.3 is
`CASSVATYSTDTQYF`, which our own screen surfaced at 4/31 AS, 0/10 healthy, 0/16 PsA —
so that candidate is a crystallographically characterised AS receptor, not a screening
artefact. Two structural facts cut against representations we have tried: recognition
runs through CDR3β position 98 and the TRAV21 CDR1α, so (a) the Y→F substitution at
the heart of our motif family is the contact residue itself, and (b) the α chain is
half the functional unit and is absent from every bulk TCRβ dataset we hold. This is
the predecessor of the [PRP paper](#deep-peptide-recognition-profiling-decodes-tcr-specificity-and-enables-disease-associated-antigen-discovery)
below — same lab, same yeast-display principle, scaled from one family to 21 TCRs.

---

## Mucosal signatures of pathogenic T cells in HLA-B*27+ anterior uveitis and axial spondyloarthritis

*Paley, Yang et al., JCI Insight 2024; 9(16):e174776.*
[10.1172/jci.insight.174776](https://doi.org/10.1172/jci.insight.174776) ·
[free at PMC11343591](https://pmc.ncbi.nlm.nih.gov/articles/PMC11343591/)

Located the motif-bearing cells' tissue of origin and built the first
antigen-specific blood assay for them.
Single-cell 5′ GEX+V(D)J on paired eye and blood from B27-associated acute anterior
uveitis (B27AAU) plus synovial fluid from 2 radiographic axSpA patients, 41 scRNA-seq
samples in total, against 1,535 expanded ocular CD8 clonotypes from non-B27 uveitis as
the comparator. Motif-bearing ocular CD8 cells share a 10-gene signature whose
standouts are **KLRB1 (CD161)**, SLC4A10, TMIGD2, HPGD and CTSH, and they match a
CD161-intermediate transcriptional state; expanded clonotypes in *other* uveitis do
not express KLRB1. In blood, motif-bearing clonotypes are confined to the CD161^int
fraction and undetectable in CD161^neg cells (*p* < 0.001, 4/4 samples). Using
HLA-B\*27 tetramers in 10 B27+ healthy controls vs 31 axSpA and/or B27AAU patients,
**YeiH-specific CD8 cells are expanded in patients** (and separately in axSpA alone
and B27AAU alone, *p* < 0.05) while cells specific for an influenza control peptide
(NP.383-391) are not. YeiH-specific cells carry **integrin α4β7, CCR6, CD103 and
CD161** — a gut-differentiated phenotype — supporting priming in the intestine
followed by trafficking to eye and joint. TRBV usage: **15 of 20 motif-bearing TCRs
(75%) use TRBV9**; the rest use TRBV5-4, TRBV5-5 or TRBV7-3, and `UV019.1`
(TRBV5-5) recognises YeiH without TRBV9 at all.

On the biomarker question the paper is explicit and its own numbers are modest:
tetramer-positive frequencies sit **below 0.1% of CD8 cells**, vary widely between
individuals, and were numerically *higher* during quiescence than after a recent
flare (not significant, and possibly confounded by NSAID/TNF-inhibitor use). Healthy
B27+ donors carry YeiH-specific cells too, and why they do not get disease is
unknown. No ROC or discrimination statistic is reported, only group comparisons. The
authors state the approach "will need refinement and adjustment for confounders",
and list as open: which intestinal compartment the cells occupy, whether they appear
before disease onset, and why TRBV9 is preferentially used. Data are on **figshare**,
open, no request needed:
[10.6084/m9.figshare.22078493.v1](https://doi.org/10.6084/m9.figshare.22078493.v1)
and [10.6084/m9.figshare.25796431.v1](https://doi.org/10.6084/m9.figshare.25796431.v1),
plus a "Supporting Data Values" file holding the per-point values behind every figure.
Correspondence: `paleym@wustl.edu`, `yokoyama@wustl.edu`.

**Against this repo.** Two things cut against our current representations. **TRBV9 is
characteristic but not required** — 5 of 20 motif-bearing receptors used TRBV5-4,
TRBV5-5 or TRBV7-3, and one recognised YeiH with no TRBV9 at all — so the planned
V-gene-conditioned similarity search must not be restricted to TRBV9. And the
invariant unit is TRAV21 + TRBJ2-3 + the `(Y/F)ST` core, two thirds of which we can
see and one third of which (the α chain) is absent from every bulk TCRβ dataset we
hold. Their blood result is also the closest existing answer to this project's own
question, measured against the right control group (B27+ donors) with an
antigen-specific reagent rather than sequence matching — but they publish no
discrimination statistic, and their figure-level source data is open, so computing
one is available to us.

---

## Supervised machine learning by combining immune receptor sequencing features yields accurate disease prediction (Mal-ID)

*Zaslavsky et al., Science 2025.*
[10.1126/science.adp2407](https://www.science.org/doi/10.1126/science.adp2407)

Reaches **0.986 multi-class AUROC over 593 people** across lupus, type 1 diabetes,
COVID-19, HIV and influenza, using a three-model ensemble combined by logistic
regression: (1) repertoire composition, (2) CDR3 edit-distance clustering, and
(3) off-the-shelf **frozen ESM-2 150M** embeddings. Models 1 and 3 performed better
on average than Model 2, and pairing them often matched the full ensemble. TCR data
contributed more for lupus and type 1 diabetes; BCR signals were clearer for
COVID-19, HIV and influenza. This is the proof that patient-level repertoire
classification works, and the template for the architecture.

**The detail that matters most is how Model 3 uses ESM-2, because it is easy to
copy wrongly.** They mean-pool ESM-2 hidden states to 640 dimensions, then train a
**separate supervised classifier per TRBV gene**, labelling each sequence with its
patient's disease status — noisy labels, as they note, since they do not know which
sequences are disease-related. Those per-sequence probabilities are aggregated
within each gene category (mean, median, trimmed mean, entropy-thresholded) into
patient-level features. They never use the embeddings as a **similarity metric**,
which is exactly what this project tried and why cosine distance returned anagrams:
a supervised head learns which of the 640 dimensions matter, while cosine weights
all of them equally and so measures amino-acid composition. Two further notes: the
authors attribute Model 2's weakness to **not modelling HLA genotype** (we have
HLA-A/B for 626 Emerson subjects, which they did not), and their known-binder
enrichment for TCRs reached only **0.56 AUROC** — a reason to expect little from
curated TCR–antigen databases.

---

## Deep peptide recognition profiling decodes TCR specificity and enables disease-associated antigen discovery

*Wang, Yeh & Garcia, Nature Biotechnology 2026.*
[10.1038/s41587-026-03128-x](https://doi.org/10.1038/s41587-026-03128-x) ·
[PubMed 42129507](https://pubmed.ncbi.nlm.nih.gov/42129507/) ·
code and data: [akds/PRP-TCR-Specificity](https://github.com/akds/PRP-TCR-Specificity)

Yeast-display selection against a **~10⁹-variant 9-mer peptide library** for **21
HLA-B\*27:05-restricted TCRs taken from AS and acute anterior uveitis patients**,
then a protein language model fine-tuned per TCR on the results. Their panel is the
**Faham `S··YSTDTQYF` family, measured** — it contains `CASSVGLYSTDTQYF` twice
(`4.1`, `19.2C5`) plus four of its one-edit neighbours, so it is the only
wet-lab-characterised sample of our exact motif family. Unique peptide ligands per
TCR ranged from hundreds to over 6,000. They report that these TCRs recognise
peptides **almost exclusively through CDR3β**, with β-centric docking geometry;
that models trained on these profiles beat AlphaFold3 and tFold-TCR at predicting
T-cell activation; and that a **PSG5**-derived peptide is a validated candidate
autoantigen, structurally confirmed against TCR 19.2.

**The concept to borrow is the metric.** A *peptide recognition profile* (PRP) is a
TCR's binding profile across a peptide panel, and PRP divergence — pairwise
Jensen–Shannon distance between position-wise binding motifs, excluding P2 (the B27
Arg anchor) and P8/P9 (fixed during library construction) — defines functional
similarity: two receptors are alike if they bind the same things, whatever they
spell. Their central generalisation result is that **predictive accuracy for a novel
TCR is governed by functional distance, not sequence similarity** (leave-one-out
AUPRC anticorrelates with PRP divergence at ρ ≈ −0.78, while CDR3 edit distance and
TCRdist both fail to predict it), and they offer **Mahalanobis distance** in the
learned embedding as a model-intrinsic confidence score for unseen TCRs (ρ ≈ 0.59
against measured PRP divergence). Two practical limits we hit: the **released
checkpoints cannot score novel CDR3βs** — the single-TCR models were fine-tuned
with β held constant, and the `19.2_allmutants` joint model learned β as a
near-pure intercept — and functional distance exists for only those 21 TCRs, so it
cannot be extended to our own candidates without their training data and a real
training run.
