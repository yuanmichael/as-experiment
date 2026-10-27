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
probably the highest-value outreach available.

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
