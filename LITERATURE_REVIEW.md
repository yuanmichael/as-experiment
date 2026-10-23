# Which omics layer actually carries diagnostic signal in inflammatory disease?

A literature review, written to answer a practical question: after the AS project reached an AUC ceiling of ~0.78 in blood, is there a modality or compartment that would do better — and does the answer transfer to AS?

**Short version: compartment matters more than modality.** The same measurement technology, applied to the same patients, performs far better in a compartment that drains the diseased organ than in blood. AS has no such compartment, and that — not the choice of omics — is the binding constraint.

---

## 1. The headline result

The Integrative Human Microbiome Project (HMP2/IBDMDB) followed 132 subjects for a year, collecting 2,965 stool, biopsy and blood specimens with metagenomics, metatranscriptomics, metabolomics, host transcriptomics and serology ([Lloyd-Price et al., *Nature* 2019](https://www.nature.com/articles/s41586-019-1237-9)).

A later analysis of that resource compared classifiers built from the same modality in two compartments, in the same people:

| Compartment | Modality | AUC for IBD |
|---|---|---|
| **Stool** | metabolomics | **0.93** (95% CI 0.87–0.99) |
| **Plasma** | metabolomics | **0.74** (95% CI 0.69–0.79) |

Same disease, same subjects, same assay chemistry. Moving from blood to stool is worth ~0.19 AUC — far more than any modelling choice we tested in the AS project, where switching from a 10-gene hand-picked score to a trained elastic net on 4,703 genes bought 0.012.

Lloyd-Price et al. themselves emphasised *complementarity* rather than a winning modality: each measurement captured different components of dysbiosis. The practical reading is that modality choice is a second-order decision; compartment choice is first-order.

---

## 2. Diagnostic performance in IBD, by compartment

| Compartment | Test | Performance | Notes |
|---|---|---|---|
| **Colon tissue** | ileocolonoscopy + biopsy | CD vs UC correct in **89%**; revised later in 4%; **7–10% remain IBD-unclassified** | the gold standard, and still imperfect |
| **Stool** | faecal calprotectin | **sens 86%, spec 92%** vs Rome-criteria IBS (17 studies, 1,956 adults); AUC 0.97 at first presentation vs functional disease | [meta-analysis](https://pubmed.ncbi.nlm.nih.gov/37823411/) |
| **Stool** | metabolomics | AUC 0.93 | HMP2 |
| **Plasma** | metabolomics | AUC 0.74 | HMP2, same subjects |
| **Blood** | CRP, ESR, serology (ASCA/pANCA) | not used as a primary diagnostic; serology mainly subtypes established disease | no reliable pooled figures found |
| **Saliva** | none clinical | GSE232891 (552 samples, CD/UC/controls) is research only | unvalidated |

**Biopsy sampling matters more than people assume**: accuracy for colitis rises from 66% to 92% when segmental biopsies are taken rather than two sites.

### The prevalence lesson, again

Faecal calprotectin at sens 86 / spec 92 — a genuinely good, genuinely used test — gives at **1% IBD prevalence** in primary care: **NPV 99.8%, PPV 9%**.

This is the same arithmetic that gave our AS NK score a PPV of 17.7%. Calprotectin is not used as a rule-in diagnosis. It is a **rule-out gatekeeper** deciding who gets a colonoscopy. That is the realistic job description for essentially every non-invasive inflammatory biomarker, and it is a useful job.

---

## 3. Modality by modality in IBD

### Transcriptome (blood)
Subject to exactly the problem the AS project hit. Blood expression differences in inflammatory disease are substantially cell-composition effects plus graded programme shifts, with effect sizes around *d* ≈ 1 and a corresponding AUC ceiling near 0.76.

### Methylation (blood)
[Ventham et al., *Nat Commun* 2016](https://www.nature.com/articles/ncomms13507) — 240 newly diagnosed IBD cases vs 190 controls, 439 differentially methylated positions, 5 regions, with RPS6KA2, VMP1, ITGB2 and TXK replicating in an independent cohort.

Two findings matter for us:

- **RPS6KA2 turned out to be a CD14+ monocyte-specific signal**, and its strongest correlation was with **CRP (Spearman ρ = −0.53 in both cases and controls)**. It is an inflammation readout, not a disease-specific one — and "correlates with CRP in controls too" is close to a definition of non-specific.
- Follow-up analysis adjusting for cell proportions left 107 DMPs in CD. So a large share of the blood methylation signal is **cell composition**, which is the same composition-versus-regulation problem we worked through in the AS deconvolution.

Blood methylation in IBD is therefore *not* a clean counterexample to the blood transcriptome. It inherits the same confounders.

### Methylation (tissue)
Intestinal mucosal methylation reflects both intrinsic and extrinsic factors, and tissue methylation does separate IBD. But it requires the biopsy, which defeats the point of a minimally invasive test.

### Metabolome
Stool AUC 0.93 vs plasma 0.74 (above). Disrupted pools include acylcarnitines, bile acids and short-chain fatty acids. The stool advantage is unsurprising: faecal metabolites are the direct output of the gut lumen and its microbiota.

### Microbiome
Captured well by stool metagenomics/metatranscriptomics; inseparable from the stool-compartment advantage. Not applicable to AS except via the gut–joint axis hypothesis.

### Proteome (plasma, AS)
[Deep plasma proteomics, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12957791/): 104 samples (48 active AS, 21 stable AS, 35 HC), 3,891 proteins, 376 differentially abundant. Four candidates validated by ELISA: SAA1, FERMT3, ILK, TLN1 — integrin pathway plus an acute-phase protein. **No AUC, sensitivity or specificity reported anywhere in the paper.** SAA1 is serum amyloid A, which rises in any inflammation and carries the same specificity problem as CRP. Data deposited at iProX IPX0011977000, currently returning an empty record (embargoed).

---

## 4. cfDNA tissue-of-origin: the idea, and the experiment that already tested it

The most promising route on paper. DNA methylation is cell-type-specific and largely binary at informative loci, so a plasma cfDNA fragment can be assigned to the tissue it came from. This underpins multi-cancer early detection, and it detects a *misplaced tissue signature* rather than a graded shift — the right shape for a diagnostic.

**Atlases.** [Moss et al. 2018](https://www.nature.com/articles/s41467-018-07466-6) covered 25 tissues; [Loyfer et al., *Nature* 2023](https://www.nature.com/articles/s41586-022-05580-6) covered **39 cell types from 205 sorted samples** — blood, epithelial, vascular endothelial, fibroblast, muscle. **Neither includes chondrocytes, osteoblasts or entheseal tissue.** cfSort extends to 29 tissues by deep learning; MetDecode reports 84.2% correct tissue calls in cancer plasma.

**It works when cell death is substantial.** In acute graft-versus-host disease, colon-specific (SESN3) and liver-specific (PTK2B) methylation markers in cfDNA gave **AUC 0.952 and 0.971**.

**It largely failed in IBD.** A [2018 study from the same methodological lineage](https://www.biorxiv.org/content/10.1101/409219v1.full) built loci unmethylated specifically in intestinal epithelial cells and looked for them in plasma. The result: *"Patients with inflammatory bowel disease have minimal amounts of intestinal cfDNA in the plasma"* — signal in **2/16 IBD samples versus 2/45 healthy**, against a strong signal in advanced colorectal cancer. The authors suggested a more sensitive assay on a larger cohort might find something.

### What that contrast means

cfDNA tissue-of-origin detects **cell death at scale**, not inflammation. GvHD and cancer kill a lot of cells quickly. Chronic mucosal inflammation apparently does not produce enough net epithelial death to clear the detection threshold — in a tissue that is large, highly vascularised and directly in contact with the circulation.

An enthesis is small, poorly vascularised, low-cellularity fibrocartilage. If the gut mucosa fails this test, the enthesis is very unlikely to pass it.

---

## 5. The principle this all points to

Minimally invasive diagnostics that work measure a **discrete, foreign or misplaced analyte released into an accessible compartment**:

| Test | Analyte | Why it works |
|---|---|---|
| NIPT (e.g. BillionToOne) | fetal cfDNA | different genome; background structurally zero |
| dd-cfDNA (e.g. CareDx) | donor cfDNA | different genome; ~0.2% stable vs >1% rejection |
| ctDNA / MCED | somatic mutations, tissue-of-origin methylation | foreign sequence or misplaced tissue signature |
| Karius | microbial cfDNA | foreign organism |
| Faecal calprotectin | neutrophil protein in lumen | misplaced cell type, undiluted compartment |
| Anti-CCP (RA) | autoantibody | discrete, ~95% specific |

Tests built on **graded shifts in host cell state** — cytokine panels, CRP, immune-cell expression signatures, including ours — sit at *d* ≈ 0.5–1.0 and AUC 0.65–0.78, because the distributions genuinely overlap. For a single axis, maximum AUC is Φ(*d*/√2); reaching AUC 0.95 requires *d* ≈ 2.33.

Two things determine whether a discrete analyte is detectable:

1. **Is it released at all?** Requires cell death or active secretion at scale.
2. **How much is it diluted?** An inflamed gut mucosa (~200 g) shedding into 5 L of blood that drains every organ is a very different proposition from the same tissue shedding into stool, a compartment it essentially owns.

---

## 6. Implications for AS — a revision

Before this review I thought cfDNA tissue-of-origin was the most promising remaining route for AS. **I now think it is probably a dead end**, for a specific reason: the closest analogous experiment has already been run in a more favourable disease and largely failed. IBD, with a far larger and better-perfused inflamed tissue, yielded intestinal cfDNA in 2 of 16 patients.

What would have to be true for AS cfDNA to work, and what is known:

| Requirement | Status |
|---|---|
| Entheseal/cartilage/bone cell death at scale | unknown, but AS is characterised by new bone *formation* more than destruction |
| Chondrocyte/osteoblast methylation markers exist | **not in any published atlas** — would have to be built |
| Analyte reaches blood | entheses are poorly vascularised; unfavourable |
| Precedent in chronic inflammation | **negative** (IBD) |
| Precedent in arthritis | one 27-vs-18 pilot on DUSP22 promoter methylation in cfDNA (36.7% vs 46.9%, age-adjusted *p* = 0.049) — weak, and not tissue-of-origin |

Three of five are unfavourable and one is unknown. The honest conclusion is that AS lacks a compartment that drains the diseased tissue, and no choice of omics fixes that.

---

## 7. What I could not verify

- Per-modality AUCs from Lloyd-Price 2019 itself; the 0.93/0.74 comparison comes from a later analysis of the same data.
- Pooled sensitivity/specificity for IBD serology (ASCA/pANCA).
- Whether the 2018 intestinal cfDNA preprint was later published, and whether a larger follow-up changed the IBD conclusion. **This is the single most load-bearing citation in section 6 and should be chased before acting on it.**
- Whether chondrocyte or osteoblast methylation reference data exists anywhere in GEO/ENCODE/Roadmap. Not searched.

---

## 8. Where this leaves the project

The AS result stands as a mechanistic finding — NK cells make less cytotoxic cargo, and it is regulation rather than abundance. It is not, and will not become, a diagnostic.

If the goal is a diagnostic, the compartment argument says to pick a disease that **has** a draining compartment: gut disease (stool), airway disease (sputum, exhaled breath), urinary tract (urine), skin (tape strip). If the goal is understanding mechanism, AS remains interesting and blood remains usable.

Those are different projects and should not be confused with each other.

---
---

# Part II — the TCR repertoire thread

A separate question from Part I. Part I asked *which compartment*; this asks *which kind of feature*. The answer that emerged: discrete sequence features behave completely differently from graded expression measurements, and T cells escape the compartment problem entirely because they circulate rather than being shed.

## 9. Why GSE277117 could not have found TRBV9

10x 3′ chemistry. Reads land at the transcript's 3′ end — the TCR **constant** region — while the V segment sits ~1,000 bp upstream and CDR3 ~600 bp upstream. Measured directly:

| | total counts | % of cells |
|---|---|---|
| TRBC2 | 225,607 | 33.4% |
| TRBC1 | 168,422 | 24.1% |
| TRAC | 130,399 | 20.7% |
| TRBV9 | 388 | 0.11% |
| all 53 TRBV genes | 17,095 | — |

**V:C ratio 0.024.** No VDJ annotation in the object. Assay limit, not power limit — sub-clustering T cells would not have helped.

Secondary point: even perfect V-gene counts would not reproduce Faham. The finding is a CDR3 **sequence** carried on TRBV9, not TRBV9 abundance. A count matrix cannot hold a sequence.

## 10. KLRB1+/CCR6+ CD8 phenotype in GSE277117

Readable (ordinary gene expression), but thin.

- 1,731 KLRB1+CCR6+ cells in CD8 lineage (2.4%)
- **1,456 of them (84%) are MAIT or SLC4A10+** — a naive KLRB1+CCR6+ gate mostly re-finds MAIT cells
- 275 non-MAIT cells across 25 patients (~11/patient)
- exhaustion score −0.036 vs +0.006 in other CD8 (cell-level *p* = 0.0014, but 275 vs 71,668 cells — pseudo-replicated, 0.04 SD, not interpretable)
- patient level: cohort B *p* = 0.289; cohort A (12 v 4) *p* = 0.058, direction = fewer in non-responders

**Verdict:** direction consistent with the published phenotype, dataset cannot test it. Right question, wrong tissue — published work used synovial fluid. GSE288581 (AS synovial fluid + paired blood, 30 samples) is the dataset for this.

## 11. TRUST4 on bulk RNA-seq — recovery measured

GSE181364, 7 of 8 samples (orchestrator timed out; children completed). ~$2.

| sample | group | productive TRBβ clonotypes | TCRβ reads assembled | largest clone |
|---|---|---|---|---|
| AS_patient_1 | AS | 87 | 146 | 4.1% |
| AS_patient_2 | AS | 206 | 335 | 4.2% |
| AS_patient_3 | AS | 146 | 311 | 13.2% |
| AS_patient_4 | AS | 224 | 456 | 9.6% |
| AS_patient_5 | AS | 116 | 179 | 3.4% |
| Healthy_1 | HC | 292 | 493 | 5.3% |
| Healthy_2 | HC | 332 | 498 | 2.6% |

**Motif hits: 0** (exact `SVGLYSTDTQ` and partials `VGLYST`, `GLYSTDTQ`, `YSTDTQ`).

**Recovery: ~300 TCRβ reads from 34M sequenced = 0.0009%.** Prior estimate was 0.02% — off by 20×. Causes: whole blood is globin/neutrophil dominated, and TRUST4 needs overlapping coverage to assemble contigs, not just junction-spanning reads.

Detection floor ≈ **1% of repertoire**. We sampled ~150 of ~10⁶ clonotypes (0.015%).

**The null is uninformative.** It rules out only that the AS clone is among the largest clones in blood. It cannot distinguish "absent from blood" from "present at 0.1%".

Do not report TRBV9 *usage* (AS 2.79% vs HC 1.41%) as a finding — n=5 v 2 over ~150 clonotypes, and it is a graded gene-level measure of exactly the kind that does not discriminate.

## 12. Faham's design — and a correction

461 people (191 B27+ AS, 43 B27− AS, 227 controls), >77M clonotypes, Adaptive immunoSEQ (multiplex PCR on gDNA — inferred from affiliation, not stated in the abstract). Method: test 50 published sequences, then train/test split to discover motifs. 6 published + 11 discovered → 15 independent motifs.

**Correction to the earlier analysis in this project.** The headline "30–39% of cases, 0% of controls" compares mostly-B27+ cases against mostly-B27− controls. The motif-vs-disease comparison (B27+ AS vs B27+ healthy) is **p = 0.049**, and p = 0.001 for the published subset. So a substantial part of the headline specificity belongs to HLA-B27, which is already routinely tested. The earlier PPV estimate of 65–80% overstates incremental value over a B27 test.

## 13. Replicating the assay

Not a barrier. Off-the-shelf TCR kits exist (Takara SMARTer Human TCR a/b, QIAseq Immune Repertoire, iRepertoire, Invivoscribe LymphoTrack), or Adaptive immunoSEQ as a service. Library prep is community-lab feasible; sequencing is outsourced. Adaptive's advantage is targeted multiplex PCR in sample prep, not custom hardware. Note Adaptive holds IP on multiplex amplification of rearranged receptor loci — fine for research, a problem for a competing commercial assay.

**The bottleneck is sample collection, IRB and controls — never the assay.**

## 14. Public repertoire data inventory (AIRR Data Commons)

`scripts/airr_inventory.py`. 10 repositories queried (iReceptor ipa1–5, VDJServer, covid19-1–4). **5,356 repertoires, 56 diagnoses.**

COVID-19 accounts for 3,356 (63%); 968 are blank; cancers make up ~600.

**Inflammatory / autoimmune: 119 repertoires total.**

| | n |
|---|---|
| multiple sclerosis (all subtypes) | 93 |
| celiac disease | 10 |
| systemic lupus erythematosus | 10 |
| uveitis | 4 |
| neurosarcoidosis | 2 |
| **ankylosing spondylitis** | **0** |
| **IBD / Crohn's / UC** | **0** |
| **rheumatoid arthritis** | **0** |

Caveats, all of which understate the true total: **immuneACCESS (Adaptive) is not in the ADC** and holds substantially more — Emerson's 666 CMV repertoires, probably the Faham AS data; needs an account. E-MTAB-11498 (Komech, 115 samples, 29 AS + 86 PsA) is in ArrayExpress/ENA, outside the ADC. The "0 healthy controls" result is a labelling artefact; controls are likely among the 968 blanks.

**Conclusion drawn at the time — since retracted, see §15:** *a multi-disease inflammatory repertoire model cannot be built from public data today.*

## 15. Raw-archive sweep — the §14 conclusion was wrong

`scripts/repertoire_inventory.py`. GEO + ENA + BioStudies/ArrayExpress instead of the ADC. **2,715 candidate studies, 1,002 mentioning an inflammatory or autoimmune disease.**

The §14 conclusion was an artefact of querying one curated index. The ADC indexes *curated* repertoire studies; the archives hold the rest, filed as ordinary sequencing. Two systematic blind spots:

1. **ArrayExpress/ENA deposits.** Komech's own data sits here, invisible to the ADC.
2. **10x 5′ libraries with V(D)J enrichment.** Indexed as "single-cell RNA-seq." These are dedicated, deep TCR sequencing *with a paired transcriptome*, and no repertoire inventory will ever surface them.

### AS / SpA — ADC said 0; the real count is 6 studies

| accession | design | size |
|---|---|---|
| **PRJNA393498** | Komech 2018 *Rheumatology* (`kex517`). 5′RACE TCRβ, PB + SF, TRBV9-sorted and CD8-depleted subsets. **Paper reports 25 AS + 108 healthy donors.** | 78 runs, 44.5 GB, 378M reads |
| PRJEB51315 = E-MTAB-11498 | Komech 2022 *Front Immunol*. 29 AS + 86 PsA, 64 SF + 51 PB. No healthy controls. | 115 runs, 47.1 GB, 631M reads |
| GSE288581 | AS, mass cytometry + scRNA-seq, clonally expanded CTL subset in SF | 30 samples |
| E-MTAB-10948 | AS scRNA-seq, PB + SF, 2 HLA-B27+ patients, >55,000 cells. No VDJ stated. | — |
| PRJNA853355 = GSE207037 | single-cell TCR, DC-stimulated CD4 T cells, AS | 1 run, 3.9 GB |
| GSE290921 | SpA, CD4 TRM Th17 driving IL-17A joint pathology | 6 samples |

**PRJNA393498 is the find.** Same lab and naming convention as E-MTAB-11498 (`Dv_PB_4`, `Kal_SF_8`, `Kos_SF_TRBV9`, `LK-BV9`), one study earlier, and its paper has the healthy-donor arm that E-MTAB-11498 lacks. Caveat: 78 runs is far short of 25 + 108, so the control repertoires are probably a previously published Chudakov-lab healthy cohort deposited separately — needs checking before it can be called an AS-vs-healthy dataset.

Adjacent and directly on-target: **E-MTAB-9492** — 10x gene expression *plus* V(D)J from blood, synovial fluid and synovial tissue in psoriatic arthritis. Only 5 patients (3 paired PB/SF, 2 synovial tissue), but it is paired sequence-and-transcriptome in a spondyloarthritis across three compartments.

### Cross-disease

| disease | dedicated repertoire | paired single-cell | standouts |
|---|---|---|---|
| IBD | 33 | 8 | PRJEB27352 — monozygotic twins, genetics-controlled |
| celiac | 11 | 0 | PRJNA650391, PRJNA678347, PRJNA1330789 — all three already ran Faham-style motif/classifier analyses |
| type 1 diabetes | 8 | 5 | GSE272431 — >2×10⁸ TCRβ from sorted CD4 subsets, patients + healthy donors; GSE178991 — 10x GEX+VDJ, n=84 |
| rheumatoid arthritis | 6 | 3 | GSE256256 — seven-chain AIRR, n=112, controls, TNFi phenotypes |
| AS / SpA | 5 | 1 | above |
| multiple sclerosis | 5 | 0 | GSE121082 is actually HTLV-I/HAM-TSP — matcher false positive |

**Revised conclusion: a multi-disease inflammatory repertoire model is buildable from public data.** The constraint is curation, not existence.

Known weakness of the sweep: `dedicated` is inferred from title text, and "AMPLICON" matches 16S microbiome studies — a substantial share of the IBD and psoriasis hits are gut-microbiome papers. The flag generates candidates; titles still need reading, and clonotype recovery per sample (≥10⁴ is the real bar) can only be measured from the data.

### Paired GEX + V(D)J — training data for a cell-level model

`scripts/vdj_paired_inventory.py`. Text queries have poor recall (GSE288581 deposits `filtered_contig_annotations.csv.gz` *and* `matrix.mtx.gz` but never says "VDJ" in its summary), so the sweep casts wide and lets GEO `filelist.txt` decide: a VDJ-shaped file and a GEX-shaped file in one series means the 10x cell barcode can join sequence to transcriptome.

**Of 450 series checked out of 1,210 candidates: 122 have VDJ files, 197 have GEX, 86 have both — 6,008 samples, 500 GB of *processed* data.** Median series is 1–5 GB, so a usable training set is ~20 series and tens of GB, not the hundreds of GB the raw FASTQs would cost. Largest: GSE272993 (360 samples, checkpoint-blockade clonal responses), GSE161918 (241, COVID), GSE270518 (201, airway memory). Mostly cancer, COVID, transplant — which is fine, because the cell-level label is the cell's own expression state, not a diagnosis.

Caveat: "both file types present" is not "barcodes overlap." That has to be checked by opening one of each, and is the first thing to verify.

## 16. Open items

- Check immuneACCESS directly (account required) for AS/IBD repertoires and the Faham deposit.
- **PRJNA393498: locate the 108 healthy-donor repertoires.** If they are retrievable, this is the first public AS-vs-healthy dedicated repertoire comparison and the Faham replication becomes possible rather than hypothetical.
- **E-MTAB-9492: confirm the GEX↔VDJ barcode pairing is intact.** If so it is training data for a cell-level `TCR sequence → cell state` model in a spondyloarthritis. n=5 is too small to train on; it is a test set.
- E-MTAB-11498 would show what frequency AS clones sit at in blood versus synovial fluid — decides whether shallow sequencing has any chance. No healthy controls.
- GSE221786 (20 AS / 8 HC, IL-17-enriched sorted T cells, 222 GB, est. $8–15): T-cell enriched, so recovery should be far better than whole blood. A positive would be informative; a null still would not. Not run.
- GSE288581 — AS synovial fluid with paired blood, for the KLRB1/CCR6 phenotype.
