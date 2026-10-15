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
