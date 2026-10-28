# Ankylosing spondylitis: what a T-cell repertoire can and cannot tell you

Independent reanalysis of public AS data, done from scratch. Written by a patient learning the methods, so the wrong turns are recorded next to the results — several of them were nearly published as findings.

The project asks one question: **can a blood draw predict ankylosing spondylitis?** It has gone through three arms. Bulk transcriptomics found a real but graded signal that makes a poor test. A cell-level receptor→expression model found that a receptor explains 1–2% of its own cell's state. The third arm, patient-level repertoire analysis, is where the useful result is, and it is the one this README leads with.

---

## The AS clonotype

`CASSVGLYSTDTQYF` is a TCRβ CDR3 reported in AS by Faham et al. (2017) and Komech et al. (2018). We tested it in every public repertoire we could reach: four cohorts, four assays, 836 blood repertoires.

| cohort | assay | group | n | present | ≥10⁻⁵ | max freq | median depth | floor |
|---|---|---|---|---|---|---|---|---|
| Komech 2018 | 5′RACE amplicon, PBMC | **AS** | 25 | **12 (48%)** | 6 | 9.3×10⁻⁵ | 5.29M junction reads | 1.9×10⁻⁷ |
| Komech 2018 | same | healthy | 10 | **0 (0%)** | 0 | 0 | 3.71M | 2.7×10⁻⁷ |
| PRJEB51315 | amplicon, blood | **AS** | 7 | **5 (71%)** | 4 | 2.2×10⁻³ | — | 6×10⁻⁸–2.4×10⁻⁵ |
| PRJEB51315 | amplicon, blood | PsA | 16 | **4 (25%)** | 1 | 5.7×10⁻⁵ | — | 1.7×10⁻⁷–9.2×10⁻⁶ |
| Emerson 2017 | gDNA immunoSEQ | healthy | 786 | **15 (1.9%)** | **0** | 9.7×10⁻⁶ | 377K templates | 2.7×10⁻⁶ |
| Emerson 2017 | gDNA, B27+ subset | healthy B27+ | 45 | **6 (13.3%)** | 0 | — | 377K | 2.7×10⁻⁶ |
| GSE288581 | 10x, sorted CD8 memory | **AS** | 4 | **3 (75%)** | — | 4.0×10⁻³ | 2,152 cells | 4.7×10⁻⁴ |
| GSE288581 | same | healthy | 4 | **0 (0%)** | — | 0 | 2,040 cells | 4.9×10⁻⁴ |

Pooled blood: **AS 20/36 (56%), PsA 4/16 (25%), healthy 15/800 (1.9%)**.

Within-study, where assay and depth are matched by construction:

- Komech AS 12/25 vs its own healthy 0/10 — **p = 0.0062**
- PRJEB51315 AS 6/7 vs PsA 6/17, synovial fluid — **p = 0.034**
- PRJEB51315 AS 5/7 vs PsA 4/16, blood — p = 0.052
- GSE288581 AS 3/4 vs healthy 0/4 — p = 0.071 (n = 4; descriptive only)

Pooled, crossing platforms: AS vs all healthy p = 2×10⁻²¹; **AS vs healthy B27+ only p = 5.9×10⁻⁵**; AS vs PsA p = 0.040. The B27+ comparison is the one that matters, because B27 is the disease-associated allele and a clonotype that merely tracks B27 adds nothing to a routine B27 test.

### Depth decides whether absence means anything

Komech's libraries are **14× deeper** than Emerson's (5.29M junction reads vs 377K templates). Presence is depth-sensitive — a clone at 10⁻⁶ gives ~5 reads at one depth and ~0.4 templates at the other — so a 48%-vs-1.9% gap between unequally deep cohorts is not trustworthy on its face. Three checks, from both directions:

1. **The within-study control is already depth-matched.** Komech's 10 healthy controls ran at the same depth as its AS patients and show 0/10, with floors of 1.4×10⁻⁷ to 8.6×10⁻⁷. Every one could have detected a single read. Those zeros are real absences.
2. **Thin the AS patients down to Emerson's depth and the signal survives.** Binomial downsampling to 377K templates gives an expected 9.0/25 present (36%) against Emerson's 1.9%; p = 5×10⁻⁹.
3. **Restrict Emerson up to its deepest repertoires.** The 32 subjects with ≥2M templates — genuinely depth-matched to Komech — carry it in **1 of 32 (3.1%)**. Against Komech AS 12/25, p = 7×10⁻⁵.

The gap is not a depth artefact in either direction.

### As a diagnostic

Presence as a binary test for AS, with PPV and NPV at four prevalences:

| control group | sensitivity | specificity | PPV @5% | PPV @20% | PPV @50% | NPV @20% |
|---|---|---|---|---|---|---|
| all healthy (n=800) | 55.6% (38–72%) | 98.1% (96.9–98.9%) | 60.9% | 88.1% | 96.7% | 89.8% |
| **healthy B27+ (n=45)** | 55.6% | **86.7% (73–95%)** | 18.0% | 51.0% | 80.6% | 88.6% |
| PsA (n=16) | 55.6% | 75.0% (48–93%) | 10.5% | 35.7% | 69.0% | 87.1% |

Read it this way. **Sensitivity ~56% means the test misses nearly half of AS patients**, so it cannot rule out disease — NPV only looks good at low prevalence, where almost anything does. Specificity against unselected healthy people is excellent (98.1%), which is why PPV reaches 88% in a population where one in five has AS.

But the clinically honest control is not an unselected healthy person. A patient being worked up for AS is typically young, with inflammatory back pain, and often already known to be B27+. Against B27+ healthy controls specificity falls to **86.7%** and PPV at 20% prevalence falls to **51%** — a coin flip. Against PsA, the differential that actually presents in a rheumatology clinic, specificity is 75%.

So: **a positive result is informative, a negative result is nearly worthless, and most of the apparent power evaporates once you compare against the patients who would actually be tested.** It is a rule-in test at best, and the B27+ and PsA control arms are built on 45 and 16 subjects respectively, so even those numbers are provisional.

Two routes could improve it. Sensitivity is depth-limited — half the AS carriers sit between 10⁻⁷ and 10⁻⁵, so deeper sequencing should recover more of them. And the clone is CD8-restricted (below), so sorting CD8 cells before sequencing concentrates it. Neither fixes specificity against B27+ controls, which is the real ceiling.

A note on the within-study numbers: Komech alone gives 48% sensitivity at 100% specificity, which would imply PPV = 100%. That is an artefact of 0 positives in 10 controls — the 95% CI on that specificity runs from 69% to 100%. Ten controls cannot establish a specificity that matters.

### The clone is CD8, which is what B27 predicts

Re-derived on the honest denominator across all 64 synovial libraries:

| group | subset | target freq | junction reads |
|---|---|---|---|
| **AS** | **CD8** | **1.00×10⁻³** | 4,378,333 |
| AS | unsorted | 2.69×10⁻⁴ | 3,569,466 |
| AS | CD4 | 2.88×10⁻⁵ | 2,533,848 |
| PsA | CD8 | 8.78×10⁻⁵ | 10,914,688 |
| PsA | CD4 | 1.56×10⁻⁶ | 6,397,279 |

**AS joint CD8 sits 35× above AS joint CD4**, and 11× above PsA joint CD8. HLA-B27 is an MHC class I molecule; class I is read by CD8 T cells and class II by CD4. A genuinely B27-driven clone therefore *must* be CD8, and it is. No technical story predicts a compartment restriction that matches the immunology, which makes this the strongest evidence that the sequence is a real B27-restricted clone. (The 35× figure survived the denominator correction unchanged — it was computed within-subset both times.)

### In the joint, the clone does reach 10⁻⁴

The withdrawn ≥10⁻⁴ claim was about blood. Synovial fluid is a different story, because frequencies there are 10–100× higher:

| group | n | present | ≥10⁻⁴ | max freq |
|---|---|---|---|---|
| **AS synovial** | 7 | 5/7 | **4/7** | 2.34×10⁻³ |
| PsA synovial | 17 | 6/17 | 2/17 | 3.08×10⁻⁴ |

Mann-Whitney on frequency, AS vs PsA joints: **p = 0.019** — threshold-free, and essentially unchanged from the pre-correction value, because both arms were inflated by the same mechanism. Every patient is evaluable (floors 1.3×10⁻⁷ to 2.2×10⁻⁵). AS joints reach 1.0×10⁻³–2.3×10⁻³ in three patients; no PsA joint exceeds 3.1×10⁻⁴.

### The clone is one member of a sequence family

[Hanson et al. 2020](https://eprints.qut.edu.au/210039/) define the AS-associated motif as **`CASSVG(V/I/L)(Y/F)STDTQYF`**, and their strongest single clone is the **F** variant, not ours: `CASSVGLFSTDTQYF` in 12/37 B27+ AS patients and 0 of all B27− patients and all controls (p = 1.6×10⁻⁵). The union of the F and Y variants reaches 18/37 (49%) B27+ AS versus 2/19 (11%) B27+ healthy (OR 8.05, p = 0.007). [Yang et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC10511244/) solved structures for three members, which they call AS4.1 (`VGLYSTDTQ`), AS4.2 (`VGLFSTDTQ`) and AS4.3 (`VATYSTDTQ`).

Counts below are presence at ≥2 reads after thinning every subject to 800,000 junction reads, as built by [`scripts/build_panel.py`](scripts/build_panel.py) — **not** the native-depth presence used in the table above, which is why the known clone reads 11 here and 12 there. Provenance is tracked separately because only literature-specified sequences can enter a leakage-free test.

**Literature-specified** — named by Hanson's motif or by Yang's solved structures, before we looked:

| cdr3 | source | Komech AS /25 | Komech healthy /10 | Emerson healthy /786 | previously reported |
|---|---|---|---|---|---|
| `CASSVGLYSTDTQYF` | Hanson motif; Yang AS4.1 | 11 | 0 | 15 (1.91%) | yes — the main result |
| `CASSVGLFSTDTQYF` | Hanson motif **and headline clone**; Yang AS4.2 | **6** | **0** | 20 (2.55%) | **no** |
| `CASSVGVYSTDTQYF` | Hanson motif | **4** | **0** | 11 (1.40%) | **no** |
| `CASSVGVFSTDTQYF` | Hanson motif | <4 | — | — | — |
| `CASSVGIYSTDTQYF` | Hanson motif | <4 | — | — | — |
| `CASSVGIFSTDTQYF` | Hanson motif | <4 | — | — | — |
| `CASSVATYSTDTQYF` | Yang AS4.3 (structure solved) | 4 | 0 | 7 (0.89%) | yes — as a candidate |

**Exploratory** — one edit from a specified member, surfaced by our own screen, so not admissible to a pre-specified test:

| cdr3 | Komech AS /25 | Komech healthy /10 | Emerson healthy /786 |
|---|---|---|---|
| `CASSLGLFSTDTQYF` | **7** | **0** | 31 (3.94%) |
| `CASSLGLYSTDTQYF` | **4** | **0** | 33 (4.20%) |

`<4` means below the panel build's four-carrier floor, i.e. not yet counted rather than absent. Three sequences (`VGLY`, `VGLF`, `LGLF`) were inside the 167-sequence burden panel, each weighted as 1/167 of a count; the rest were below its ≥5-carrier cut. The 3-mer composition model never represented any of them individually. **The family has never been scored as a unit.**

### The union of the family clears significance where one sequence does not

[`scripts/motif_family_test.py`](scripts/motif_family_test.py) scores the seven literature-named sequences as a unit, reusing `clono_sweep`'s loader, thinner and carrier matrix so presence means exactly what it means in `build_panel.py`. The panel was fixed by Hanson and Yang before we looked, so no membership decision saw our labels.

| test (thinned to 800K) | AS /25 | healthy /10 | sensitivity | one-sided Fisher p |
|---|---|---|---|---|
| `CASSVGLYSTDTQYF` alone — the previous headline | 11 | 0 | 44.0% | 0.0107 |
| **union of the 7 pre-specified sequences** | **14** | **0** | **56.0%** | **0.0019** |
| union, at native depth | 15 | 0 | 60.0% | 0.0010 |

The union adds **three AS patients the known clone misses**, and it is the first open-ended result in this project to clear α = 0.01 in the within-study comparison. That matters because of the power arithmetic above: with 25 cases and 10 controls a marker needs ≥12 AS carriers and 0 controls before significance is even reachable, and the single clone sat at 11 after thinning. Threshold-free, the carrier count (0–7) gives Mann-Whitney **p = 0.0020**, rank-biserial AUC **0.780** — AS median 1 (range 0–3), healthy median 0 (range 0–0, no control carries any member).

Per sequence, thinned:

| cdr3 | AS /25 | healthy /10 | source |
|---|---|---|---|
| `CASSVGLYSTDTQYF` | 11 | 0 | Hanson motif; Yang AS4.1 |
| `CASSVGLFSTDTQYF` | 6 | 0 | Hanson motif **and headline clone**; Yang AS4.2 |
| `CASSVGVYSTDTQYF` | 4 | 0 | Hanson motif |
| `CASSVATYSTDTQYF` | 4 | 0 | Yang AS4.3 (structure solved) |
| `CASSVGVFSTDTQYF` | 2 | 0 | Hanson motif |
| `CASSVGIYSTDTQYF` | 0 | 0 | Hanson motif |
| `CASSVGIFSTDTQYF` | 0 | 0 | Hanson motif |
| `CASSLGLFSTDTQYF` | 7 | 0 | *exploratory* — one edit, found by us |
| `CASSLGLYSTDTQYF` | 4 | 0 | *exploratory* — one edit, found by us |

Neither `VGIY` nor `VGIF` occurs in any of the 35 subjects, so Hanson's six-member grid is effectively a four-member family in this cohort. `CASSLGLFSTDTQYF` at 7 carriers has more support than any member except the known clone, but we found it, so it stays out of the primary test and is reported here only.

**What this does not yet establish.** The control arm is 10 Komech donors of unknown HLA-B27 status, which cannot estimate specificity: 0/10 is consistent with anything up to ~26% carriage. The Emerson arm is the specificity test and needs per-subject repertoires, which are not currently on disk — per-sequence carriage bounds the union at 2.5–6.7% of 786 healthy donors, against 1.9% for the known clone alone, so the union certainly costs specificity and the question is how much. Hanson's B27-matched measurement of the two commonest members, 18/37 versus 2/19, is the closest existing answer.

---

## The best blood discrimination anyone has measured is 0.84, and it is not a sequence test

[Paley et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11343591/) stained blood with an HLA-B27 tetramer loaded with **YeiH**, the bacterial peptide [Yang 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC10511244/) showed this receptor family recognises, and counted the T cells that bind it — no sequencing, no string matching. They report group comparisons but no discrimination statistic. Their figure-level source data is open, so [`scripts/paley_tetramer_roc.py`](scripts/paley_tetramer_roc.py) computes one from the per-subject values behind Figures 4B and 4C:

| measurement | n | AUC (95% CI) | Mann-Whitney p |
|---|---|---|---|
| **YeiH-specific, % of CD8** | 31 patients vs **10 B27+ healthy** | **0.837 (0.673–0.961)** | 0.0008 |
| YeiH-specific, per 10⁶ PBMC | same | 0.829 (0.645–0.974) | 0.0010 |
| influenza NP-specific, % of CD8 — *specificity control* | same | 0.642 (0.473–0.805) | 0.093 |
| influenza NP-specific, per 10⁶ PBMC — *same control* | same | 0.719 (0.545–0.871) | 0.020 |

It holds in each disease separately against the same controls: B27+ uveitis AUC 0.841 (p = 0.0046), B27+ axSpA AUC 0.850 (p = 0.0045).

**0.84 against B27-matched controls is the number to beat**, and it beats everything in this repo — the bulk NK score (0.77), the composition model (0.787 against unmatched controls), the motif-family carrier count (0.780 against 10 donors of unknown B27 status). It is also measuring the right thing: antigen specificity rather than sequence identity, so it is immune to depth, batch and the exact-match problem that has cost this project three representations.

**Two caveats that keep it provisional.** The specificity control is not clean — influenza NP tetramers reach AUC 0.642–0.719 and p = 0.020 on one of the two scales, so some of the YeiH separation may be a general difference in the patients' CD8 compartment rather than YeiH specificity itself. The per-person YeiH/NP ratio would settle it and the source data may support it if the two sheets are ordered by subject, which is unverified. And 10 controls cannot estimate specificity: the best Youden cut gives 97% sensitivity at 60% specificity, which is a wide-CI guess, not an operating point.

---

## Is anything else there? One candidate, and the control arm is everything

The clonotype above came from someone else's hypothesis. The harder question is whether an **open-ended** sweep of millions of clonotypes finds anything separating AS from healthy. Komech 2018 is the right cohort: 25 AS and 10 healthy, one lab, one assay, 6.6M (subject, clonotype) rows, 5.3M distinct clonotypes.

### The within-study sweep finds nothing, because it cannot

13 clonotypes reach p ≤ 0.01 against **26.9 expected under 200 label shuffles** — fewer than chance, with the null ranging 14 to 80. The expansion sweep gives 1 against 2.8 expected. No aggregate repertoire feature differs: clonality, Shannon entropy and distinct-clonotype count are flat (p ≥ 0.25) after depth-thinning.

That negative is almost entirely a power statement, and the arithmetic says so exactly. With 25 cases and **10 controls**, the smallest p-value attainable by a clonotype carried by *a* AS patients and zero controls is fixed:

| AS carriers | best possible p | clears α = 0.01? |
|---|---|---|
| 4 | 0.242 | no |
| 8 | 0.046 | no |
| 11 | 0.0107 | no |
| **12** | **0.0062** | **yes** |

So **a clonotype must sit in ≥12 of 25 AS patients (48%) and 0 of 10 controls before significance is even reachable.** Nothing I chose imposed that; ten controls do. The funnel:

| stage | clonotypes |
|---|---|
| distinct clonotypes observed | 5,329,665 |
| after thinning to 800K reads/subject | 2,678,691 |
| in ≥1 AS patient, 0 controls | 1,685,436 |
| in ≥4 AS patients, 0 controls | 7,019 |
| in ≥6 AS patients, 0 controls | 1,081 |
| **in ≥12 AS patients, 0 controls — the only testable zone** | **6** |

**1,075 clonotypes sit in 6–11 AS patients and no control, and the 10-control arm can never call any of them significant.** That is what the within-study null is actually reporting.

The sweep's positive control confirms it: our known clone is recovered at rank 160 of 33,846 — top 0.5% — at p = 0.046, which does not clear α. **This cohort is underpowered to have discovered its own clonotype de novo.**

### Every apparent hit is a public clonotype

Checking the 14 top within-study hits against the 786 Emerson healthy repertoires:

| cdr3 | AS /25 | Komech healthy /10 | Emerson healthy /786 |
|---|---|---|---|
| CASSLAGYGYTF | 12 | 0 | **664 (84.5%)** |
| CASSFGTGGYEQYF | 16 | 1 | 559 (71.1%) |
| CASSPRDSYEQYF | 12 | 0 | 551 (70.1%) |
| CASSYSSNQPQHF | 15 | 1 | 548 (69.7%) |
| CASSLGLEQYF | 18 | 1 | 466 (59.3%) |
| CASSPGGPYEQYF | 13 | 0 | 497 (63.2%) |
| …9 more | 12–15 | 0–1 | 278–552 (35–70%) |
| **CASSVGLYSTDTQYF** | 12 | 0 | **15 (1.9%)** |

Every one of the fourteen is present in 35–85% of healthy people. They look AS-specific only because ten healthy controls are too few to exclude a clonotype that two-thirds of humanity carries. Our target is the sole exception at 1.9%.

### Swapping in 786 controls: one clonotype survives out of 1,081

The fix for both problems is the same — discover candidates in Komech's 25 AS patients, then control against Emerson's 786 healthy repertoires instead of Komech's 10. Screening all **1,081** clonotypes carried by ≥6 AS patients and no Komech control:

| how common in 786 healthy | candidates | share |
|---|---|---|
| absent from all 786 | 0 | 0.0% |
| <1% | **1** | 0.1% |
| 1–5% | 20 | 1.9% |
| 5–10% | 81 | 7.5% |
| **>10%** | **979** | **90.6%** |

**90.6% of the candidates are carried by more than a tenth of healthy humanity.** Exactly one is rare, and it is new:

| | AS present | healthy present | max freq | Fisher vs healthy | ×1,081 |
|---|---|---|---|---|---|
| CASSVGLYSTDTQYF (known) | 12/25 (48%) | 15/796 (1.9%) | 9.3×10⁻⁵ | 4.6×10⁻¹³ | 5.0×10⁻¹⁰ |
| **CASSFGTGGINEQYF (new)** | **7/25 (28%)** | **5/796 (0.6%)** | **2.2×10⁻⁴** | 7.1×10⁻⁹ | 7.7×10⁻⁶ |

Both survive Bonferroni correction for all 1,081 tests. The new clonotype reaches a *higher* peak frequency than the known one (2.2×10⁻⁴ in patient P15).

**It is not a sequencing artefact.** It sits 2 edits from `CASSFGTGGYEQYF`, which 71% of healthy people carry, so the obvious worry is a PCR error variant of an abundant parent — which would favour deeper libraries, and AS libraries are deeper. But it reaches 1,068 reads in P15 and 535 in P17 against parent counts of 103 and 6, i.e. **10× and 89× *above* the supposed parent.** An error variant cannot outnumber its parent 89-fold.

Treat it as a candidate, not a finding: it is one cohort, 7 patients, discovered by a screen whose control arm is a different platform (Emerson is gDNA immunoSEQ, Komech is RNA amplicon), and its frequency distribution does not separate from Komech's own controls (Mann-Whitney p = 0.11). It needs replication in PRJEB51315 and the synovial arm before it means anything.

The methodological result is the durable one: **a ten-control arm cannot distinguish a disease-specific clonotype from a public one, and a 786-control arm does it in one pass.** Every future candidate is screened against Emerson before anything else is done with it.

It partially replicates in the joint: `CASSFGTGGINEQYF` is present in 2 of 7 AS synovial samples, both above 10⁻⁴ (1.1×10⁻⁴ and 1.0×10⁻⁴), against 1 of 17 PsA joints at 2.6×10⁻⁶. Mann-Whitney p = 0.061 — the right direction, not yet significant.

### Two control arms, not one

Screening wider surfaced a second artefact class that one control arm cannot catch. Of 21,206 clonotypes carried by ≥5 AS patients, only 9 are rare (≤1%) in Emerson — and five of those are a **batch family**:

| cdr3 | AS /25 | Komech healthy /10 | Emerson /786 |
|---|---|---|---|
| CRVSGL**W**AKNIQYF | 7 | 3 | 0 |
| CRVSGL**G**AKNIQYF | 6 | 3 | 0 |
| CRVSGL**C**AKNIQYF | 6 | 4 | 0 |
| CRVSGL**R**AKNIQYF | 5 | 1 | 0 |
| CSWGTHDALCSNQPQHF | 6 | 5 | 0 |

Absent from all 786 external repertoires, but carried by 30–50% of this deposit's *own healthy controls*. A one-residue ladder at a fixed position, deposit-specific and disease-indifferent. Tested against Emerson alone they reach p = 10⁻¹¹; tested against Komech's healthy controls they reach p = 0.71. **A candidate must clear both arms** — the external arm kills public clonotypes, the internal arm kills batch artefacts, and neither does the other's job.

### A multi-clonotype burden score: retracted

Thresholding one clonotype at a time is arbitrary — "≥6 carriers and 0 controls" is a cut with no principle behind it — and it ignores that weak signals might combine. The alternative tried here was a burden score: build a panel of candidate sequences, then count how many of them a subject carries.

The panel is [`data/panel/as_burden_panel_167.tsv`](data/panel/as_burden_panel_167.tsv), rebuilt by [`scripts/build_panel.py`](scripts/build_panel.py). A CDR3 is kept if it appears in **≥5 of the 25 Komech AS patients** (presence = ≥2 reads, after thinning every subject to 800,000 junction reads) **and** in **≤5% of Emerson 2017's 786 healthy repertoires**. That yields 167 sequences; the score is the fraction a subject carries. A wider panel at ≥4 carriers gives 1,014 sequences ([`as_burden_panel_1014.tsv`](data/panel/as_burden_panel_1014.tsv)).

**This reached AUC 0.962 within Komech, 0.991 against an independent healthy cohort, and 84% sensitivity at 100% specificity. All of those numbers were inflated by selection leakage and are withdrawn.**

The leak: the candidate pool was defined as clonotypes appearing in ≥5 of **all 25** AS patients, including the subject being held out. Only the second-stage carrier filter was fold-restricted. So for a held-out AS patient, the sequences they personally carry were preferentially in the pool, and their score was inflated by construction. The same flaw affects the external comparison, where in-sample AS scores were compared against out-of-sample healthy donors.

Rebuilt with the pool constructed from the 24 training patients only:

| panel rule | median size | leak-free LOO AUC | leak-free AS vs Britanova |
|---|---|---|---|
| ≥4 carriers | 907 | 0.528 | 0.680 |
| **≥5 carriers** | 150 | **0.672** | 0.803 |
| ≥6 carriers | 37 | 0.620 | 0.803 |

Permutation-tested with the carrier grid inside the null: **observed 0.672, null mean 0.561, 95th percentile 0.744, p = 0.17.** Not significant. The best operating point is 60% sensitivity at 80% specificity.

The cross-deposit 0.803 is not clean either. Komech's healthy controls score higher than Britanova's (0.054 vs 0.036) because they share deposit-specific sequences with Komech's AS patients — so part of that number is deposit discrimination, batch happening to work in the result's favour. The within-deposit figure is the one to trust, and it is null.

**There is no demonstrated multi-clonotype signal.** Two things are worth keeping from the attempt. The score *is* continuous rather than binary — AS patients carry a median of 38 of the 167 and healthy controls 6–9, with nobody above 56 — so "do you have the panel" was never the question; and the control arm dominates everything, since the same score built against Komech's 10 controls instead of Emerson's 786 lands at AUC 0.44.

What survives untouched is every result for `CASSVGLYSTDTQYF`, which was specified from the literature before any of this data was examined, so no selection leak is possible. `CASSFGTGGINEQYF` was discovered by screening and is Bonferroni-corrected, but its carrier count was used both to select and to test it, so it needs the independent replication noted above.

### Where the panel fails regardless: psoriatic arthritis

| test | AUC | p |
|---|---|---|
| PRJEB51315 AS vs PsA, blood | 0.43–0.67 | 0.087–0.700 |
| PRJEB51315 AS vs PsA, synovial fluid | 0.63–0.64 | 0.157–0.222 |

Even taking the leaky panel at face value, it never separated AS from psoriatic arthritis in either blood or joint. PsA is also chronic inflammatory arthritis, and the single known clonotype has the same weakness, sitting in 25% of PsA patients. With 16 PsA subjects in blood and 17 in synovial fluid from one deposit, the limit is clear in direction and poorly estimated in size — and it is the clinical blocker, since PsA is the differential that actually presents.

#### The lesson

Selection and evaluation have to be separated at *every* stage, not just the last one. Here the final scoring step was correctly fold-restricted, which made the procedure look cross-validated, while an earlier filter upstream had already seen all the labels. A permutation test does not catch this on its own: shuffling labels and re-running the *same leaky pipeline* reproduces the leak in the null only if the leaked stage is inside the shuffle. It was not, which is why the null looked clean at p < 0.01 and the honest answer is p = 0.17.

---

## V genes: recoverable after all, and the clone is TRBV9

Read-level extraction pattern-matches the junction and never compares the read to a reference, so it kept the CDR3 and discarded the germline V call. That is the wrong half for an HLA-driven disease. CDR3 contacts the **peptide**; CDR1 and CDR2 come entirely from the V segment and contact the **MHC molecule**. HLA-B27 is MHC class I, and this clone sits 35× higher in CD8 than CD4 cells — i.e. it is class-I restricted — so the V segment is the part doing the B27 reading.

The V sequence was never missing from the data: 95% of junction-spanning reads carry ≥60 nt of V-segment sequence immediately 5′ of the conserved cysteine. [`scripts/vgene_pilot.py`](scripts/vgene_pilot.py) assigns V by comparing that region against every germline TRBV allele, using OLGA's IMGT-derived reference (GPL). No aligner is needed, since V segments are not spliced here. MiXCR would be more thorough but requires a **paid licence for commercial use**.

**`CASSVGLYSTDTQYF` is TRBV9** — in all three AS patients tested (P18 ×23, P22 ×35, P19 ×13; 71 reads, none assigned elsewhere), 35 of 39 at distance 0 with a median margin of 3 over the runner-up gene. That reproduces the published TRBV9 association in AS from data previously thought untestable, and it is a fourth independent line pointing at one mechanism: CD8 → class I → B27 is class I → the V segment reads class I → the V segment is TRBV9. None of those is a thing a denominator bug or batch effect produces.

Three things a fresh session needs to know before using this:

- **TRBV9 is common and is not a marker.** It is 3.61% of V-assigned reads in an AS patient and 2.53% in a healthy control — roughly 1 in 30 of anyone's T cells, 7th–10th most used gene in both. What is rare is the *pair* (TRBV9 + that junction), not either part. V-gene usage on its own should be expected to discriminate weakly.
- **Assignment needs a margin rule, not a distance cutoff.** Germline probes are unevenly spaced: TRBV9 sits 3 mismatches from TRBV13 and TRBV7-4, inside an absolute cutoff of 4. Requiring the best gene to beat the runner-up by ≥2 refuses a third of reads (90% → 53% of junctions) and visibly shifts the usage profile — TRBV20-1 moves 18.2% → 29.8% — because genes with close neighbours are refused more often. Irrelevant for naming one clone's V gene; essential for comparing profiles.
- **Read length is a group confound.** 24 of 27 AS libraries are 2×100 bp against 8 of 12 healthy, so 11% of AS libraries are long-read against 33% of healthy (p = 0.091). Longer reads raise assignment rates gene-specifically, so `--truncate 100` is mandatory for any cross-group usage comparison. Probe length is capped by read length: at 150 bp the cysteine sits a median 78 nt in and a 45 nt probe works; at 100 bp the median is 31 nt and only 0.2% of reads clear 45, which is why the probe is 25 nt.

### The representation question this reopens

Two earlier failures share a cause. Exact CDR3 matching is too strict — 5.3M distinct strings, with near-identical receptors counted as unrelated. Edit-distance widening is too loose — 74% of healthy repertoires contain something within 1 edit of the clone, taking p from 0.0033 to 0.58. **Both ignored the V gene.** "Within 1 edit" sweeps in receptors built on different V segments that engage different MHC molecules; "within 1 edit *and* TRBV9" requires the same MHC-binding scaffold as well as a similar peptide-contact loop.

| representation | size | status |
|---|---|---|
| exact CDR3 strings | 5.3M | too large for 36 patients; produced a leaked AUC of 0.96 |
| V gene alone | 59 | too coarse — everyone carries every V gene |
| V/J usage levels | 624 | testable, expect weak signal |
| **(V gene + similar CDR3) families** | not yet computed | the biologically motivated unit; untested |

The promising next step is re-running the similarity search *conditioned on V gene*, which repairs an already-diagnosed failure and needs no million-feature selection step to leak through. Testing V/J usage levels is the weaker hypothesis.

### A third representation, with wet-lab ground truth: functional distance

[Wang, Yeh & Garcia](https://doi.org/10.1038/s41587-026-03128-x) (*Nat Biotechnol* 2026) profiled 21 HLA-B\*27:05-restricted TCRs from AS and acute anterior uveitis patients against a ~10⁹-variant 9-mer peptide library. **Their panel is the Faham `S··YSTDTQYF` family, measured** — it contains our clone twice (`4.1`, `19.2C5`) plus four of its one-edit neighbours. A *peptide recognition profile* gives the functional similarity metric the two failures above lacked: receptors are alike if they bind the same peptides, whatever they spell. On their cached profiles, [`scripts/prp_divergence.py`](scripts/prp_divergence.py) reproduces their published result — **motif Jensen-Shannon distance is uncorrelated with CDR3β edit distance** (ρ = +0.113, *p* = 0.10, flat from 0 to 7 edits), with their engineered 19.2 neighbourhood recovered as a positive control (JS 0.158 vs an all-pairs median of 0.224). The clone's nearest functional relatives sit **three** edits away, so radius-1 widening was searching the wrong neighbourhood, not too wide a one.

It does not transfer to disease. Scanning all 16 panel sequences against 31 AS, 10 healthy and 16 PsA subjects ([`scripts/prp_candidates_scan.py`](scripts/prp_candidates_scan.py); unsorted blood, uniform ≥2-read floor, carriage at ≥10⁻⁵) **inverts the expectation**: AS carriage tracks **edit** distance (ρ = −0.644, *p* = 0.007) and not functional distance (ρ = −0.028, *p* = 0.92). Both closest functional relatives fail — `CASSLEVYSTDTQYF` is 1/31 AS and 1/10 healthy, `CASSLGLYSTMEQYF` absent from all 57 subjects — and the reason **independently supports the V-gene diagnosis above**: `CASSLGLYSTMEQYF` is the only panel sequence with a different J region, `TMEQYF` (TRBJ2-5) against every other's `TDTQYF` (TRBJ2-3). It may well bind the same peptides; it is a different rearrangement, so it is not in these repertoires. **Functional similarity and disease association are different questions** — a clone must also be generatable and selected — and the signal stays with the TRBV9/TRBJ2-3 scaffold. Read the negative weakly: functional distance exists for **16 sequences**, 9 with zero AS carriers, so it cannot separate a weak real relationship from none. Two by-products worth keeping: their released checkpoints **cannot** score novel CDR3βs — single-TCR models held β constant while fine-tuning, and the `19.2_allmutants` joint model learned β as a near-pure intercept (pairwise *r* = 0.9996 between β logit vectors, all 14 held-out βs at AUC 0.784–0.789, so [`scripts/prp_score.py`](scripts/prp_score.py) records a failed control rather than a result, having nearly produced a convincing ranking with our clone first); and the cached profiles sit inside 3.2 GB archives that are 99.96% model weights, so range-fetching one zip member cost **32 MB instead of 70.6 GB**.

**Candidate worth revisiting:** `CASSVATYSTDTQYF` (2 edits, motif JS 0.237) is **4/31 AS, 0/10 healthy, 0/16 PsA** at ≥10⁻⁵ and present in GSE288581's AS cells. It is in the wider 1,014-sequence panel but not the 167, so the burden work surfaced it independently at ≥4 carriers. *p* = 0.27 against PsA across 16 tests, so a candidate, not a finding; the untested arms are the 64 synovial-fluid runs and Emerson for specificity.

---

## Data-handling rules, earned the hard way

Every error in this project so far belongs to one of two families. Check both before believing any result.

**Denominators.** Compute a numerator and its denominator in the *same code path*. A 34× frequency inflation survived review because the motif was counted on both strands and its denominator on one. Normalise depth before any diversity, clonality or presence comparison, and thin to a depth near the *shallowest* library, not a round number far below it. Never sum across sorted fractions (CD4/CD8/unsorted) when groups have different fraction mixes — it mixes denominators, and here it diluted a real effect.

**Selection leakage.** Separate selection from evaluation at *every* stage, not just the last one. A fold-restricted final step makes a pipeline look cross-validated while an upstream filter has already seen all the labels; that alone turned a null result (p = 0.17) into AUC 0.96. A permutation test does not catch this unless the leaked stage sits inside the shuffle. Any threshold chosen by looking at results is post-hoc — report the threshold-free test alongside it.

**Two control arms, always.** An external cohort kills public clonotypes; the study's own controls kill batch artefacts. Neither does the other's job: `CRVSGLWAKNIQYF` is absent from all 865 external repertoires (p = 2×10⁻⁵) and carried by 3 of 10 controls in its own deposit.

**Every sweep reports a positive control.** A screen that cannot recover a signal known to be present is uninformative, not negative.

**Persist the artefact.** The 167-sequence panel existed only in memory for several sessions, which made every number computed from it unreproducible. Write it to the repo, then score from the file.

**Beware the glob.** Scan JSONs were duplicated across directories, so a recursive glob counted 153 runs as 306. Ratios survived; the detection floors derived from them were 2× optimistic. Count distinct accessions before aggregating.

---

## Errors that each nearly became a finding

### A reverse-strand bug inflated every frequency by ~34×

The read-level scan searched **both** strands for the target motif but only the **forward** strand for the `CASS` control used as the denominator. In these 5′RACE libraries, read 2 carries the CDR3 reverse-complemented: 47.2% of R2 reads contain `CASS` on the reverse complement and 0.06% on the forward strand. Measured across both mates, the scan counted `CASS` in 0.72% of reads where the truth is 24.3%.

So `motif_reads / control_reads` overstated frequency ~34× — and not by a constant, since the inflation tracks each library's R1/R2 mix, which varies 7–8× across subjects. Confirmed independently against the CDR3 tables, where 58.6% of junction reads carry a CDR3 beginning `CASS`, as expected biologically.

**What it cost.** The previous headline — "AS 9/25 expanded ≥10⁻⁴ versus healthy 0/785" — compared inflated amplicon frequencies against Emerson frequencies computed correctly from Adaptive's own clone tables. Exactly the asymmetry that manufactures a positive result. On the honest denominator **no AS patient in either cohort reaches 10⁻⁴**, and that claim is withdrawn.

**What survived.** Presence was never affected — the motif search did try the reverse complement — so every presence number in this README stands. The result is now about presence and relative magnitude rather than clearing a threshold, which is a weaker claim but does not depend on a cutoff chosen after the fact. Threshold-free Mann-Whitney on frequency: p = 0.0051 (Komech AS vs healthy), p = 0.013 (AS vs PsA blood).

The general lesson: **a numerator and its denominator must be computed by the same code path.** They were not, and nothing about the output looked wrong.

### Thinning to a round number destroyed the sweep's sensitivity

The first sweep thinned every subject to 30,000 junction reads to equalise depth. The libraries hold ~5 million. At 30K reads only clones above ~10⁻⁴ are visible, so that sweep could not have found our own clonotype, and its null result meant nothing. Thinning is still right — AS libraries are deeper than healthy ones, and under label shuffling that difference would manufacture enrichment — but the target depth belongs near the **shallowest** library (~800K here), not far below it. Every sweep now reports a positive control, and a sweep that cannot recover it is reported as uninformative rather than negative.

### A permutation null that could not be computed

The sweep is a Fisher test per clonotype under 200 label shuffles: hundreds of millions of `fisher_exact` calls, which does not finish. Two facts make it cheap. With group sizes fixed, the one-sided p-value depends only on the carrier count among cases and the total carrier count — so one small precomputed table holds every p-value any shuffle can need. And the per-clonotype case count is a single sparse matrix-vector product. Each shuffle becomes a matvec over ~7M non-zeros plus integer lookups: milliseconds instead of hours. The null is cheap and therefore non-negotiable.

### Earlier near-misses, kept short

- **A `CASS`-fraction denominator audit that almost sank a sound result.** Across 39 amplicon libraries `CASS` spans 0.022%–0.966%, a 44-fold range, which looked like a broken denominator. It was: see the strand bug above. The audit found the symptom and the wrong cause.
- **Bulk RNA-seq cannot test this at all.** GSE221786 gives a detection floor of 3.6–3.8×10⁻⁴ per sample — about 4× too coarse, needing ~190M reads/sample to fix. An early pilot also reported 0.50% of reads as `CASS`, which was pure chance: `CASS` is 12 nt, matching ~1 position in 58,000, and 139 positions × 2 strands predicts 0.476%. Observed: 0.5055%. The assay, not the analysis, is the limit.
- **Widening the search by sequence similarity fails twice.** Edit distance saturates — 74% of healthy repertoires contain something within one edit, taking p from 0.0033 to 0.58. And ESM-2 embeddings return anagrams as nearest neighbours, because mean-pooling discards order.
- **Batch effects here are enormous.** Cohort membership is predictable from repertoires at **AUC 0.994**, and HLA-A\*02 at 0.953. Any cross-deposit comparison has to be read against that. For calibration: CMV serostatus transfers across cohorts at 0.851, while age+sex reaches only 0.584.

---

## The earlier arms, in brief

### Bulk transcriptomics: NK cells are quiet, not absent

Five cohorts, **97 AS patients vs 70 healthy controls**, random-effects meta-analysis, then deconvolution against a 312,000-cell reference to separate cell proportion from per-cell expression.

| | pooled *g* | 95% CI |
|---|---|---|
| NK cytotoxic cargo (granzymes, perforin, granulysin) | **−0.98** | −1.33 to −0.63 |
| …after adjusting for NK cell proportion | **−0.80** | −1.17 to −0.44 |
| NK cell proportion | −0.54 | −0.91 to −0.17 |

**82% of the cytotoxic deficit survives adjustment for cell number.** The two numbers use **disjoint gene sets** — proportion from NK identity genes with every effector gene deleted from the signature matrix, outcome from the effector genes — because otherwise the analysis is circular and returns "it's all composition" regardless of truth. Count the soldiers by their uniforms, then ask how much ammunition they carry.

A [meta-analysis of 95 flow-cytometry studies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8291033/) (4,020 patients, 3,065 controls) found **no difference in NK percentage**. Flow counts cells by surface protein and is blind to transcription; our method infers counts from transcripts and mistakes a quiet cell for an absent one. Opposite blind spots, same conclusion.

Not a diagnostic: AUC 0.77, and in a back-pain clinic a positive result means a **17.7% chance of AS**. CRP manages 15.6%, HLA-B27 37.2%. All three are rule-out tests. And every control is a healthy volunteer — the comparison that matters clinically is AS versus *mechanical* back pain, and no public dataset has it.

Also there: platelet and integrin genes are the strongest upward signal (10 of 11 markers up, mean *g* = +0.44), matching a [2025 plasma proteomics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12957791/) that independently found integrin-pathway proteins as its top AS biomarkers. Never followed up.

Validation that changed the answer: deconvolution was tested on held-out pseudobulk from 13 unseen patients. NK recovered at *r* = 0.83; **CD8 T failed at 0.57 and was discarded** despite a nominally significant result (p = 0.029), because cytotoxic T cells share GNLY, NKG7 and PRF1 with NK cells. The genome-wide median effect is −0.13, not zero, so against a 50/50 null every downward gene set looks significant.

### Cell-level model: a receptor explains 1–2% of its own cell's state

**1,862,114 cells** across 291 libraries with both a productive TRB and a transcriptome; 1,421,133 distinct clonotypes, **91.3% singletons**. Stage 1 compresses 2,009 genes to a latent using expression alone; stage 2 predicts that frozen latent from CDR3 + V/J. Every measurement runs twice, once with receptors permuted within subject — **the gap is the signal**.

| | real | shuffled |
|---|---|---|
| exhaustion, per cell | +0.247 | −0.003 |
| exhaustion, **per clonotype** | **+0.126** | +0.010 |
| exhaustion, held-out cohort | +0.012 | +0.005 |

Half the per-cell signal was pseudo-replication: cells in one expanded clone share a receptor *and* a history. Counting each clone once took 0.247 → 0.126.

The pipeline is sound — positive controls recover known receptor→phenotype links at literature magnitudes (CD4-vs-CD8 0.658, MAIT 0.802 on the AS cohort). Two findings matter more than the model:

- **Plain features match the transformer.** Logistic regression on V/J one-hot, CDR3 length and amino-acid composition ties a 1.5M-parameter transformer to within 0.02 AUC and beats it on one control. The information in play is gene usage and composition, not sequence order.
- **Nothing transfers across cohorts, not even known biology.** CD4-vs-CD8 falls 0.658 → 0.528 on a held-out study; MAIT 0.612 → 0.507. Both arms, both feature sets. That is a larger obstacle than the weak programme signal and plausibly its cause.

A systematic inventory blind spot worth knowing: **10x 5′ libraries with V(D)J enrichment are indexed as "single-cell RNA-seq."** They are deep TCR sequencing with a paired transcriptome, and no repertoire inventory surfaces them. Of 450 series checked, 86 have both — 6,008 samples. An AIRR Data Commons sweep returned zero AS repertoires and we wrongly concluded none existed; the archives hold them filed as ordinary sequencing.

---

## What is left in the data inventory

The shortlist holds 119 candidate studies. For AS specifically it is exhausted: `PRJNA853355` is a single CD4 sample, `E-MTAB-11498` is the ArrayExpress mirror of `PRJEB51315` (same 115 runs), `GSE207037`/`GSE290921` are CD4-only, and `E-MTAB-10948` has no VDJ. Worth pursuing, in order:

1. **HLA-B27-associated acute anterior uveitis** — three untouched datasets (`GSE274471`, 88 samples of aqueous humour + blood scTCR; `GSE232802`; `PRJNA973988`). Uveitis is on the spondyloarthritis spectrum and shares the B27 association. If the clone is B27-driven rather than AS-driven it should appear in B27+ uveitis patients too — which would reframe it as a B27 marker and largely kill its diagnostic value over a B27 test. This is the single most informative experiment left.
2. **Steel 2019** (immuneACCESS `/pub/steel-2019-ar`) — IL-17A+CD8+ cells in PsA, assay-matched to Emerson. The most on-target control available for the PsA arm, which currently rests on 16 subjects. Needs a user-performed export.
3. **HLA-B27 status for the Komech patients** is absent from all metadata. It gates the only comparison that matters diagnostically, and the 45 Emerson B27+ controls are the sole substitute. Worth emailing the authors.
4. **TRBV9 usage is untestable from our tables.** Read-level extraction recovers the junction but not the germline V call, so the published TRBV9 association cannot be checked without alignment — a separate job.
5. Longitudinal sampling exists in `PRJEB51315` but almost entirely in PsA patients (5, plus AS patient `GE` at two timepoints), too thin to test clonal persistence in AS.

---

## Reports

| | |
|---|---|
| [BULK_META.md](BULK_META.md) | The meta-analysis and deconvolution: datasets, methods, every caveat. |
| [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md) | Why compartment beats modality, and the repertoire inventory in full. |
| [NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) | Handoff: the patient-level plan and the traps that cost time. |
| [REPORT.md](REPORT.md) | Earlier single-cell work — GSE277117 (Song et al., *Nat Commun* 16:5531, 2025), TNF-inhibitor response in 25 patients. |
| [NEXT_STEPS.md](NEXT_STEPS.md) | Project plan and the running list of traps. |

## Running it

```bash
python -m venv .venv && .venv/Scripts/pip install -r requirements.txt
```

The bulk meta-analysis runs locally in minutes. Repertoire extraction and the single-cell reference need [Modal](https://modal.com). Every script re-fetches its own data; nothing raw is committed.

```bash
# extract CDR3s from raw reads (Modal, caches per run on the volume)
modal run --detach scripts/komech2018.py::cdr3_all
modal run --detach scripts/as_trust4.py::cdr3_all

# open-ended sweep, with a positive control
python scripts/clono_sweep.py --tables <dir> --case as --ctrl healthy \
    --depth 800000 --min-case 4 --control-cdr3 CASSVGLYSTDTQYF

# rebuild the burden panel; assign V genes (both local, no Modal)
python scripts/build_panel.py --tables <dir>
python scripts/vgene_pilot.py --runs SRR5812652 SRR5812641 --truncate 100
```

Provenance of numbers quoted above, for anything not obvious from the script name:

| script | what it produced |
|---|---|
| [`bulk_score.py`](scripts/bulk_score.py) | the NK score's AUC 0.77 and the back-pain-clinic PPVs |
| [`rnaseq_pilot.py`](scripts/rnaseq_pilot.py) | the bulk RNA-seq ruling-out (detection floor 3.6–3.8×10⁻⁴) |
| [`esm_embed.py`](scripts/esm_embed.py) | the ESM-2 negative result (nearest neighbours are anagrams) |
| [`as_clonotypes.py`](scripts/as_clonotypes.py) | `data/as_vdj/clonotypes.parquet`, the GSE288581 sorted-CD8 numbers |

### Never let a last-step failure discard a finished run

This cost us three times, most expensively when a training run completed all 12 epochs then died computing an optional diagnostic, throwing away the weights.

**`--detach` does not help, and it is worth being clear why.** It detaches the app from the local client, so the run survives a closed terminal or a lost network — the run that threw away 12 epochs *was* detached. It does not checkpoint, persist partial results, or catch exceptions.

What works, in order of value:

1. **Persist the expensive artefact before anything optional.** Save and `vol.commit()` immediately after training, then compute diagnostics.
2. **Wrap every optional diagnostic in `try/except`** and record the exception into the results dict. A failed metric should appear as `"library_probe_error": "TypeError: ..."`, not as a dead run.
3. **Cache per unit of work.** Every extraction function checks for its own output file first, so a rerun recomputes only what is missing. This turned a stalled 30-of-39 run into a 9-library job.
4. **Use `.map(..., return_exceptions=True)`.** A container timeout or OOM is killed *outside* the interpreter, so `try/except` inside the worker cannot catch it and it propagates through `.map`, killing the dispatch.
5. **Smoke-test the whole path with one epoch.** The bug that killed stage 1 was reachable in 25 seconds and found after 5 minutes of GPU time.
6. **Pin library versions in the image.** One failure was `LogisticRegression(multi_class=...)`, removed in recent scikit-learn.

`scripts/vdj_expr_vae.py` is the worked example of all of these. Modal also preempts containers (`Container terminated due to preemption`), which costs wall time, not correctness.
