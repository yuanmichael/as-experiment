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

### Combining clonotypes works where single clonotypes fail

Thresholding one clonotype at a time is the wrong test, and arbitrary: "≥6 carriers and 0 controls" is a cut with no principle behind it. The alternative is a burden score — count how many candidate clonotypes a subject carries — evaluated leave-one-subject-out, selecting candidates from only the other 34 subjects.

| selection control arm | best LOO AUC |
|---|---|
| Komech's 10 healthy controls | **0.44** (at or below chance) |
| **Emerson's 786 healthy controls** | **0.962** |
| single best clonotype, presence alone | 0.720 |

Selecting on "absent from 9 training controls" picks public clonotypes that happen to be missing from nine people, and held-out healthy subjects carry them just as often — so that score is worthless. Swap in 786 controls and the same score reaches **AUC 0.962**, well above the 0.720 of the best single clonotype.

Permutation-tested with the nine-configuration grid search *inside* the null, so the p-value pays for the search: observed best-of-grid 0.962 against a null mean of 0.626 and a maximum of 0.838 over 100 shuffles, **p < 0.01**.

#### It is not the batch effect

Cohort membership is predictable from repertoires at AUC 0.994, so the first question is whether the score reads library prep rather than disease. Three checks say it does not:

- The comparison is **within one deposit** — Komech AS against Komech's own healthy controls, same lab and protocol — the opposite setup from the cross-deposit 0.994 result.
- Library prep batches are **interleaved, not disease-segregated**: 7 of Komech's batches (`MK-108/135/142/181/218/227/248`) contain both AS patients and healthy controls. Restricted to those — 8 AS, 7 healthy, batch held constant — **AUC = 0.991**. Fixing batch does not degrade the score.
- Within the healthy controls the score barely moves across sequencing lanes (medians 0.036–0.066), and 22 of 25 AS patients score above the highest healthy control (0.090).

#### Operating points, and why PPV is not yet meaningful

| threshold | sensitivity | specificity | TP/FN | FP/TN |
|---|---|---|---|---|
| 0.07 | 92.0% (74–99%) | 80.0% (44–97%) | 23/2 | 2/8 |
| 0.09 | 92.0% (74–99%) | 90.0% (55–100%) | 23/2 | 1/9 |
| **0.10** | **84.0% (64–95%)** | **100.0% (69–100%)** | 21/4 | 0/10 |
| 0.15 | 76.0% (55–91%) | 100.0% (69–100%) | 19/6 | 0/10 |
| 0.20 | 56.0% (35–76%) | 100.0% (69–100%) | 14/11 | 0/10 |

At threshold 0.10 the panel reaches **84% sensitivity at 100% specificity**, against the single clonotype's 56% at 98%. The gain is in sensitivity, which was the single clonotype's fatal weakness.

**But PPV cannot be estimated from 10 controls.** At threshold 0.10 there are zero false positives, so PPV is 100% at every prevalence — meaningless. Drop to 0.09, one false positive, and PPV at 20% prevalence falls from 100% to **69.7%**. Specificity moves only in 10-point steps here and its confidence interval runs 69–100%. Any PPV from this is a statement about ten people, not about a test. The thresholds are also chosen on the same data that evaluates them.

#### It transfers to an independent healthy cohort

The panel was built with no Britanova data. Scored on [Britanova's](https://doi.org/10.5281/zenodo.826447) 80 healthy donors (clonotype tables, same amplicon platform), thinned to the same 800K depth:

| group | n | median score | range |
|---|---|---|---|
| Komech AS | 25 | **0.228** | 0.060–0.335 |
| Komech healthy | 10 | 0.054 | 0.024–0.090 |
| **Britanova healthy** | **41** | **0.036** | 0.000–0.096 |

| comparison | AUC | p |
|---|---|---|
| AS vs Komech healthy (10) | 0.962 | 1.3×10⁻⁵ |
| **AS vs Britanova healthy (41)** | **0.991** | **1.4×10⁻¹¹** |
| AS vs pooled 51 controls | 0.985 | 3.8×10⁻¹² |

Independent healthy donors score like healthy controls, not like AS patients. Operating points against all 51 controls:

| threshold | sensitivity | specificity | TP/FN | FP/TN |
|---|---|---|---|---|
| 0.09 | 92.0% (74–99%) | 96.1% (87–100%) | 23/2 | 2/49 |
| **0.10** | **84.0% (64–95%)** | **100.0% (93–100%)** | 21/4 | **0/51** |
| 0.15 | 76.0% (55–91%) | 100.0% (93–100%) | 19/6 | 0/51 |

NPV 96.2% at 20% prevalence. PPV remains 100% at every prevalence because there are no false positives in 51 controls — read the specificity interval instead: the true value could be as low as 93%, which would put PPV near 75% at 20% prevalence.

Caveats: Britanova is Chudakov's group, as Komech is, so this controls for cohort and population but not fully for protocol. The 0.10 threshold was chosen by looking at Komech's controls, then applied here unchanged. Only 41 of 80 donors reached 800K depth; at 200K (73 donors) healthy median is 0.006, the same conclusion at a lower scale.

#### What it does not do: separate AS from psoriatic arthritis

Scoring a *different* deposit's patients on the same 167-sequence panel:

| test | AUC | p |
|---|---|---|
| AS vs healthy, pooled 51 controls | **0.985** | 3.8×10⁻¹² |
| PRJEB51315 AS vs PsA, blood | 0.59–0.67 | 0.087–0.299 |
| PRJEB51315 AS vs PsA, synovial fluid | 0.63–0.64 | 0.157–0.222 |

This was the open question and it is now resolved in one direction only. Two readings were possible: the panel might be **cohort-specific**, or **AS-versus-PsA** might simply be a harder contrast than AS-versus-healthy. Britanova settles it — an independent healthy cohort scores at 0.991, so the panel is not a cohort fingerprint, and the PsA failure is a genuine specificity limit rather than evidence against the panel.

That limit is the clinical blocker. A test separating AS from healthy at 0.99 and from PsA at 0.63 is not usable in a rheumatology clinic, where PsA is the differential that actually presents. The single known clonotype has the same weakness, sitting in 25% of PsA patients. PsA controls currently number 16 in blood and 17 in synovial fluid, from one deposit, so the limit is firmly established in direction and poorly estimated in size.

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
```

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
