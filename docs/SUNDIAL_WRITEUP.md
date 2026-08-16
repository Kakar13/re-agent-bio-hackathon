# Human HLA binding and presentation profiling for de novo protein design

re:AGENT hackathon, August 15–16 2026. Track C: biological design to a spec, plus
evidence it could hold up.

## The problem

A standard campaign — RFdiffusion, ProteinMPNN, then AlphaFold — can leave
thousands of structurally plausible candidates. Folding and target contact do
not say whether the protein will later be processed into peptides that human
MHC molecules present. Those presentation events are one framework through
which immunogenicity can begin: cytosolic MHC-I processing and HLA-A*02:01
presentation on the nucleated cell that made the protein, including immune
cells. TCR recognition, danger sensing, and animal or human outcomes are
later steps. They are why a construct that must function in a human still
has to be tested in mice. They are not what this model scores.

This MHC-I-only filter is the primary match when the designed protein is
expressed inside the cell from a plasmid, mRNA, or viral vector. That is
not limited to secreted drugs. An RFdiffusion binder can be an
intracellular tool — a biosensor, or a binder to a phosphorylation site —
in which case the payload is delivered so the cell manufactures the
protein. Extracellular protein can still reach MHC-I by cross-presentation,
so an injected binder is not MHC-I-invisible. For that product class,
MHC-II, CD4 help, and antibodies remain the larger unmodeled gap.

The job is downselection. Remove the bulk of campaign-scale HLA-A*02:01
processing and presentation liabilities so wet-lab and animal capacity is
spent on a smaller, inspectable set. A flagged window is not a clinical
immunogenicity diagnosis. A clean window is not deimmunization.

We built a per-residue MHC-I processing profile fast enough to score every
9-mer of every design, with explicit separation between cleavage, TAP, BA
binding, EL presentation, and downstream recognition that we do not model.

The training corpus is the Protein Design Archive (PDA): a recent curated
collection of synthetic designed proteins, mostly from the last several
years. It is not the Protein Data Bank. The Chronowska et al. paper title
refers to 40 years of the protein-design field, not to the age of the
archive.

## System

| Stage | Tool | What it produced here |
| --- | --- | --- |
| Backbone generation | Proto `rfdiffusion3-design` on Modal | 500 binder backbones across 5 targets |
| Literature grounding | Paperclip | target hotspots, benchmark selection, teacher choice |
| HLA profiling | our trained model | 289,335-row MHC-I profile over the PDA corpus |
| Evaluation | held-out CV + an external teacher-fidelity cohort | evaluations below |

The design lane is reported in `results/rfd3_binders/report.md`. The PD-L1
positive control is a geometric check on RFdiffusion3 complexes, not a
refolded AlphaFold or ESMFold model: 15 of the first 16 backbones contacted
all three foundry hotspots (Y56, M115, Y123) within 4.5 Å, and all 16 did so
within 6 Å. The one miss, `pdl1_bb_0015.pdb`, is 5.87 Å from Y56 and still
contacts the other two hotspots inside 4.5 Å. ProteinMPNN sequence design and
AlphaFold/ESMFold validation were scoped out of the weekend. Everything below
concerns the HLA-A*02:01 binding and presentation lane, which is where our own
modeling work went.

## How to read the benchmark

This is an evidence ladder, not one leaderboard. Each comparison answers a different
question:

1. **Held-out teacher imitation:** can the student reproduce NetMHCpan on peptide
   components withheld from training? This measures distillation fidelity, not
   biological accuracy.
2. **Legacy-head replacement:** on identical PDA peptides, does the student surface
   rare NetMHCpan strong binders better than the chao1 MHCflurry head? Average
   precision is primary because only 2.5% of peptides are positive.
3. **External generalization:** on 95 unrelated NY-ESO-1 designs, does fidelity survive
   a distribution shift? This remains a teacher-fidelity test.
4. **Screening calibration:** given a good ranking model, which threshold matches the
   asymmetric cost of design-time screening? Missing a risky window is costlier than
   resampling a false flag.
5. **Operational and design value:** latency tests whether the model is cheap enough
   for an inner loop, while the structure join tests whether flagged windows are
   actually editable.

The 160,012-row cross-validation result comes directly from the training corpus. The
159,038-row head comparison uses the unique peptides that also materialized in the
occurrence-level profile. Those denominators are kept separate throughout. Likewise,
“weak-binder AUPRC” means teacher EL rank ≤ 2.0, while “strong-binder AP” means rank
≤ 0.5; reporting either without its threshold would mix two different tasks.

## Why we replaced the core MHC model

We started from `chao1`, a frozen three-head MHC-I processing surrogate: N/C-terminal
cleavage, TAP transport, and an MHC presentation head calibrated against MHCflurry.
The cleavage and TAP heads are independently supervised and we kept them. The MHC head
we did not keep, and the decision was evidence-driven rather than aesthetic.

Two problems. First, its training corpus had been sampled by MHCflurry quantile, so the
peptide distribution was shaped by the very model being imitated. Second, MHCflurry and
NetMHCpan disagree, and NetMHCpan 4.1 is the stronger presentation predictor in the
published comparisons we surveyed with Paperclip.

Rather than patch the existing corpus, we rebuilt it from raw parent proteins under a
teacher we could defend, and measured the difference.

The deeper rationale is a missing-label problem. NetMHCpan 4.1 is the field-standard
HLA-A\*02:01 model, trained mainly on natural peptide binding and eluted-ligand
evidence. Measured MHC data on de novo designed proteins are almost unavailable,
with only sparse exceptions, so we cannot supervise a campaign filter on
experimental designed-protein HLA labels. Designed 9-mers are out of distribution
for that original natural-peptide teacher corpus. We generalize NetMHCpan for de
novo use by a pseudo-de novo distillation: label every PDA 9-mer with the teacher,
embed it with frozen ESM-2, and train a compact student to reproduce those scores
in designed-protein representation space. After that training, PDA is not an
out-of-distribution test for the student. Held-out PDA folds measure teacher
imitation on designed peptides. A newly generated non-PDA sequence is the
remaining shift, and the labels themselves are still teacher scores, not measured
de novo HLA outcomes.

## Methods

### Corpus

Every parent in the Protein Design Archive de novo set — 1,602 recent synthetic
proteins, not PDB naturals — was tiled into all 9-mer windows with 4-residue
flanks and labelled by NetMHCpan 4.1 through the IEDB Tools API, on both the EL
(presentation) and BA (binding affinity) channels, for HLA-A\*02:01. Labelling
happened *before* any sampling, so no model's opinion selected the training
distribution. 290,937 windows collapsed to 160,012 unique peptides.

Parents sharing any exact 9-mer were joined into connected components (608 of them),
and whole components were assigned to folds. This matters more than it sounds: de novo
designs reuse motifs heavily, and per-peptide splitting would have leaked 40,125 shared
peptides across the train/test boundary.

### Student

A small head over frozen ESM-2 `esm2_t33_650M_UR50D` embeddings, mean-pooled across
the nine peptide residues, predicting three channels:

| Channel | Target | Why this target |
| --- | --- | --- |
| EL presentation | monotonic transform of percentile rank | raw EL scores pile up at zero and give the optimizer almost no gradient where the ranking matters |
| BA binding | monotonic transform of percentile rank | same reasoning as EL |
| BA affinity | NetMHCpan's `1 - log(IC50)/log(50000)` | exactly invertible to nanomolar, and only 0.02% of the corpus is censored at the 50,000 nM ceiling, so the raw scale is well behaved |

The affinity channel is what lets the profile report "predicted 180 nM" instead of only
"predicted rank 0.4%". It was added after the rank-only version, once we checked that
the affinity target was already well conditioned rather than assuming it shared EL's
pile-up problem.

Three tasks do not fit in the rank-only model's 256-unit hidden layer. At that width the
affinity channel costs 0.014 Spearman on EL; at 512 units it costs 0.005, which is the
deployed configuration:

| Model | EL Spearman | EL AUPRC | IC50 Spearman | Median fold error on binders |
| --- | --- | --- | --- | --- |
| 256 units, rank only | 0.948 | 0.901 | — | — |
| 256 units, with affinity | 0.934 | 0.874 | 0.930 | 2.47× |
| **512 units, with affinity** | **0.943** | **0.892** | **0.938** | **2.32×** |

Five heads were trained, each holding out one fold for test and the next for
validation. Training ran on Modal A10G and takes about three minutes once the
embeddings are cached, which is why re-deriving the width was cheap.

### Profile

The deployed profile keeps each axis in its own lane rather than collapsing to one
number:

| Axis | Source | Field |
| --- | --- | --- |
| Cleavage | chao1 | N- and C-terminal probabilities |
| Transport | chao1 | TAP score and bootstrap uncertainty |
| Binding | NetMHCpan student | BA propensity and predicted IC50 in nM |
| Presentation | NetMHCpan student | EL propensity |
| Binder class | derived | strong / weak / nonbinder from predicted EL rank |
| Screening flag | derived | recall-calibrated risk flag; see Results §4 |
| Accessibility | deposited structure | window mean relative SASA and buried fraction |
| Immunogenicity | — | **not provided** |

Accessibility is a design-actionability axis, not an immunological one. MHC-I peptides
come from proteasomal degradation of unfolded protein, so burial does not gate
presentation. What burial decides is whether the designer can act on a flag: a risky
9-mer on a surface loop can be resampled freely, while one packed into the core cannot
be changed without risking the fold.

EL is deliberately excluded from the inspectable composite processing score, which
combines N-cleavage, C-cleavage, and BA only. EL already encodes cleavage and
transport information, so averaging it with BA would count the same MHC event twice.
Affinity is also excluded: it is a monotone transform of the BA score
(`IC50 = 50000^(1 - score)`), not a third independent measurement.

Campaign ranking does not use that equal-weight composite. It uses a pathway score

`S = EL^0.70 * (sqrt(N * C))^0.30`

because peptide–MHC binding/display is the selective step (NetMHCpan EL is the
presentation endpoint) and proteasomal generation is a necessary but more
promiscuous soft gate. The 0.70/0.30 split is a biological prior, not a fitted
weight. Protein-level rank is the mean of the top-5 window `S` values.

The materialized dataset carries one row per 9-mer occurrence rather than per unique
peptide, because cleavage depends on flanking context: the same 9-mer in two different
parents has two different cleavage probabilities. Student values in that table are
out-of-fold — each peptide is scored by the one head that held its component out — so
the table contains no in-sample optimism.

## Results

### 1. Held-out teacher imitation

Pooled over all 160,012 out-of-fold predictions:

| Channel | Spearman vs rank | Weak-binder AUPRC (rank ≤ 2.0) |
| --- | --- | --- |
| EL, pooled out-of-fold | 0.943 | 0.892 |
| BA, pooled out-of-fold | 0.948 | 0.859 |

Held-out fold Spearman was 0.943 ± 0.004 for EL and 0.948 ± 0.004 for BA (population
standard deviation). The small fold-to-fold variance is the useful part:
component-grouped splitting did not produce one lucky fold.

The affinity channel recovers the teacher's nanomolar values closely across the corpus
and less tightly where it matters most:

| Measure | All 159,038 peptides | 5,976 binders at IC50 ≤ 500 nM |
| --- | --- | --- |
| Median fold error | 1.11× | 2.32× |
| Within 2× of teacher | 90.2% | 42.7% |
| Within 5× of teacher | — | 79.1% |
| Spearman vs teacher IC50 | 0.938 | — |

The overall figures are flattered by the 96% of peptides that are non-binders, where
both the teacher and the student agree on "very weak" and being off by 30% of a large
number is easy. The binder column is the honest one: for peptides that actually bind,
the student lands within about 2.3× of NetMHCpan's IC50, and within 5× four times out
of five. That is useful for triage and ordering candidates, and it is not a substitute
for running the teacher when a specific number matters.

### 2. Old head vs new head, on the same rows

Both scores evaluated against the same NetMHCpan EL teacher on 159,038 unique PDA
9-mers, of which 3,956 are strong binders. Neither model saw these labels for the rows
it is scored on.

| Score | Spearman vs teacher | Strong-binder AUROC | Strong-binder AP |
| --- | --- | --- | --- |
| chao1 MHCflurry head | 0.354 | 0.731 | 0.060 |
| NetMHCpan student, EL | **0.943** | **0.993** | **0.820** |
| NetMHCpan student, BA | 0.829 | 0.985 | 0.654 |
| NetMHCpan student, affinity | 0.823 | 0.985 | 0.651 |

The AUROC column understates the gap and the AP column tells the truth. Strong binders
are 2.5% of the corpus, so chao1's 0.731 AUROC translates to an average precision of
0.060 — for the actual task of surfacing the rare risky peptide, it is close to
unusable, and the student is roughly 14× better on the metric that matters under this
class imbalance. This is the evidence that justified replacing the head.

Three-way framing matters here, because two of these comparisons are not the same kind
of claim. Against **chao1** the student is straightforwardly better, and the table above
is the evidence. Against **NetMHCpan itself** the student cannot win and we do not claim
it does: it is distilled from that teacher, so 0.943 Spearman is a ceiling being
approached, not a bar being cleared. What the student offers over the teacher is not
accuracy but availability — a local forward pass with no API, no rate limit, and no
network round trip, which is what makes per-residue scanning affordable inside a
generation loop. Section 6 measures that at 13.6×.

### 3. External teacher-fidelity cohort: de novo NY-ESO-1 designs

Our external test contains 95 designed HLA-A\*02:01 9-mers from an unrelated de novo
campaign (Visani et al. 2025, PNAS). Those peptides were designed by HERMES against the
1G4 TCR and selected on TCRdock/AlphaFold3 PAE, a structural criterion with no
connection to our PDA training corpus. We called NetMHCpan fresh on all 95 and scored
every model on identical inputs.

**Out-of-corpus imitation (valid).** The deployed v4 student holds up away from its
training distribution: Spearman 0.909 (EL) and 0.898 (BA) against freshly called
NetMHCpan, and 97.9% agreement on the weak-binder call at rank ≤ 2.0. That is a modest
drop from the 0.943/0.948 measured in-corpus, which is the honest expected degradation
rather than a collapse.

The student is systematically conservative at the strong-binder threshold — it called
73 strong where the teacher called 88, and the same direction appears in the PDA
profile. For a risk-screening tool this is the wrong direction of error. Section 4
quantifies it properly and fixes it.

The cohort's T-cell activation labels are not used as a model benchmark because they
depend on downstream TCR recognition. Direct biological validation requires a
different cohort: quantitative human HLA-A*02:01 IC50 or \(K_d\) for BA, and
monoallelic human HLA-A*02:01 immunopeptidomics for EL and the processing stack.

### 4. The default cutoff is the wrong operating point

The 95-peptide external cohort hinted at a conservative bias. Measuring it on all
159,038 out-of-fold PDA peptides shows it is much larger than that hint suggested:
**the default strong-binder rule catches only 52.3% of the teacher's strong binders.**
It misses nearly half of them, at 87.1% precision on the ones it does flag.

The cause is a threshold choice, not a broken model. The class boundary inverts a
predicted percentile rank, and because the student's predictions are slightly regressed
toward the middle, a fixed rank cutoff systematically under-calls the tail. Ranking
quality is fine — AUROC is 0.993 — so the ordering is right and only the cut is wrong.

Choosing the cutoff by target recall instead fixes it without retraining:

| Operating point | Predicted-rank cutoff | Recall | Precision | Peptides flagged |
| --- | --- | --- | --- | --- |
| Default (rank ≤ 0.5) | 0.50 | 52.3% | 87.1% | 2,376 |
| 80% recall | 0.88 | 80.0% | 67.0% | 4,727 |
| 90% recall | 1.24 | 90.0% | 53.7% | 6,636 |
| **95% recall** | **1.65** | **95.0%** | **43.9%** | **8,563** |
| 99% recall | 3.03 | 99.0% | 27.4% | 14,277 |

The deployed profile ships the 95% row. The precision cost is real, but the asymmetry
favours recall heavily: a missed strong binder ships a predicted HLA liability into a
design, while a false flag costs one resampled window during generation. At that cutoff
the screen flags 5.4% of the corpus, which is a tractable review load, and it recovers
**1,688 teacher strong binders that the conventional rank rule misses**.

The profile reports both calls rather than replacing one with the other, because they
answer different questions:

| Field | Question | Rule | Recall | Precision |
| --- | --- | --- | --- | --- |
| `binder_class` | what would NetMHCpan call this? | predicted EL rank ≤ 0.5 | 52.3% | 87.1% |
| `screening_flag` | should a designer look at this? | EL propensity ≥ 0.7892 | 95.0% | 43.9% |

Keeping `binder_class` on NetMHCpan's conventional thresholds matters: "strong binder"
has a standard field meaning at rank ≤ 0.5, and quietly redefining it to hit a recall
target would make our outputs incomparable with everyone else's. The screen gets its own
field and its own documented cutoff, and the ranking model is untouched — this is a
threshold choice, not a retrain.

### 5. Are flagged HLA-binding windows reachable by redesign?

Joining per-residue relative solvent accessibility from the PDA parents' own deposited
structures onto the HLA profile answers the follow-up question a designer asks after
seeing a flag: can I do anything about it?

The answer is mostly unwelcome. Across 186,485 windows in 1,184 parents with a matched
deposited structure, high-risk windows are preferentially **buried**, not exposed:

| Student binder class | Windows | Mean relative SASA | Buried share |
| --- | ---: | ---: | ---: |
| strong | 2,652 | 0.171 | 80.2% |
| weak | 9,620 | 0.176 | 77.1% |
| nonbinder | 174,213 | 0.200 | 67.7% |

Of the 12,272 windows flagged strong or weak, only **22.2% are surface-exposed** and can
be handed to ProteinMPNN for free resampling. The remaining 77.8% are buried and need a
backbone-aware fix or a fold-stability check. Risk and exposure correlate negatively
overall (Spearman −0.130).

The effect is directionally clear and consistent across both classes, though modest in
size: a 12.5-point difference in buried share between strong binders and non-binders.
It is also biophysically unsurprising in hindsight — HLA-A\*02:01 prefers hydrophobic
anchor residues at P2 and P9, and hydrophobic residues are exactly the ones that pack
into a core. The practical consequence is sharp: the predicted binder windows most worth
reviewing sit in positions least tolerant of mutation, so an HLA screen that ignores
structure will keep proposing edits that break folds.

Coverage is a real limitation. 390 of 1,602 parents could not be matched to a deposited
chain sequence and 28 had no retrievable structure, so this covers 74% of the corpus.
Burial is measured on the deposited model, so residues buried at a crystallographic
interface count as buried, and windows with fewer than five resolved residues are
excluded rather than averaged over partial coverage. Full breakdown in
`results/benchmarks/epitope_accessibility/REPORT.md`.

### 6. What the student actually buys: scanning cost

The student cannot beat its teacher on accuracy, so the case for it rests entirely on
being cheap enough to run inside a generation loop. Scoring every 9-mer window of a
156-residue design against HLA-A\*02:01:

| Path | Wall clock | Windows/second |
| --- | ---: | ---: |
| Student, CPU, end to end | 0.76s | 193 |
| NetMHCpan via IEDB API | 10.34s | 14 |

**13.6× faster**, and the qualitative difference matters more than the factor: scanning a
full protein drops from ten seconds to sub-second, on a laptop CPU, with no API key, no
rate limit, and no network dependency. That is what makes per-residue scanning affordable
during generation rather than as a batch job afterwards.

The student timing includes the ESM-2 forward pass, which dominates it — the head itself
is negligible — so this is the honest end-to-end cost rather than a flattering
measurement of the small model alone. Two caveats bound the claim: this compares local
compute against a network round trip rather than model against model, and a local
NetMHCpan binary would close most of the gap. In the other direction, the single timed
API call does not capture IEDB rate limiting, which is what actually made corpus-scale
labelling a multi-hour job.

## What we do not claim

- **No whole-response axis.** Every number above measures agreement with NetMHCpan.
  Nothing here is trained on TCR recognition, cytokines, antibody responses, danger,
  or animal outcomes.
- **No direct mechanistic validation yet.** BA still requires an independent measured
  HLA-A*02:01 affinity cohort. EL and the processing stack require independent human
  HLA-A*02:01 immunopeptidomics.
- **PDA is training data, not a benchmark.** Once the full PDA corpus was used for
  training, it stopped being available for out-of-distribution claims. The NY-ESO
  cohort is the only external evaluation we report.
- **One allele.** HLA-A\*02:01 and 9-mers only. Population-level risk needs a panel.
- **Teacher-bounded.** A distilled student cannot exceed its teacher. The value is
  speed and per-residue integration with cleavage and TAP, not accuracy beyond
  NetMHCpan.
- **Predicted affinity is an imitation, not a measurement.** The IC50 the profile
  reports is what NetMHCpan would predict, reproduced within about 2.3× on binders. It
  is not an experimental binding constant, and the fold error is too wide to compare two
  candidates that sit within roughly threefold of each other.
- **Accessibility is a design metric, not a biological one.** Burial does not gate MHC-I
  presentation. It predicts only whether a flagged epitope can be resampled safely.
- **Cleavage and TAP are inherited.** Those heads come from the legacy checkpoint and
  carry its MHCflurry-era training distribution. We replaced the MHC lane, not all of
  chao1.
- **The speed advantage is 13.6×, not orders of magnitude.** Measured, not asserted (§6),
  and it compares local compute against a network round trip rather than model against
  model. A local NetMHCpan binary would close most of it.

## Reproduce

```bash
# 1. corpus: label all PDA parents with NetMHCpan 4.1 (detached CPU job)
uv run --extra proto modal run --detach scripts/build_netmhcpan_corpus_modal.py \
  --run-name a0201-netmhcpan-pda-full-v1 --pda-only-full --pda-challenge-rows 0

# 2. embeddings + five-fold student, three channels (detached A10G job)
uv run --extra proto modal run --detach scripts/train_mhci_student_modal.py \
  --corpus-run-name a0201-netmhcpan-pda-full-v1 \
  --run-name a0201-netmhcpan-pda-cv5-v4 --epochs 60 --hidden-dim 512

# 3. occurrence-level profile with chao1 and out-of-fold student lanes
uv run --extra proto modal run --detach scripts/build_mhci_profile_modal.py \
  --source-corpus-run a0201-netmhcpan-pda-full-v1 \
  --student-run a0201-netmhcpan-pda-cv5-v4 \
  --run-name a0201-pda-mhci-profile-v4

# 4. evaluations (local)
uv run --extra ml python scripts/compare_mhc_heads_pda.py
uv run --extra ml python scripts/benchmark_external_nyeso.py
uv run --extra ml python scripts/calibrate_screening_thresholds.py
uv run --extra ml python scripts/build_epitope_accessibility.py
uv run --extra ml python scripts/benchmark_screening_latency.py
```

Artifacts and metrics:

```text
data/processed/profiles/a0201-pda-mhci-profile-v4/   profile parquet + manifest
models/a0201-netmhcpan-pda-cv5-v4/checkpoint/        five fold heads + metrics.json
results/benchmarks/pda_mhc_head_comparison/          old head vs new head
results/benchmarks/nyeso_a0201/                      external cohort
results/benchmarks/screening_calibration/            recall-driven cutoffs
results/benchmarks/epitope_accessibility/            burial vs epitope risk
results/benchmarks/screening_latency/                scanning cost vs the teacher
```

Every manifest records the checkpoint SHA-256, the corpus SHA-256, the teacher version
and endpoint, and the caveats that travel with the numbers.

## Citations

- Visani et al. (2025) T cell receptor specificity landscape revealed through de novo
  peptide design. *PNAS*. doi:10.1073/pnas.2504783122 — source of the external cohort.
- Reynisson et al. (2020) NetMHCpan-4.1 and NetMHCIIpan-4.0. *Nucleic Acids Research*.
- Butcher et al. (2025) De novo Design of All-atom Biomolecular Interactions with
  RFdiffusion3. bioRxiv. doi:10.1101/2025.09.18.676967
- Lin et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with
  a language model. *Science*. (ESM-2)
