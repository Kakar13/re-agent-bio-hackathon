# Code and model synthesis

## What the system is

The repository contains an inspectable design-to-screen workflow with two
scientifically distinct layers:

1. An agent and campaign layer that researches a target, builds typed design
   calls, requires approval before compute, validates generated candidates, and
   preserves provenance.
2. An immune-evidence layer that keeps MHC-I processing, MHC-II presentation,
   CD4-response evidence, tolerance, and structure-linked editability as
   separate channels.

The system's defensible product claim is campaign triage. It helps decide which
candidate or peptide window deserves follow-up. It does not produce a calibrated
probability of human immunogenicity.

## Runtime path

### Agent orchestration

- `src/re_agent/agent/graph.py` selects screening profiles and injects the
  available model checkpoints into candidate, preflight, and campaign tools.
- `src/re_agent/agent/tools.py` resolves only repository-owned model artifacts,
  instantiates the MHC-I adapters, screens structurally eligible candidates,
  and attaches residue tracks.
- `src/re_agent/design/campaign.py` materializes the generation and validation
  plan.
- `src/re_agent/agent/run_store.py` preserves run and artifact state.
- `docs/ORCHESTRATION.md` states the fail-closed policy: only candidates with
  structural status `pass` are screened, and missing CD4-response evidence
  keeps the MHC-II combined rank null.

The deployed graph points to
`models/a0201-netmhcpan-pda-cv5-v4/checkpoint` when that directory is present.
If it is absent, the runtime degrades to the legacy Chao1-only lane. This
fallback is operationally useful, but it means every demo artifact must record
which checkpoint path was actually used.

### Hybrid MHC-I adapter

`src/re_agent/immuno/e2e_pls_pickle.py` loads the restricted Chao1 artifact and,
when configured, the NetMHCpan student ensemble.

The hybrid lane is not one jointly trained three-head network:

- Chao1 supplies N-terminal cleavage.
- Chao1 supplies C-terminal cleavage.
- Chao1 supplies TAP transport and bootstrap uncertainty.
- The NetMHCpan student replaces the legacy MHCflurry-derived MHC output.
- The student emits EL rank propensity, BA rank propensity, and BA affinity.

This distinction is central. The latest deck still depicts the legacy unified
Pepsickle, TAP, and MHCflurry training story. The manuscript's headline result
is the replacement MHC head trained on PDA.

## Exact score semantics

The code exposes two MHC-related risk fields and one protein summary. They must
not be collapsed into one unnamed "immunogenicity score."

### `overall_mhci_risk`

Defined in `src/re_agent/e2e_pls/netmhcpan_student.py` by
`format_mhci_profile`.

- Value: student EL presentation propensity.
- Binder class: derived from predicted EL rank.
- Conventional thresholds: strong at rank at most 0.5 percent, weak at rank at
  most 2 percent.
- Screening flag: a separate EL-propensity threshold of 0.7892, selected for
  95 percent recall of teacher strong binders on the out-of-fold PDA profile.

This field is a teacher-like MHC-I presentation score. It is not the composite
processing score.

### `composite_processing_risk`

Defined in `src/re_agent/immuno/e2e_pls_pickle.py`.

For student-backed inference:

```text
composite_processing_risk =
    (N-cleavage probability
     * C-cleavage probability
     * BA binding propensity)^(1/3)
```

EL and TAP are displayed separately. BA is used in the geometric composite so
that explicit cleavage is not combined with the more presentation-oriented EL
channel.

For the legacy Chao1-only path, the third factor falls back to Chao1's
MHCflurry-derived presentation propensity. Therefore the exact composite
definition is checkpoint-dependent and must be recorded with the result.

The geometric mean is an inspectable processing heuristic. The three factors are
not independent biological probabilities, and the result is not a joint
probability.

### `pathway_rank_score`

Campaign ranking does not average EL, BA, and affinity. Affinity is a monotone
transform of BA, and EL already includes presentation information. Ranking uses

```text
pathway_rank_score =
    EL^0.70 * (sqrt(N-cleavage * C-cleavage))^0.30
```

EL is the presentation endpoint. Proteasomal generation is a soft gate. The
0.70/0.30 split is a prior, not a fitted coefficient.

### Protein-level `overall.score`

The adapter sorts all 9-mer pathway-rank scores and reports the mean of the five
largest values. This top-five mean reduces dependence on one maximum window, but
it remains an uncalibrated campaign-ranking heuristic.

## Model lineage

### Legacy Chao1

- Checkpoint: `models/chao1/cv5_heads.pkl 2`.
- Encoder: frozen ESM-2 650M.
- Associated corpus: 5,000,613 teacher-labeled 9-mers from 28,741 natural and
  de novo parents.
- Labels: Pepsickle cleavage, MHCflurry presentation, and a TAP ridge associated
  with 613 measured TAP peptides.
- Validation boundary: the repository model card states that the artifact does
  not contain a held-out metrics artifact, a dataset manifest, or a training
  snapshot sufficient to reproduce claimed performance.

The legacy checkpoint remains useful for inherited cleavage and TAP features.
Its MHC head is not the current headline model.

### NetMHCpan PDA student

- Deployment model version: `mhci-netmhcpan-student-v3`.
- Repository bundle: `models/a0201-netmhcpan-pda-cv5-v4/checkpoint`.
- Encoder: frozen `esm2_t33_650M_UR50D`, 1,280 dimensions.
- Input representation: mean of nine peptide-residue vectors, excluding
  special tokens.
- Head: layer normalization, a 256-unit GELU hidden layer, dropout, and three
  sigmoid outputs.
- Outputs: EL rank propensity, BA rank propensity, and BA affinity score.
- Deployment: arithmetic mean of five fold heads.

The 650M-parameter encoder still runs at inference. The task head is compact,
but the deployed system is not a fully compressed small model.

## Data lineage and leakage control

The current PDA student uses:

- 1,602 de novo parent proteins.
- 290,937 overlapping 9-mer occurrences before exact-peptide deduplication.
- 160,012 unique peptides in the completed corpus.
- 608 connected parent components.

Parents that share any exact 9-mer are connected before fold assignment. Whole
components stay in one fold. This is a strong control against direct peptide
reuse across cross-validation boundaries.

The downstream v4 profile contains 289,335 occurrences and 159,038 unique
peptides after profile-specific availability and filtering. These lower counts
are not a contradiction if they are labeled as profile counts rather than full
corpus counts.

The control does not remove every possible dependency:

- Near-identical but non-exact peptides can span folds.
- Shared structural families can span folds.
- All PDA folds contribute to the deployed five-head ensemble.
- ESM-2 pretraining overlap with individual peptide fragments was not audited.

PDA therefore supports grouped teacher-imitation evaluation, not an external
out-of-distribution claim.

## Quantitative evidence

### Held-out teacher imitation

Across five component-grouped held-out test folds:

- EL Spearman mean 0.9431, population SD 0.0043.
- BA Spearman mean 0.9478, population SD 0.0043.
- EL binder AUPRC mean 0.8914, population SD 0.0090.
- BA binder AUPRC mean 0.8595, population SD 0.0028.
- EL binder prevalence mean 0.0665.
- BA binder prevalence mean 0.0341.

These results establish strong NetMHCpan fidelity under the chosen split. They
do not establish measured presentation or immune-response accuracy.

### Legacy-head comparison

On the same 159,038 profile peptides against NetMHCpan EL:

- Chao1 MHCflurry head: Spearman 0.354 and strong-binder average precision
  0.060.
- Student out-of-fold EL: Spearman 0.943 and strong-binder average precision
  0.820.

This establishes that the replacement is a better NetMHCpan surrogate. It is a
cross-teacher comparison, not proof that MHCflurry is biologically inferior.

### Independent NY-ESO-1 peptide cohort

The cohort contains 95 independently designed HLA-A*02:01 peptides and 32
reported T-cell responders.

- Student EL versus fresh NetMHCpan EL: Spearman 0.909.
- Student BA versus fresh NetMHCpan BA: Spearman 0.898.
- Weak-binder call agreement: 97.9 percent.
- NetMHCpan EL versus activation: AUROC 0.728.
- Student EL versus activation: AUROC 0.672.

The teacher and student confidence intervals overlap. Also, 93 of 95 peptides
were already NetMHCpan weak binders or better, creating severe range
restriction. This is useful out-of-corpus fidelity evidence and a negative
biological limits result, not superiority evidence.

### Operating point

The conventional predicted strong-binder rule recovers 52.3 percent of teacher
strong binders. The separate 0.7892 EL-propensity screening threshold recovers
95.0 percent at 43.9 percent precision and flags 8,563 of 159,038 peptides.

This is a product threshold calibrated to the teacher, not a biological
sensitivity estimate.

### Structure-linked editability

For 186,485 windows with matched structure evidence:

- 80.2 percent of strong windows were buried.
- 77.1 percent of weak windows were buried.
- 67.7 percent of nonbinders were buried.
- 22.2 percent of strong or weak windows were surface-exposed.

Solvent accessibility indicates mutation difficulty. It does not gate antigen
processing, because peptides arise after protein degradation.

### Latency

For a 156-residue sequence with 147 windows:

- Local CPU hybrid inference: 0.761 seconds best of three.
- IEDB NetMHCpan API request: 10.336 seconds.
- Operational ratio: 13.57 times.

The local timing includes ESM-2, Chao1 cleavage and TAP, and the student
ensemble. The comparison includes network transport for the teacher and is not
a compute-equivalent local-binary benchmark.

## Mechanistic scope

The current student is restricted to human HLA-A*02:01 MHC-I.

Its outputs map to two different experimental targets:

- BA estimates the peptide-HLA binding event and should be validated against
  quantitative IC50 or \(K_d\).
- EL estimates ligand presentation and should be validated against human
  HLA-A*02:01 eluted-ligand or monoallelic immunopeptidomics data.

BA does not establish that a protein is cleaved into the peptide. The inherited
cleavage and TAP outputs model upstream steps, and the combined processing
claim requires naturally presented ligand evidence. T-cell recognition,
cytokines, antibody responses, danger, and animal outcomes are downstream of
the modeled mechanism and are excluded from primary model evaluation.

## Generation evidence

The preserved main-branch summary reports 500 RFdiffusion3 backbone complexes:
100 for each of five targets.

The defensible status is narrower than a completed design-to-screen campaign:

- The backbone PDB files are not in the review snapshot.
- Sequence design, refolding, interface confidence, binding, and immune
  screening were not completed on those exact 500 backbones.
- The PD-L1 QC supports 15 of 16 designs contacting all three hotspots within
  4.5 angstrom and 16 of 16 within 6 angstrom.
- The five reported per-target GPU times sum to 6,147.7 seconds, or 1.71 hours.

The full agent workflow is implemented in code, but a retained, explicitly
approved, low-budget live generation-to-screen trace remains the cleanest proof
that all deployed stages work together.
