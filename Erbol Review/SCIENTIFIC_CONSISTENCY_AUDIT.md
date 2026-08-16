# Scientific consistency audit

## Review verdict

**Recommendation: major revision before external presentation.**

The central machine-learning technique is valid for the stated narrow task:
knowledge distillation from NetMHCpan into a small head over frozen ESM-2
representations. The grouped split, out-of-fold scoring, class-imbalance-aware
metrics, external teacher-fidelity test, and explicit claim boundaries are
strong.

The work does not validate immunogenicity prediction. It supports a fast,
inspectable MHC-I liability screen for campaign downselection.

The manuscript is substantially more defensible than the latest slide deck.
The deck currently mixes the legacy Chao1 model, the replacement NetMHCpan
student, and generation evidence into one story. Several deck claims contradict
the preserved evidence and must be corrected.

## Submission-blocking findings

### 1. The deck and manuscript describe different MHC models

**Evidence**

- `slides/chao_final_deck/index.html:1003-1019` describes three heads over one
  ESM-2 representation trained on 5M+ rows.
- `slides/chao_final_deck/index.html:1244-1285` names Pepsickle, MHCflurry, and
  the TAP ridge as the teachers.
- `sundial/main.tex:91-100` describes the current hybrid: Chao1 cleavage and TAP
  plus a PDA-trained NetMHCpan EL, BA, and affinity student.
- `src/re_agent/agent/graph.py:60-75` explicitly states that the NetMHCpan
  student replaces Chao1's MHCflurry-derived MHC lane.

**Scientific consequence**

A reviewer cannot tell which model produced the headline results. The
5,000,613-row legacy dataset has no reproduced held-out metric artifact in the
Chao1 model card, while the 160,012-peptide PDA student does.

**Required correction**

Choose the PDA NetMHCpan student as the headline model. Describe Chao1 only as
the inherited cleavage and TAP source. If the legacy data story remains, label
it as model lineage rather than current student training.

### 2. The PD-L1 4.5-angstrom control claim is false as written

**Evidence**

- `slides/chao_final_deck/index.html:1478` and the preserved report claim all 16
  controls contacted all three hotspots within 4.5 angstrom.
- `evidence/rfd3_binders/pdl1_control_qc.json:163-170` shows
  `pdl1_bb_0015.pdb` at 5.87 angstrom from hotspot 56 and `n_le_4.5: 2`.
- The complete QC supports 15 of 16 at 4.5 angstrom and 16 of 16 at 6 angstrom.
- `sundial/main.tex:161` already uses the correct values.

**Required correction**

Use: "15 of 16 PD-L1 controls contacted all three configured hotspots within
4.5 angstrom; all 16 did so within 6 angstrom."

### 3. The deck's 5.6-hour GPU claim lacks supporting evidence

**Evidence**

- `slides/chao_final_deck/index.html:1468-1478` reports 5.6 GPU hours.
- `evidence/rfd3_binders/summary.json` reports per-target elapsed times of
  847.1, 1,024.1, 738.2, 1,157.6, and 2,380.7 seconds.
- Those values sum to 6,147.7 seconds, or 1.71 hours.

**Required correction**

Use 1.71 reported GPU hours. Use 5.6 only if a separate immutable run artifact
defines what the extra time includes and reconciles the per-target summary.
Also replace "train all 500 backbones" with "generate 500 backbones."

### 4. The manuscript says the current student has two outputs, but v3 has three

**Evidence**

- `sundial/main.tex:98` describes a two-output head.
- `src/re_agent/e2e_pls/netmhcpan_student.py:20-29` defines three v3 outputs.
- `evidence/model/deployment_manifest.json:7-12` lists EL rank propensity, BA
  rank propensity, and BA affinity score.
- `sundial/PROJECT_PLAN.md:95-101` correctly lists all three.

**Required correction**

Change the methods to a three-output head. State that the first two outputs
target transformed percentile ranks and the third targets the bounded BA
affinity score that can be inverted to predicted IC50.

### 5. "Immunogenicity score" is not an acceptable name for the output

**Evidence**

- The student predicts teacher-like EL, BA, and affinity values.
- The hybrid adapter computes a geometric processing heuristic from Chao1
  cleavage and student BA.
- The system lacks a calibrated MHC-II/CD4 response model and withholds the
  corresponding combined rank.

**Scientific consequence**

Calling the number an immunogenicity score implies calibration to a measured
human endpoint that does not exist.

**Required correction**

Use one of these exact terms:

- "MHC-I presentation propensity" for student EL.
- "MHC-I binding propensity" for student BA.
- "composite processing risk" for the N-cleavage, C-cleavage, and BA geometric
  mean.
- "protein-level top-five processing score" for the sequence summary.
- "immune-liability triage profile" for the complete output vector.

### 6. The current generation evidence is not an end-to-end campaign

**Evidence**

- `evidence/rfd3_binders/report.md:92-95` defers ProteinMPNN, AF2/ESMFold
  validation, and disordered targets.
- The preserved branch contains summaries, but not the 500 generated PDBs.
- `sundial/main.tex:161` correctly says those exact backbones were not taken
  through sequence design, refolding, binding, or immune screening.
- `sundial/PROJECT_PLAN.md:402-409` still lists one approved live LangGraph
  campaign as submission-critical.

**Required correction**

Describe two separate achievements:

1. A reported 500-backbone RFdiffusion3 generation run with limited geometric
   control QC.
2. An implemented inspectable LangGraph design-to-screen workflow that still
   needs one retained low-budget live trace.

Do not imply that all 500 backbones passed the complete agent pipeline.

### 7. The paper's architecture figure is stale and visibly overlaps

**Evidence**

- `sundial/main.tex:66-68` renders Chao1 as implemented and the NetMHCpan
  student as an offline dotted node.
- The graph uses the v4 student automatically when the checkpoint is present.
- `evidence/rendered-architecture-overlap.png` shows overlapping Chao1 and
  student labels and an obscured edge annotation.

**Required correction**

Render one "Hybrid MHC-I profile" node with two incoming sublabels: "Chao1
cleavage/TAP" and "NetMHCpan student EL/BA/affinity." Mark the student as a
frozen local artifact, not incomplete. Increase horizontal spacing and move the
"displayed, not fused" annotation below the two MHC lanes.

## Major scientific comments

### The benchmark must stop at the modeled mechanism

The current model does not target a whole immune response. Its human
HLA-A*02:01 outputs have two direct mechanistic interpretations:

- BA estimates peptide-HLA binding and requires measured IC50 or \(K_d\).
- EL estimates ligand presentation and requires human eluted-ligand or
  monoallelic immunopeptidomics evidence.

BA does not validate peptide generation. Cleavage and TAP are separate upstream
model outputs. T-cell activation, cytokines, antibody responses, danger, and
animal outcomes are downstream and must not be used as primary accuracy labels.

### The grouped split is good, but it is not a complete OOD test

Connecting parents that share an exact 9-mer is a meaningful leakage control.
It prevents 40,125 reused peptides from crossing folds. It does not prevent
near-neighbor sequence, structural-family, or ESM-pretraining overlap.

Use "component-grouped held-out teacher imitation" for PDA. Use "independent
out-of-corpus cohort" for NY-ESO-1. Do not describe PDA as OOD after it has been
used to train all five deployment heads.

### The external cohort establishes teacher fidelity only

The student preserves teacher rankings on NY-ESO-1: EL Spearman is 0.909, BA
Spearman is 0.898, and weak-binder agreement is 97.9 percent. This is useful
out-of-corpus evidence.

The cohort's T-cell activation labels are excluded from the primary benchmark.
They add downstream TCR-recognition biology that the model does not represent.
The valid claim is external teacher fidelity, not experimental mechanistic
accuracy or whole-response prediction.

### The geometric score needs an explicit non-probabilistic interpretation

The N-cleavage, C-cleavage, and BA outputs are not demonstrated to be
independent or jointly calibrated. Their geometric mean is defensible as a
monotonic triage heuristic because it requires all three factors to be high.
It is not a measured pathway probability.

The repository also exposes `overall_mhci_risk`, which is EL propensity alone.
Every user-facing artifact should name both fields explicitly to avoid silent
semantic drift.

## Strengths worth emphasizing

- Exact-9-mer connected-component split control.
- Out-of-fold predictions for threshold calibration and legacy comparison.
- AUPRC reporting against explicit no-skill prevalence.
- Separate EL, BA, affinity, cleavage, TAP, confidence, and accessibility
  fields.
- A high-recall screening threshold kept separate from conventional binder
  terminology.
- Immutable manifests and checkpoint hashes.
- Honest reporting that the encoder remains 650M parameters.
- Fail-closed structural and MHC-II/CD4 orchestration.
- External out-of-corpus teacher-fidelity evidence.
- Explicit separation of BA, EL, cleavage, TAP, and downstream recognition.

## Required revision order

### Before the judging demo

1. Apply every correction in `DECK_REDLINE.md`.
2. Fix the manuscript's two-output description and architecture figure.
3. Use one canonical score vocabulary.
4. Recompile the PDF and visually inspect every figure and table.
5. Record which checkpoint hashes produced the displayed example.
6. Capture one approved low-budget live agent trace, or explicitly label the
   end-to-end execution as not yet demonstrated.

### Before a journal submission

1. Add an independent quantitative HLA-A*02:01 IC50 or \(K_d\) cohort.
2. Add monoallelic HLA-A*02:01 immunopeptidomics on expressed de novo
   constructs.
3. Audit overlap between the measured benchmark and NetMHCpan teacher training.
4. Compare frozen ESM-2 with simpler peptide baselines and an encoder-ablation
   baseline.
5. Add near-neighbor and structural-family split sensitivity analyses.
6. Add distance-to-training-manifold diagnostics and an abstention policy.
7. Revalidate or replace the inherited Chao1 cleavage and TAP heads.
8. Benchmark against a licensed local NetMHCpan binary for compute-equivalent
   latency.

## Final defensible conclusion

The project demonstrates a valid and useful engineering-science contribution:
an inspectable agent can place a fast, teacher-distilled human HLA-A*02:01
binding and presentation screen inside a de novo protein campaign and preserve
the evidence needed for review. The results justify downselection and
experimental prioritization. They do not establish measured HLA affinity or
natural presentation until those direct mechanistic benchmarks are completed.
