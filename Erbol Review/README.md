# Erbol Review

External-review bundle for the CHAO / re:AGENT submission.

## Bottom line

The core machine-learning result is scientifically defensible as a fast,
HLA-A*02:01 NetMHCpan teacher surrogate for design-time liability triage. It is
not a validated predictor of human immunogenicity, anti-drug antibodies, or
clinical safety.

The current Sundial manuscript and project plan largely preserve that boundary.
The latest slide deck from `origin/main` does not yet match the manuscript and
should not be presented unchanged.

Three corrections are submission-critical:

1. The deck describes the legacy 5,000,613-row Pepsickle, MHCflurry, and TAP
   three-head system. The current manuscript's headline model is the separate
   PDA-trained NetMHCpan student: 1,602 parents, 160,012 unique 9-mers, and
   component-grouped five-fold teacher-imitation evaluation.
2. The deck says all 16 PD-L1 controls contacted all three hotspots within
   4.5 angstrom. The preserved QC artifact supports 15 of 16 at 4.5 angstrom
   and 16 of 16 at 6 angstrom.
3. The deck reports 5.6 GPU hours for 500 backbones. The preserved per-target
   summary totals 6,147.7 seconds, or 1.71 hours. Use 1.71 hours unless a
   separate auditable timing artifact supports 5.6 hours.

## Recommended reviewer framing

Use this one-sentence claim:

> A frozen ESM-2 adapter reproduces NetMHCpan 4.1 HLA-A*02:01 rankings with high
> held-out fidelity and enables local, sub-second liability triage inside a de
> novo design workflow, while remaining explicitly separate from biological
> immunogenicity prediction.

Use these boundaries every time:

- The model imitates NetMHCpan. It does not beat its teacher on the teacher's
  own target.
- PDA is training data, not an external out-of-distribution benchmark.
- The 95-peptide NY-ESO-1 cohort is independent of PDA, but it was heavily
  pre-screened for NetMHCpan binding and cannot establish model superiority.
- Animal outcomes, PBMC cytokines, T-cell activation, and antibody responses
  are outside the primary benchmark because the model does not represent those
  downstream mechanisms.
- BA must be benchmarked against measured human HLA-A*02:01 affinity. EL and
  the processing stack must be benchmarked against human HLA-A*02:01
  immunopeptidomics or eluted-ligand evidence.
- BA does not measure peptide generation. EL is the closer endpoint for
  naturally processed and presented ligands.
- The reported speed comparison is local CPU inference versus a network API,
  not local model versus local NetMHCpan compute.

## Bundle contents

- `slides/chao_final_deck/index.html`: exact latest CHAO deck from
  `origin/main`, deck commit `0f1a410`.
- `sundial/main.pdf`: compiled Sundial manuscript.
- `sundial/main.tex`, `sundial/refs.bib`: manuscript source and bibliography.
- `sundial/PROJECT_PLAN.md`: consolidated narrative, demo script, gaps, and
  self-rubric.
- `references/ARCHITECTURE.md`: what is frozen, what is trained, and exactly which
  chao1 head was replaced. Read this before the manuscript if the legacy/replacement
  boundary is unclear.
- `references/`: methods, orchestration, and model-card snapshots.
- `evidence/model/`: deployed student manifest and cross-validation metrics.
- `evidence/benchmarks/`: external, calibration, latency, accessibility, and
  qualitative case-study outputs.
- `evidence/rfd3_binders/`: main-branch generation summary and PD-L1 control QC.
- `SCIENTIFIC_CONSISTENCY_AUDIT.md`: prioritized validity and consistency
  review.
- `CODE_SYNTHESIS.md`: implementation map and exact score semantics.
- `DECK_REDLINE.md`: required slide corrections.
- `MOCK_REVIEWER_QA.md`: journal-style questions and defensible answers.
- `MANIFEST.json`: source revisions and SHA-256 checksums for the bundle.

## Ten-minute review order

1. Read this file and `SCIENTIFIC_CONSISTENCY_AUDIT.md`.
2. Skim `references/ARCHITECTURE.md` — it is the fastest way to see that the legacy
   and replacement models are distinct, which is correction 1 below.
3. Open `sundial/main.pdf`.
4. Use `CODE_SYNTHESIS.md` to connect each paper claim to code and evidence.
5. Review `DECK_REDLINE.md` before rehearsing the slides.
6. Rehearse `MOCK_REVIEWER_QA.md` with one teammate challenging every answer.

## Decision

The submission is credible if it is presented as an inspectable
design-to-screen system with a teacher-distilled MHC-I liability lane. It is
not credible if the deck merges the legacy and replacement models, reports
unsupported timing or control-gate values, or uses "immunogenicity score" as if
it were a calibrated human outcome probability.
