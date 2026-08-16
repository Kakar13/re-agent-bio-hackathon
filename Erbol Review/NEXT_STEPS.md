# Erbol Review next steps

## Purpose

This document is the handoff for updating the presentation and the live Sundial paper. It records what must change without applying additional edits to Sundial or the slide deck.

## Baseline synchronized from main

- `origin/main` was fetched and merged into `chao2`.
- Main revision reviewed: `d6feb4f9fb8af3990cb951e47d7efbcf342b723d`.
- The latest deck is `deck/chao_final_deck/index.html`.
- The review snapshot at `Erbol Review/slides/chao_final_deck/index.html` has the same Git blob hash as the latest deck: `ac3e92adbb6f5f5c1d01c201a4210c1c82679e15`.
- The current scientific reference is `Erbol Review/sundial/main.tex`, supported by `PROJECT_PLAN.md`, `CODE_SYNTHESIS.md`, and the preserved benchmark evidence.

Repeat the hash comparison after any later merge from main. Do not redline an older deck snapshot.

## Scientific statement that all artifacts must share

The current local model is a teacher-distilled human HLA-A*02:01 MHC-I screening module for campaign downselection.

The outputs have different meanings:

1. BA estimates peptide-HLA binding.
2. EL estimates ligand presentation and contains information learned from eluted-ligand evidence.
3. Cleavage and TAP are separate upstream predictions inherited from Chao1.
4. The geometric processing composite combines N-cleavage, C-cleavage, and BA.
5. TCR recognition, danger, cytokines, antibody responses, animal outcomes, and whole immune response are not model targets.

Current quantitative results establish NetMHCpan teacher fidelity. They do not establish measured human HLA affinity or natural ligand presentation.

## Sundial paper changes

Use `Erbol Review/sundial/main.tex` as the proposed source and apply the following changes to the live Sundial paper.

### Abstract and introduction

- Lead with de novo campaign downselection and the cost of direct HLA experiments.
- Describe the model as HLA-A*02:01 binding and presentation screening.
- State that BA, EL, cleavage, and TAP are separate outputs.
- Remove language suggesting prediction of a whole human immune response.
- State that the student is a local NetMHCpan surrogate, not a model that beats NetMHCpan.

### Methods

- Keep the full-PDA teacher-distillation lineage and exact-9-mer connected-component split.
- Correct the student output description to EL, BA, and BA-derived affinity.
- Define the processing composite as the geometric mean of N-cleavage, C-cleavage, and BA propensity.
- State explicitly that BA does not validate peptide generation.
- Define direct validation targets:
  - BA against quantitative human HLA-A*02:01 IC50 or \(K_d\).
  - EL and the processing stack against human HLA-A*02:01 monoallelic immunopeptidomics or confidently assigned eluted ligands.

### Results

- Retain the five-fold teacher-fidelity values from the deployed v4 checkpoint.
- Retain the same-row Chao1 versus student comparison only as evidence for replacing the legacy MHC head.
- Treat NY-ESO-1 only as out-of-corpus teacher fidelity.
- Remove Neo-2/15, animal outcomes, Latent-X2 PBMC outcomes, and downstream T-cell activation as model-validation cases.
- Keep 9S14 only as a traceable sequence-to-structure visualization, not an external validation case.

### Discussion and limitations

- Present the product as a fast, inspectable filtering layer before direct HLA experiments.
- Name the missing independent affinity and immunopeptidomics benchmarks.
- Keep the one-allele limitation and the teacher-bounded limitation visible.
- Do not imply population coverage, measured immunogenicity, T-cell recognition, ADA prediction, or clinical safety.

### Sundial execution and verification

1. Obtain a fresh Sundial connection prompt from the workspace owner.
2. Do not store the bearer token or a credential-bearing prompt in Git.
3. Connect, read the current live `main.tex`, `refs.bib`, and `PROJECT_PLAN.md`, and compare them with this review snapshot.
4. Apply all edits to each file in one batched edit request.
5. Compile with Tectonic in Sundial.
6. Resolve undefined references, missing citations, overflow, and figure-placement errors.
7. Visually inspect the final PDF.
8. Export the new `main.tex`, `refs.bib`, `PROJECT_PLAN.md`, and `main.pdf` into `Erbol Review/sundial/`.
9. Commit the refreshed snapshot separately from unrelated model or workbench changes.

## Slide changes

Use `Erbol Review/DECK_REDLINE.md` for detailed slide-by-slide language. The required high-level changes are:

### Model lineage

- Replace the legacy-only 5-million-row Pepsickle/MHCflurry story with the current hybrid:
  - Chao1 supplies cleavage and TAP.
  - The deployed v4 student supplies NetMHCpan-like EL, BA, and affinity outputs.
- Keep the legacy corpus only as lineage, not as the evidence for the replacement MHC head.

### Claims and metrics

- Label all EL and BA metrics as teacher fidelity.
- Do not say the student beats NetMHCpan.
- Do not call PDA an external OOD benchmark.
- Use the deployed v4 values consistently:
  - EL Spearman: 0.943.
  - BA Spearman: 0.948.
  - External NY-ESO-1 teacher fidelity: EL 0.909 and BA 0.898.
- Correct the PD-L1 control to 15 of 16 at 4.5 angstrom and 16 of 16 at 6 angstrom.
- Use 1.71 hours for the preserved per-target GPU summaries, not 5.6 hours.
- Describe 500 designs as reported backbone summaries, not a completed end-to-end sequence-to-screen campaign.

### Required new slide

Add one mechanistic-boundary slide with four stages:

1. Proteolytic generation and TAP transport.
2. Peptide-HLA binding, represented by BA.
3. Natural ligand presentation, represented by EL.
4. TCR recognition and downstream response, not modeled.

The slide should end with the two direct next experiments: measured HLA-A*02:01 affinity and monoallelic HLA-A*02:01 immunopeptidomics.

### Closing language

Use:

> A fast, local, inspectable HLA-A*02:01 screen that helps scientists downselect de novo campaigns before direct HLA experiments.

Avoid:

- Predicts human immunogenicity.
- Predicts danger or immune recognition.
- Validated in animals or humans.
- Beats NetMHCpan.
- Completes the in silico to in vivo bridge.

## Acceptance checklist

- [ ] Latest main has been merged before editing the deck.
- [ ] Review deck hash matches the latest main deck before redlining.
- [ ] Sundial paper and project plan use the same BA, EL, cleavage, TAP, and composite definitions.
- [ ] Animal, PBMC, T-cell activation, and whole-response cases are absent from primary validation claims.
- [ ] All displayed metrics come from the deployed v4 manifest.
- [ ] Slide numerical corrections match preserved evidence.
- [ ] The mechanistic-boundary slide is present.
- [ ] The live Sundial PDF compiles and is visually inspected.
- [ ] The refreshed PDF and source are copied back into `Erbol Review/sundial/`.
- [ ] No Sundial token or other credential is committed.
