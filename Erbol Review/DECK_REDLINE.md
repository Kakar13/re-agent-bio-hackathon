# Slide deck redline

Source preserved at `slides/chao_final_deck/index.html`.

Do not overwrite the preserved source. Apply these changes to the presentation
copy, then export a new revision with a new checksum.

## Global language

Replace:

- "immunogenicity score" with "immune-liability profile" or the exact score
  field.
- "predicts immunogenicity" with "predicts HLA-A*02:01 MHC-I teacher-like
  processing and presentation liabilities."
- "small model" with "small task head over a frozen 650M-parameter encoder."
- "validation" with "teacher-imitation evaluation" unless the endpoint is
  experimentally measured.

Add a persistent footer:

> HLA-A*02:01 MHC-I triage only. Teacher-distilled and not validated for human
> immunogenicity or ADA.

## Slide 01, CHAO

### Current problem

The slide says "Three heads, one vector" and "5M+ 9-mers." That describes the
legacy Chao1 training story, not the current NetMHCpan student.

### Replace with

**What it is**

> A hybrid, inspectable MHC-I profile. Chao1 supplies cleavage and TAP; a
> frozen-ESM student supplies NetMHCpan-like EL, BA, and affinity outputs.

**What the current MHC head is built on**

> 160,012 unique PDA 9-mers from 1,602 de novo parent proteins, grouped into 608
> connected components before five-fold assignment.

**What it does**

> Flags HLA-A*02:01 liabilities for review and resampling. It does not predict
> clinical immunogenicity.

Remove "0 of 29,353 clusters leaking across splits" unless the deck explicitly
labels it as a legacy-dataset result. For the current student, say:

> Exact shared 9-mers never cross component-grouped folds.

## Slide 02, blind spot

The overall motivation is strong, but the slide should not imply that the
current model closes every immunological blind spot.

Add:

> For extracellular therapeutics, MHC-II/CD4 and B-cell biology remain primary
> gaps. This completed local model is one MHC-I evidence lane.

## Slide 03, three projections

### Current problem

The MHC diagram depicts a 64-dimensional projection and cosine to a binder
centroid. That is the legacy MHCflurry-derived head.

### Replace with

- Mean-pool the nine ESM-2 residue vectors.
- Layer normalization.
- 256-unit GELU hidden layer with dropout.
- Three outputs: EL rank propensity, BA rank propensity, BA affinity score.
- Arithmetic mean across five deployment heads.

Keep Chao1 cleavage and TAP as visually separate inherited components. Do not
draw all outputs as one jointly trained model.

## Slide 04, what taught it

### Current problem

The entire slide describes the 5M-row legacy Pepsickle, MHCflurry, and TAP
corpus. It cannot serve as evidence for the current student metrics.

### Recommended replacement

Use the current PDA lineage:

1. 1,602 PDA parents.
2. 290,937 9-mer occurrences before deduplication.
3. 160,012 unique peptides.
4. NetMHCpan 4.1 EL, BA rank, and BA affinity labels.
5. Exact-9-mer parent graph with 608 components.
6. Five component-grouped folds.
7. Out-of-fold teacher-imitation metrics.

If legacy lineage is important, add one small note:

> Chao1 cleavage and TAP are inherited from an earlier 5M-row
> teacher-labeled system and have not been revalidated on the PDA corpus.

## Slide 05, reading the score

Show three named outputs:

1. `overall_mhci_risk`: EL presentation propensity.
2. `composite_processing_risk`: geometric mean of N-cleavage, C-cleavage, and
   BA binding propensity.
3. Protein summary: mean of the five highest pathway scores `S_i`.

Add:

> These are unitless triage heuristics and teacher-like propensities, not immune
> response probabilities.

Do not display the geometric composite without its three component values.

## Slide 06, steer to safety

Keep the steering workflow only if the demo actually uses the legacy
`src/re_agent/e2e_pls/steer.py` path.

Add:

- Function and binding are not preserved by ESM likelihood alone.
- Mutations require structure-aware and experimental follow-up.
- Lower model score is not equivalent to deimmunization.

If the steering output was not produced by the v4 hybrid checkpoint, label it
as a legacy-model demonstration.

## Slide 07, what it unlocks

### Correct the title

Replace "five targets nobody can drug" with:

> 500 generated backbone complexes across five challenging targets and one
> positive-control context.

PD-L1 is druggable and is explicitly the positive control.

### Correct the timing

Replace "5.6 h" with:

> 1.71 reported GPU hours across the five per-target summaries.

Replace "GPU time to train all 500 backbones" with:

> Reported GPU time to generate 500 backbone complexes.

### Correct the control result

Replace:

> All 16 contacted all three hotspots within 4.5 angstrom.

With:

> 15 of 16 contacted all three configured hotspots within 4.5 angstrom; all 16
> did so within 6 angstrom.

### Add the evidence boundary

> Backbone counts are preserved in summary artifacts, but the 500 PDB files are
> absent from this review snapshot. Those exact candidates were not completed
> through sequence design, refolding, binding, and immune screening.

## Slide 08, system architecture

Distinguish implementation from demonstrated execution.

Use these states:

- **Implemented and tested:** graph state, typed tool adapters, approval
  interrupt, structural gates, MHC-I adapter, MHC-II EL/BA lane, deterministic
  review, workbench rendering.
- **Demonstrated separately:** 500-backbone RFdiffusion3 summary, PD-L1
  geometric control, local MHC-I benchmarks, one-sequence preflight.
- **Still required:** one retained approved low-budget generation-to-screen
  trace over the same candidate.
- **Unavailable:** calibrated CD4-response artifact; combined MHC-II/CD4 rank
  remains null.
- **Not implemented:** Benchling handoff.

Do not use one "built" badge for both code existence and end-to-end execution
evidence.

## Slide 09, student and teacher

### Current problem

The slide again names Pepsickle, MHCflurry, and TAP ridge as the student
teachers and says the student is trained on 5M windows. That is not the current
NetMHCpan student.

### Replace teacher branch

> PDA 9-mers to NetMHCpan 4.1 EL rank, BA rank, and BA affinity.

### Replace student branch

> Frozen ESM-2 mean-pooled 9-mer to a 256-unit head with three outputs,
> ensembled across five folds.

### Add result strip

- EL Spearman: 0.9431 plus or minus 0.0043.
- BA Spearman: 0.9478 plus or minus 0.0043.
- EL AUPRC: 0.8914 against prevalence 0.0665.
- BA AUPRC: 0.8595 against prevalence 0.0341.

Caption:

> Held-out component-grouped teacher imitation. Direct human HLA-A*02:01
> affinity and ligand-presentation validation remains required.

## Add one external-evidence slide

Include:

- 95 independent NY-ESO-1 designed peptides.
- Student versus teacher EL Spearman 0.909.
- Student versus teacher BA Spearman 0.898.
- Weak-binder call agreement 97.9 percent.

Headline:

> The student transfers teacher ranking out of corpus. This is still teacher
> fidelity, not measured HLA validation.

## Add one mechanistic-validation slide

Show four separate stages:

1. Proteolytic generation and TAP transport.
2. Peptide-HLA binding affinity, represented by BA.
3. Natural ligand presentation, represented by EL.
4. TCR recognition and downstream immune response, not modeled.

State prominently:

> BA must be benchmarked against measured human HLA-A*02:01 IC50 or \(K_d\).
> EL and the processing stack must be benchmarked against human HLA-A*02:01
> immunopeptidomics. Animal and whole-response outcomes are not model labels.

## Closing slide

Use:

> CHAO makes one early evidence lane fast, local, and inspectable so scientists
> can spend wet-lab time on the most defensible designs.

Avoid:

- Predicts immunogenicity.
- Deimmunizes proteins.
- Beats NetMHCpan.
- Validated in humans.
- Completes the full in silico to in vivo bridge.
