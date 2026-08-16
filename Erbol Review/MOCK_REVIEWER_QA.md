# Mock journal review and defense questions

## Reviewer summary

This work presents an inspectable agent that integrates literature-grounded
protein-design planning, approval-gated scientific tools, structural triage,
and a local MHC-I processing profile. Its main quantitative contribution is a
small multi-output head over frozen ESM-2 representations that imitates
NetMHCpan 4.1 on HLA-A*02:01 peptides from the Protein Design Archive.

The engineering integration is useful and the teacher-imitation results are
strong. The work is not yet a biological immunogenicity model. The current
evidence supports design-time liability triage and experiment prioritization.

## Mock recommendation

**Major revision.**

The central result is publishable in principle as a reproducible scientific
workflow and surrogate-model study, but the claims, model lineage, generation
evidence, and biological validation require clearer separation.

## Major reviewer comments with proposed responses

### Comment 1: The manuscript must define the exact modeled endpoint

**Reviewer concern**

The language sometimes moves from HLA binding and presentation to whole immune
response, although those are different endpoints.

**Defensible response**

We agree. The target is human HLA-A*02:01 peptide binding and
presentation-related processing. BA estimates the peptide-HLA binding event.
EL estimates ligand presentation. Cleavage and TAP remain separate upstream
predictions. TCR recognition, cytokines, antibodies, danger, and whole immune
response are not model targets.

**Revision**

Use the exact BA, EL, cleavage, TAP, and composite field names throughout.
State which measured HLA endpoint would validate each output.

### Comment 2: Why train a surrogate rather than call NetMHCpan directly?

**Reviewer concern**

The student cannot exceed a teacher it is trained to copy. The scientific value
of distillation is unclear.

**Defensible response**

The value is operational and integrative, not teacher-beating accuracy. The
local student avoids per-candidate network calls, returns all overlapping
windows in about 0.76 seconds on CPU for a 156-residue protein, preserves
versioned local artifacts, and supports inner-loop scanning. We explicitly do
not claim that the student is more accurate than NetMHCpan.

**Revision**

Report the comparison as local hybrid inference versus the IEDB network API.
Do not call it a compute-equivalent NetMHCpan speedup.

### Comment 3: The deck and paper appear to use different training datasets and teachers

**Reviewer concern**

The deck describes 5M Pepsickle, MHCflurry, and TAP-labeled windows, while the
paper reports a 160,012-peptide NetMHCpan PDA student.

**Defensible response**

These are two model generations. Chao1 is the legacy checkpoint and remains the
source of cleavage and TAP features. The current MHC head is a separate
NetMHCpan student trained on PDA and replaces Chao1's MHCflurry-derived head in
the hybrid runtime.

**Revision**

Rebuild the deck around the current hybrid. Retain the 5M corpus only as legacy
lineage and clearly state that its held-out performance is not reproduced by
the supplied Chao1 artifact.

### Comment 4: How was peptide leakage controlled?

**Reviewer concern**

Overlapping windows from related proteins can make peptide-level random splits
overly optimistic.

**Defensible response**

We tile every parent first, then connect any parents that share an exact 9-mer.
Whole connected components remain in one fold. This prevents more than 40,000
shared peptides from crossing the train-test boundary.

**Limitation**

The split does not guarantee separation of near-neighbor peptides, structural
families, or pretraining overlap. Those are planned sensitivity analyses.

### Comment 5: Is PDA actually out of distribution?

**Reviewer concern**

The paper motivates de novo distribution shift but trains on all PDA proteins.

**Defensible response**

PDA is not an external OOD test for the final model. It is the training
distribution. The shift exists relative to ESM-2's mostly natural-sequence
pretraining and relative to many conventional MHC datasets. The independent
NY-ESO-1 cohort provides out-of-corpus evidence, not a complete protein-level
OOD benchmark.

### Comment 6: Why is frozen ESM-2 a valid representation for isolated 9-mers?

**Reviewer concern**

ESM-2 was pretrained mostly on natural full-length proteins, not isolated de
novo 9-mers.

**Defensible response**

Frozen protein-language representations are a valid transfer-learning
strategy, and the high held-out and external teacher correlations show that
the representation contains sufficient information for this teacher task.
That does not prove calibrated extrapolation to every designed sequence. We
therefore retain a domain-shift limitation and plan distance-to-training and
abstention diagnostics.

### Comment 7: Why use only HLA-A*02:01?

**Reviewer concern**

One allele cannot support population-level conclusions.

**Defensible response**

We chose one common allele to produce one coherent, reproducible end-to-end
model during the project. The output is allele-specific. It is not population
coverage. Production use requires a predeclared multi-allele panel and
population-weighted evaluation.

### Comment 8: Does BA measure the probability that a protein is cleaved and presented?

**Reviewer concern**

The manuscript appears to conflate peptide-HLA binding with peptide generation
and downstream recognition.

**Defensible response**

No. BA models only peptide-HLA binding. Proteolytic cleavage and TAP transport
are separate upstream steps. EL is the closer presentation endpoint because its
training evidence includes naturally eluted ligands. TCR recognition is a
separate downstream step.

**Revision**

Present BA, EL, cleavage, TAP, and recognition as separate stages. Validate BA
with quantitative affinity and EL with human HLA-A*02:01 immunopeptidomics.

### Comment 9: Is the geometric composite biologically justified?

**Reviewer concern**

The score multiplies three model outputs without evidence that they are
independent or jointly calibrated.

**Defensible response**

It is a monotonic prioritization heuristic, not a mechanistic joint
probability. The geometric mean forces all three factors to contribute and
prevents a single large factor from dominating. We display each component
alongside the composite and do not interpret the result as calibrated risk.

### Comment 10: Why does the composite use BA rather than EL?

**Reviewer concern**

NetMHCpan EL is usually the preferred presentation channel.

**Defensible response**

EL is reported directly as the presentation-oriented output.
The separate composite already includes explicit cleavage features, so adding
EL could count related upstream presentation evidence twice. BA is the
binding-focused factor and is therefore the cleaner third component. This is a
design choice, not a demonstrated optimal biological fusion rule.

### Comment 11: Have the inherited cleavage and TAP heads been validated with the new student?

**Reviewer concern**

Replacing only the MHC head may create an incoherent composite.

**Defensible response**

No. The inherited Chao1 cleavage and TAP heads were not revalidated on the PDA
corpus or on the external cohort. We report them as separate legacy-derived
evidence and list their revalidation as a priority. The external NY-ESO-1 assay
uses synthetic pulsed peptides, so cleavage and TAP are out of scope there.

### Comment 12: Is the high-recall screening threshold independently calibrated?

**Reviewer concern**

The threshold is selected on the same set of out-of-fold predictions used to
report recall.

**Defensible response**

The predictions are out of fold, so no peptide was scored by a head trained on
its label. However, threshold selection and performance reporting use the same
pooled OOF set. The 95 percent recall is therefore an internal product
operating point against the teacher, not an independently validated sensitivity
estimate.

### Comment 13: Why report AUPRC?

**Reviewer concern**

Binder classes are highly imbalanced.

**Defensible response**

AUPRC is appropriate because EL and BA binder prevalence are approximately
6.7 and 3.4 percent. We report prevalence as the no-skill AUPRC baseline. EL
AUPRC 0.891 and BA AUPRC 0.860 are therefore meaningfully above baseline, but
only for teacher-defined binder labels.

### Comment 14: What does the NY-ESO-1 cohort actually validate?

**Reviewer concern**

The cohort includes a downstream activation endpoint that the model does not
target.

**Defensible response**

It validates out-of-corpus teacher fidelity only: student-versus-teacher
Spearman is 0.909 for EL and 0.898 for BA, with 97.9 percent weak-binder
agreement. We exclude T-cell activation from the primary model benchmark
because it depends on TCR recognition and other downstream mechanisms.

### Comment 15: What direct human benchmark is still missing?

**Reviewer concern**

Teacher imitation is not direct evidence that the model matches measured human
HLA binding or natural ligand presentation.

**Defensible response**

We need two matched benchmarks. First, an independent canonical 9-mer cohort
with quantitative human HLA-A*02:01 IC50 or \(K_d\) for BA. Second,
monoallelic or confidently assigned human HLA-A*02:01 immunopeptidomics with
source-protein context for EL and the processing stack. Animal and whole-PBMC
outcomes are not substitutes for either benchmark.

### Comment 16: Are the 500 generated designs auditable?

**Reviewer concern**

The deck presents a completed large-scale campaign.

**Defensible response**

The preserved branch contains per-target summaries for 500 RFdiffusion3
backbone complexes and limited PD-L1 geometric QC. It does not contain the 500
PDBs. Those exact candidates were not completed through sequence design,
refolding, binding, and immune screening. We now call this reported generation
scale, not an end-to-end campaign result.

### Comment 17: What exactly did the PD-L1 positive control establish?

**Reviewer concern**

Contact to configured hotspots does not establish binding.

**Defensible response**

Correct. Fifteen of 16 sampled backbones contacted all three configured
hotspots within 4.5 angstrom, and all 16 did so within 6 angstrom. This checks
that generation was geometrically directed toward the intended interface. It
does not establish sequence viability, folding, affinity, specificity, or
function.

### Comment 18: What is agentic rather than a fixed pipeline?

**Reviewer concern**

The workflow may be conventional orchestration.

**Defensible response**

The graph converts a natural-language objective into cited target research,
discovers exact external schemas, writes a versioned design specification,
materializes a tool plan, pauses for human approval, executes available tools,
validates outputs, routes eligible candidates into pathway-specific screens,
and runs deterministic review. The scientific gates are explicit and
inspectable rather than delegated to free-form model judgment.

### Comment 19: Has the complete agent path been demonstrated live?

**Reviewer concern**

Implemented stages and end-to-end evidence are conflated.

**Defensible response**

The components and bounded preflight paths are implemented and tested.
The separate 500-backbone run is preserved only as summaries. A retained,
approved, low-budget generation-to-screen trace over one candidate remains a
submission-critical demonstration. We will label that gap if it is not
completed.

### Comment 20: Is the latency comparison fair?

**Reviewer concern**

The comparison includes network overhead only for NetMHCpan.

**Defensible response**

It is an operational comparison: 0.761 seconds for local CPU hybrid inference
versus 10.336 seconds for one batched IEDB API request. It is not a local
model-compute comparison. A licensed local NetMHCpan binary is required for
that benchmark.

### Comment 21: What does structure accessibility contribute?

**Reviewer concern**

Buried peptides can still be processed and presented.

**Defensible response**

Accessibility is an editability field, not an immunogenicity field. It tells a
designer whether a flagged window is likely to tolerate surface resampling or
requires a backbone-aware stability check. It never gates antigen presentation.

### Comment 22: What uncertainty is calibrated?

**Reviewer concern**

The interface may imply calibrated confidence.

**Defensible response**

The TAP lane carries bootstrap uncertainty from the legacy ridge. The composite
confidence is decomposed but is not a calibrated probability of correctness or
coverage. The student does not yet provide calibrated epistemic uncertainty or
an OOD abstention rule.

### Comment 23: What baselines are missing?

**Reviewer concern**

High student fidelity may not require ESM-2.

**Defensible response**

We need amino-acid composition, one-hot peptide, motif/PWM, smaller language
model, and trainable-versus-frozen encoder baselines. These will quantify the
incremental value of ESM-2 rather than assuming it.

### Comment 24: What is the most important missing result?

**Defensible response**

There is no independent measured human HLA-A*02:01 affinity or natural-ligand
benchmark yet. The current quantitative results establish teacher fidelity,
not direct mechanistic accuracy.

### Comment 25: What single experiment would most strengthen the paper?

**Defensible response**

Prospectively select high, medium, and low student-scored peptides from held-out
de novo proteins and measure quantitative HLA-A*02:01 binding affinity. In
parallel, express the source proteins in a monoallelic HLA-A*02:01 system and
measure naturally presented ligands by immunopeptidomics.

## Fast judge questions

### What are you claiming in one sentence?

A frozen-ESM student provides fast local imitation of NetMHCpan HLA-A*02:01
binding and ligand rankings for inspectable de novo campaign triage.

### Did you beat NetMHCpan?

No. We reproduced it locally. Superiority requires an independent measured
human HLA affinity or ligand-presentation cohort.

### Why is this useful if it is only a surrogate?

It makes per-window screening local, sub-second, versioned, and easy to place
inside a generation loop without repeated network calls.

### What prevents leakage?

Parents sharing any exact 9-mer are assigned to the same fold. This blocks
direct peptide reuse across train and test.

### What remains vulnerable to leakage?

Near-neighbor peptides, structural-family similarity, and ESM-2 pretraining
overlap were not fully controlled.

### What is the score?

The hybrid composite is the geometric mean of predicted N-cleavage,
C-cleavage, and BA binding propensity. EL presentation and TAP are reported
separately. The protein summary is the mean of the five highest window
composites.

### Is it a probability?

No. It is an uncalibrated prioritization heuristic.

### Why only one allele?

Scope and data. It proves one coherent end-to-end lane, not population
coverage.

### Does BA include cleavage and presentation?

No. BA is the peptide-HLA binding event. Cleavage and TAP are separate upstream
predictions. EL is the presentation-oriented channel.

### What is your strongest evidence?

Component-grouped five-fold teacher fidelity, out-of-corpus fidelity on 95
independent designed peptides, and an inspectable runtime with preserved model
and artifact provenance.

### What is your weakest evidence?

Direct mechanism. There is no independent measured HLA-A*02:01 affinity or
monoallelic eluted-ligand benchmark.

### What would make you stop or abstain?

Unsupported allele, missing model artifact, structurally invalid candidate, or
a future distance-to-training threshold outside the calibrated domain.

## Publication-grade experiments still required

1. Quantitative HLA-A*02:01 peptide-binding measurements.
2. Monoallelic HLA-A*02:01 immunopeptidomics on full expressed constructs.
3. Teacher-training overlap audit for both measured cohorts.
4. Multi-allele training and population-weighted evaluation.
5. Near-neighbor and structural-family split sensitivity.
6. Simple peptide and smaller-encoder baselines.
7. Prospective operating-point calibration on independent data.
8. OOD distance and abstention evaluation.
9. Compute-equivalent latency against a local NetMHCpan installation.
10. Revalidation of Chao1 cleavage and TAP outputs.
11. One retained, approved end-to-end generation-to-screen agent trace.
