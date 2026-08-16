# re:AGENT project plan, write-up, demo script, and rubric

**Status:** Submission working document, 16 August 2026
**Track:** Track C, biological design to a specification, with a Track A style inspectable agent wrapper
**Project:** Design-to-screen agent for de novo protein campaign down-selection

This is the single source of truth for the project narrative, current evidence, demo video, judge questions, remaining gaps, and submission rubric. The scientific paper remains in `main.tex`; this document is the operational guide for presenting and finishing the project.

## 1. The project in 30 seconds

A standard RFdiffusion → ProteinMPNN → AlphaFold campaign can leave thousands of structurally plausible candidates. Folding and target contact do not say whether the protein will be processed into peptides that human MHC molecules present. Those presentation events are one immunogenicity framework: cytosolic MHC-I processing and HLA-A*02:01 presentation on the nucleated cell that made the protein, including immune cells. TCR recognition, danger sensing, and animal or human outcomes are later steps — they are why a construct that must function in a human still has to be tested in mice, and they are not what this model scores.

This MHC-I-only filter is the primary match when the designed protein is expressed inside the cell from a plasmid, mRNA, or viral vector. That includes intracellular binders and biosensors, not only secreted drugs. An injected extracellular protein can still reach MHC-I by cross-presentation; MHC-II, CD4 help, and antibodies remain the larger gap for that product class.

We built an inspectable design-to-screen agent that adds an early HLA-A*02:01 processing and presentation filter after structural triage. The agent grounds a target in literature, prepares typed design calls, gates GPU execution behind explicit approval, validates generated structures, and then scores every 9-mer in passing candidates. The output is a per-window liability profile used to remove the bulk of high-presentation candidates before wet-lab and animal work.

The current model is not a predictor of clinical immunogenicity. It asks three narrower questions: how the protein can be proteolytically processed into peptides, how likely those peptides are to bind HLA-A*02:01, and how likely they are to be presented to CD8 T cells. Danger recognition is the biological reason those questions matter. It is not a model output.

### One-sentence claim

> A frozen ESM-2 student reproduces NetMHCpan HLA-A*02:01 rankings well enough to downselect a de novo campaign for cytosolic MHC-I processing and presentation liabilities — the framework that matters most when a designed protein is expressed from a genetic payload inside a nucleated cell — while remaining separate from TCR recognition, danger, MHC-II/ADA, and animal or clinical immunogenicity.

## 2. Why this matters

The bottleneck is moving. Protein generation is becoming cheap, but experimental characterization remains slow and expensive. A scientist may generate 10,000 structurally gated candidates and express only a handful. Those proteins still need animal testing before human use. Existing structure and binding scores help answer whether a design is likely to fold or contact its target. They do not answer whether degradation of the protein may produce strongly presented human MHC peptides.

Our product decision is therefore down-selection, not diagnosis. A candidate with a flagged window is not declared immunogenic. It is removed from the bulk or moved into a higher-priority review bucket, where a scientist can consider mutation tolerance, direct binding assays, immunopeptidomics, MAPPs, or donor-cell experiments.

This framing sits between in silico generation and in vivo testing. It does not claim that sequence-only computation can replace mice or humans.

## 3. What the system does

```text
natural-language design objective
          |
          v
Paperclip literature and protein evidence
          |
          v
versioned design specification with cited hotspots
          |
          v
Proto campaign plan: RFdiffusion3 -> ProteinMPNN -> AlphaFold2
          |
          v
explicit LangGraph approval interrupt before GPU execution
          |
          v
structural gates: uniqueness, perplexity, pLDDT, RMSD, PAE, clashes
          |
          v
HLA-A*02:01 9-mer processing and presentation screen
          |
          v
per-window profile, screening flags, structure mapping, provenance
          |
          v
deterministic scientific review and campaign down-selection
```

### Agent architecture

1. **Research:** Paperclip supplies literature and protein evidence with citations.
2. **Tool discovery:** The graph discovers exact Proto schemas instead of inventing payloads.
3. **Design specification:** Target chain, hotspots, binder length, and budget are versioned and hashed.
4. **Campaign planning:** RFdiffusion3, ProteinMPNN, AlphaFold2, RMSD, and clash calls are materialized before execution.
5. **Human approval:** A LangGraph interrupt exposes the expected GPU call count and requires explicit approval.
6. **Structural validation:** Missing metrics fail closed. Only candidates with a `pass` status enter screening.
7. **Immunological profiling:** MHC-I, cleavage, TAP, and accessibility remain separate evidence lanes.
8. **Review:** Deterministic checks enforce citations, candidate counts, structural gates, null handling, and residue alignment.
9. **Workbench:** The UI shows the reasoning trace, campaign stages, model outputs, and structure-linked residue tracks.

### Current design-lane status

A separate design run produced 500 RFdiffusion3 binder backbones across five targets, and the PD-L1 positive-control campaign passed its contact gate. The fully approved LangGraph campaign still needs one low-budget live run that retains the final trace and generated manifests. This distinction must remain visible in the demo.

## 4. Data and machine-learning method

### 4.1 Corpus construction

The Protein Design Archive supplied 1,602 de novo parent proteins. PDA is a recent curated archive of synthetic designed proteins, mostly from the last several years. It is not the Protein Data Bank, and the Chronowska et al. "40 years" title refers to the protein-design field, not the age of the database. Every parent was tiled into all overlapping 9-mers before any sampling. NetMHCpan 4.1 then labeled each peptide for HLA-A*02:01 using separate EL presentation and BA binding channels.

- 290,937 peptide occurrences before deduplication
- 160,012 unique 9-mers
- 608 parent connected components
- Five folds balanced by unique-peptide count
- 40,125 shared peptides that would have leaked under a naive peptide split

Parents sharing any exact 9-mer were connected, and each whole connected component stayed in one fold. This prevents motif reuse across related designs from creating train-test leakage.

PDA is training data in the final model. It must not be described as an external or out-of-distribution benchmark. The OOD statement that is allowed is narrower: designed proteins are out of distribution for NetMHCpan's original natural-peptide training data, and measured MHC labels on de novo proteins are almost unavailable. Distilling NetMHCpan onto PDA ESM embeddings is how we generalize that teacher for de novo campaign scoring. After training, held-out PDA folds measure teacher imitation, not biological OOD. External teacher-fidelity evidence comes from the independent NY-ESO-1 design cohort. Direct human HLA-A*02:01 affinity and ligand-presentation validation remains outstanding.

### 4.2 Student model

The student uses frozen `esm2_t33_650M_UR50D` embeddings, mean-pooled over each 9-mer, followed by a compact 512-unit head. Five heads were trained so each connected-component fold was held out once for testing.

The deployed student predicts:

- NetMHCpan EL rank propensity
- NetMHCpan BA rank propensity
- NetMHCpan-like BA affinity, invertible to predicted IC50 in nM

This is knowledge distillation. The student learns to imitate NetMHCpan, not to infer immune responses directly. Frozen ESM-2 embeddings make the learner small and cheap, but do not prove that every de novo sequence lies within a safe representation domain.

### 4.3 Processing profile

Each 9-mer receives separate fields:

- N-terminal cleavage probability from the legacy chao1 checkpoint
- C-terminal cleavage probability from the legacy chao1 checkpoint
- TAP transport score and uncertainty from chao1
- EL presentation propensity from the NetMHCpan student
- BA binding propensity and predicted IC50 from the NetMHCpan student
- conventional strong, weak, or nonbinder class from predicted EL rank
- a separate recall-calibrated screening flag
- optional residue accessibility from a matched deposited structure

The inspectable processing composite for window `i` is:

```text
R_i = (p_N,i * p_C,i * p_BA,i)^(1/3)
```

Campaign ranking uses a separate pathway score. EL, BA, and affinity are not averaged: affinity is a monotone transform of BA, and EL already includes presentation information.

```text
S_i = EL^0.70 * (sqrt(p_N,i * p_C,i))^0.30
```

EL is the presentation endpoint. Proteasomal generation is a soft gate. The 0.70/0.30 split is a biological prior, not a fitted coefficient. TAP remains separately reported.

At protein level, the displayed summary is the mean of the five highest window-level pathway scores `S_i`. This is a prioritization heuristic, not an immune-response probability.

## 5. Results that can be defended

### 5.1 Held-out teacher imitation

| Result | EL | BA | Interpretation |
| --- | ---: | ---: | --- |
| Pooled out-of-fold Spearman | 0.943 | 0.948 | Strong rank fidelity to NetMHCpan |
| Weak-binder AUPRC | 0.892 | 0.859 | Good recovery under class imbalance |
| Fold-to-fold Spearman SD | about 0.004 | about 0.004 | Stable across component-grouped folds |

These are teacher-fidelity metrics. They are not evidence that the student predicts immunogenicity.

### 5.2 Legacy head replacement

On 159,038 identical PDA peptides, including 3,956 NetMHCpan strong binders:

| Model output | Spearman vs NetMHCpan EL | Strong-binder AUROC | Strong-binder average precision |
| --- | ---: | ---: | ---: |
| Legacy chao1 MHCflurry head | 0.354 | 0.731 | 0.060 |
| NetMHCpan student EL | 0.943 | 0.993 | 0.820 |

The new head is much better than the old MHCflurry-derived head at reproducing the chosen NetMHCpan teacher. It does not beat NetMHCpan itself and cannot use this comparison to claim biological superiority.

### 5.3 External de novo peptide teacher fidelity

The external test contains 95 independently designed NY-ESO-1 HLA-A*02:01 peptides.

- Student versus fresh NetMHCpan: EL Spearman 0.909
- Student versus fresh NetMHCpan: BA Spearman 0.898
- Weak-binder call agreement: 97.9 percent

This supports out-of-corpus teacher imitation. The cohort's T-cell activation labels are not used as the model benchmark because activation is downstream of peptide generation, HLA binding, presentation, and TCR recognition.

### 5.4 Screening operating point

The conventional predicted strong-binder cutoff recovered only 52.3 percent of teacher strong binders. For design-time triage, that miss rate is unacceptable. The deployed screening flag therefore uses a separate threshold selected for 95 percent teacher recall:

- Recall: 95.0 percent
- Precision: 43.9 percent
- Peptides flagged: 8,563 of 159,038
- Teacher strong binders recovered beyond the conventional rule: 1,688

The standard `binder_class` field is unchanged. The high-recall `screening_flag` is a separate product decision, not a redefinition of NetMHCpan terminology.

### 5.5 Structure and actionability

Among 186,485 windows with matched structure evidence:

- 80.2 percent of strong windows were buried
- 77.1 percent of weak windows were buried
- 67.7 percent of nonbinders were buried
- Only 22.2 percent of strong or weak windows were surface-exposed

This matters for redesign. A surface window may be easy to resample. A buried window requires a backbone-aware change or stability check. Accessibility is an editability field, not an immunogenicity field.

### 5.6 Latency

For a 156-residue design:

- Local student, CPU, end to end: 0.76 seconds
- NetMHCpan through the IEDB API: 10.34 seconds
- Measured ratio: 13.6 times faster

This compares local inference with a network service, not model compute with a local NetMHCpan binary. The valid value claim is local, sub-second, key-free scanning inside a generation loop.

## 6. Human mechanistic validation target

The primary biological benchmark must match the model output rather than a downstream whole immune response.

### BA endpoint

Use an independent human HLA-A*02:01 peptide cohort with quantitative measured IC50 or \(K_d\). Report Spearman correlation, error in log affinity, median fold error, binder AUROC, and binder AUPRC at a predeclared threshold. Exclude teacher-training overlap where provenance permits.

### EL and processing endpoint

Use monoallelic or confidently assigned human HLA-A*02:01 immunopeptidomics with source-protein context and explicit negative sampling. This tests whether peptides are naturally presented. It is the appropriate endpoint for evaluating whether cleavage, transport, and binding together identify presented ligands.

### Boundary

BA models only peptide-HLA binding. It does not model peptide generation. EL includes presentation-related information but does not establish T-cell recognition or danger. Animal outcomes, PBMC cytokines, T-cell activation, and anti-drug antibodies are outside the primary benchmark.

## 7. What we claim and what we do not claim

### We claim

1. The student reproduces NetMHCpan rankings on component-grouped held-out data.
2. Fidelity remains high on an independent de novo peptide cohort.
3. The new student is a better NetMHCpan surrogate than the legacy MHCflurry-derived head.
4. Local inference is fast enough for per-candidate inner-loop screening.
5. The agent preserves tool versions, inputs, hashes, intermediate artifacts, and claim boundaries.
6. The output is useful for campaign down-selection and deciding which liabilities deserve follow-up.

### We do not claim

1. The system predicts clinical immunogenicity or ADA.
2. The student beats NetMHCpan on accuracy.
3. PDA is an external OOD benchmark after being used for training.
4. Teacher agreement is equivalent to measured HLA binding or presentation.
5. HLA-A*02:01 represents population-level risk.
6. A predicted IC50 is an experimental binding measurement.
7. Buried residues are less likely to be presented.
8. The 13.6 times latency result compares equivalent local binaries.
9. The inherited cleavage and TAP heads have been revalidated on the new corpus.
10. The current score predicts T-cell recognition, danger, cytokines, antibodies, or whole immune response.

## 8. Demo video script

**Target length:** 4 minutes 15 seconds
**Format:** one narrator, screen recording of the workbench, with the paper and evidence artifacts ready in adjacent tabs

### 0:00-0:20, hook

**Visual:** Start on a campaign view containing many de novo candidates. Zoom to a ranked shortlist.

**Narration:**

> De novo protein design is becoming a generation problem at scale. We can create hundreds of plausible binders faster than we can decide which ones deserve wet-lab time. Folding and binding scores are common. Immune-liability screening is still late, slow, and disconnected from the design loop.

### 0:20-0:45, mission

**Visual:** Show the workbench title and design-to-screen stages.

**Narration:**

> We built an inspectable agent that moves immunological screening earlier. It grounds the design objective in literature, prepares and validates a protein-design campaign, and then scores every 9-mer for HLA-A*02:01 processing and presentation liabilities. The goal is down-selection, not clinical prediction.

### 0:45-1:15, agent architecture

**Visual:** Follow the pipeline from Paperclip to Proto to the approval interrupt and structural gates.

**Narration:**

> Paperclip supplies cited target evidence. The graph discovers the exact Proto tool schemas, writes a versioned design specification, and materializes RFdiffusion3, ProteinMPNN, and AlphaFold2 calls before running them. GPU execution pauses at an explicit human approval gate. Generated candidates must pass uniqueness, pLDDT, RMSD, PAE, perplexity, and clash checks before immunological screening.

### 1:15-1:50, model and data

**Visual:** Show the corpus-lineage figure, then the ESM-2 student block.

**Narration:**

> For the MHC lane, we labeled every 9-mer from 1,602 Protein Design Archive parents with NetMHCpan 4.1. Parents sharing any exact peptide stayed in one cross-validation fold, preventing more than forty thousand shared peptides from leaking across train and test. A small head over frozen ESM-2 embeddings learns separate EL presentation, BA binding, and affinity channels.

### 1:50-2:15, risk profile

**Visual:** Open one candidate heatmap and hover the top window.

**Narration:**

> Each window keeps cleavage, TAP, EL, BA, predicted IC50, binder class, confidence, and accessibility separate. The processing score combines N-cleavage, C-cleavage, and BA propensity. EL remains separate so we do not count presentation evidence twice. The result is a profile a designer can inspect residue by residue.

### 2:15-2:50, evidence

**Visual:** Show the result table and latency comparison.

**Narration:**

> On held-out component folds, the student reaches 0.943 Spearman on EL and 0.948 on BA against its NetMHCpan teacher. On an external de novo NY-ESO cohort, EL fidelity remains 0.909. The student does not beat NetMHCpan and we do not claim it does. Its value is local availability: a 156-residue protein scans in 0.76 seconds on CPU, compared with 10.34 seconds through the IEDB API.

### 2:50-3:30, exact validation boundary

**Visual:** Show BA, EL, cleavage, and downstream immune recognition as separate stages. Highlight BA and EL.

**Narration:**

> Our benchmark stops at the mechanism we model. BA asks whether a peptide binds human HLA-A*02:01. EL asks whether a peptide is likely to become a presented ligand and captures more of the upstream processing pathway. Cleavage and TAP remain separate predicted components. We do not benchmark against animal outcomes, PBMC cytokines, or whole immune responses because those add mechanisms the model was never designed to predict. The next direct validation is measured HLA-A*02:01 affinity and monoallelic human immunopeptidomics.

### 3:30-4:00, inspectability and impact

**Visual:** Show provenance, hashes, citations, and deterministic review results.

**Narration:**

> Every stage preserves its inputs, provider versions, hashes, outputs, and caveats. Missing evidence stays missing, and candidates that fail structural validation never receive a final campaign rank. This lets a scientist use the screen as a traceable filter rather than a black-box verdict.

### 4:00-4:15, close

**Visual:** Return to the shortlist and final project title.

**Narration:**

> The future of protein design is more candidates, not fewer. Our contribution is a fast, inspectable screen that helps scientists spend wet-lab time on the right designs and brings immune-liability reasoning into the generation loop.

## 9. Likely judge questions and concise answers

### Are you predicting immunogenicity?

No. We predict an HLA-A*02:01 MHC-I processing and presentation liability. Immunogenicity also depends on MHC-II/CD4 help, TCR recognition, B-cell epitopes, tolerance, aggregation, route, dose, formulation, and patient state. The product is a triage filter.

### Can your student beat NetMHCpan?

Not on the task it was distilled to imitate. NetMHCpan is the teacher, so the student's accuracy is teacher-bounded. The student improves local availability, latency, integration, and per-residue screening. Biological superiority would require an independent experimental cohort and a model trained on those outcomes.

### Why train on PDA and then call this out-of-distribution?

We do not call PDA OOD in the final evaluation because PDA is training data. The allowed OOD claim is about the teacher, not the student test set: NetMHCpan was trained mainly on natural peptides, measured MHC data on designed proteins are almost unavailable, and designed 9-mers are therefore out of distribution for that original corpus. Distilling NetMHCpan onto PDA ESM embeddings is a pseudo-de novo adaptation so a local student can score new campaign sequences in designed-protein representation space. Out-of-corpus teacher-fidelity evidence comes from the independent NY-ESO design cohort. The ESM-2 encoder was pretrained on natural sequence, but that fact alone does not make the final PDA evaluation OOD.

### How did you prevent leakage?

We connected any parents sharing an exact 9-mer and assigned whole connected components to folds. This prevented 40,125 reused peptides from crossing the train-test boundary.

### What exactly is the risk score?

There are two numbers. The inspectable processing composite is the geometric mean of N-cleavage, C-cleavage, and BA. Campaign ranking uses the pathway score `S = EL^0.70 * (sqrt(N * C))^0.30`. The protein summary is the mean of the top five pathway scores. BA, TAP, binder class, screening flag, confidence, and accessibility remain separate.

### Why is BA inside the composite instead of EL?

NetMHCpan EL already contains presentation-related information, including upstream effects. Combining EL with BA in one average would count the same MHC event twice. BA is the binding-only lane, so it is the cleaner third factor in the inspectable composite. Ranking instead uses EL as the presentation endpoint and cleavage as a soft gate.

### Does MHC binding affinity measure peptide generation or immune recognition?

No. BA measures the peptide-HLA binding event. Proteolytic generation and transport are separate upstream processes, and TCR recognition is downstream. EL is the closer presentation endpoint because its training evidence includes naturally eluted ligands.

### Why only HLA-A*02:01?

It gave us one coherent allele and enough teacher labels to finish a defensible end-to-end model. It is not population coverage. A production screen needs a multi-allele panel weighted by the intended population.

### Why MHC-I for therapeutic proteins when ADA often depends on MHC-II?

That is the largest biological scope limitation, and it depends on the product. MHC-I is on nearly all nucleated cells, not only non-immune cells. The current lane is the primary match for intracellular expression from a genetic payload — plasmid, mRNA, or viral vector — including intracellular binders and biosensors that the cell manufactures. An injected extracellular protein drug can still reach MHC-I by cross-presentation, but protein-therapeutic ADA requires MHC-II/CD4 and B-cell modeling, which we do not score. The orchestration preserves MHC-I and MHC-II as separate lanes so a future class-II model can be added without relabeling this output.

### Is the speedup fair?

It is an operational comparison between local CPU inference and the IEDB network API. It is not a compute-equivalent comparison with a local NetMHCpan binary. We report 13.6 times, the hardware, sequence length, and network caveat together.

### What is agentic about the system?

The graph researches the target, discovers external schemas, drafts a cited design specification, plans a campaign, pauses for approval, executes tools, validates structures, screens candidates, and runs deterministic review. Every step produces an inspectable artifact rather than hiding the workflow in one model response.

### What would actually validate the biology?

Direct HLA-A*02:01 peptide-binding assays validate BA. Monoallelic HLA-A*02:01 immunopeptidomics on expressed constructs validates presentation-related EL and the processing stack. Whole immune-response assays are outside the model target.

### What is the most important missing result?

We do not yet have an independent measured human HLA-A*02:01 affinity or eluted-ligand benchmark. Current quantitative results establish teacher fidelity, not direct mechanistic accuracy.

### Why use ESM-2 embeddings?

They provide a strong frozen representation that makes the task head small and training fast. This is a valid transfer-learning strategy, but the frozen representation does not remove domain shift. External cohorts, distance diagnostics, and abstention remain necessary.

### What happens when a tool or model is unavailable?

The evidence lane returns an explicit unavailable status. Missing structural metrics fail the candidate closed, and missing response evidence keeps the combined score null. The agent never fabricates a replacement value.

## 10. Gaps and next experiments

### Submission-critical gaps

1. **Record the final video:** Use the script above, keep it under 4 minutes 30 seconds, and show the workbench rather than slides alone.
2. **Capture one approved live LangGraph campaign:** Use a minimal generation budget, retain the approval payload, tool trace, manifests, and final review artifact.
3. **Align every displayed metric to the deployed v4 checkpoint:** Remove older v1 values or label them explicitly. The final headline values are EL 0.943/0.892 and BA 0.948/0.859 for Spearman/AUPRC.
4. **Keep the score definition consistent:** Composite uses N-cleavage, C-cleavage, and BA. EL and TAP remain separate.
5. **Show the limitations before judges ask:** Put the one-allele, teacher-bounded, non-immunogenicity statement on the result screen.
6. **Freeze the demo artifacts:** Record checkpoint hashes, corpus hash, external cohort source, threshold manifest, and code revision.
7. **Verify the final PDF:** Check figures, table widths, references, and the mechanistic validation section after the latest Sundial compile.

### Mechanistic validation gaps

1. Independent quantitative HLA-A*02:01 IC50 or \(K_d\) measurements for high, medium, and low predicted peptides.
2. Monoallelic HLA-A*02:01 immunopeptidomics on full expressed de novo constructs.
3. Source-protein-aware negative construction for EL evaluation.
4. Explicit teacher-training overlap audit for the measured affinity cohort.
5. Revalidation of the inherited cleavage and TAP heads against presentation-relevant evidence.

### Machine-learning gaps

1. Multi-allele training and population-weighted evaluation.
2. Distance-to-training-manifold and calibrated abstention on novel designs.
3. Direct comparison with a local NetMHCpan binary for model-only latency.
4. Revalidation or replacement of inherited cleavage and TAP heads.
5. Prospective external protein-level presentation evaluation, not only peptide-level teacher imitation.
6. Calibration against measured affinity and eluted-ligand outcomes rather than teacher recall.
7. Baselines that separate the contribution of ESM-2 from simpler peptide features.

### Product and agent gaps

1. One fully recorded low-budget generation-to-screen live run.
2. Stable artifact packaging so another judge can reproduce the run without hidden cloud state.
3. Multi-candidate shortlist controls in the workbench.
4. Exportable FASTA, CSV, JSON, and structure-annotated reports.
5. A scientist-reviewed mutation handoff that respects structural constraints without claiming automated deimmunization.

## 11. Working rubric

This is a team self-rubric based on the event's stated judging bar: cite or measure every claim, make reasoning inspectable, keep results reproducible, and demo one focused result worth trusting. It is not an organizer-published weighted score.

Legend: **PASS** is demo-ready, **PARTIAL** needs a named follow-up, **MISSING** blocks the claim.

| Rubric item | Status | Evidence | Final action |
| --- | --- | --- | --- |
| One focused problem and claim | PASS | Campaign down-selection through MHC-I liability triage | Keep the same opening and closing claim |
| Scientific framing | PASS | Presentation is separated from immunogenicity throughout | Put limitation on screen, not only in paper |
| Data provenance | PASS | PDA lineage, teacher version, hashes, manifests | Freeze final artifact index |
| Leakage control | PASS | Exact-9-mer connected-component folds | Show corpus-lineage figure |
| ML evaluation | PASS | Five-fold fidelity, old-head comparison, external cohort | Use deployed v4 values only |
| Independent mechanistic validation | MISSING | No independent measured HLA-A*02:01 BA or EL cohort yet | Add affinity and monoallelic ligand benchmarks |
| Human relevance | PARTIAL | Human HLA-A*02:01 teacher and allele-specific outputs | Validate on human HLA measurements, not whole immune responses |
| End-to-end agent | PARTIAL | Research, planning, interrupt, execution, review are implemented | Record one low-budget approved live run |
| Structural gating | PASS | Six fail-closed checks and candidate-count reconciliation | Show one passing and one blocked candidate |
| Inspectable reasoning | PASS | Workbench, artifacts, citations, deterministic review | Keep provenance panel visible in video |
| Reproducibility | PARTIAL | Scripts and manifests exist; some artifacts depend on Modal volumes | Freeze hashes and a small downloadable example |
| Sponsor-tool integration | PASS | Paperclip, Proto, Modal, Anthropic, Sundial | Name each tool only where it materially contributed |
| Quantitative result | PASS | Fidelity, threshold, accessibility, and latency measurements | Avoid unsupported accuracy claims |
| Product usefulness | PASS | High-recall screen and structure-aware editability | Demo shortlist decision, not only a heatmap |
| Demo readiness | PARTIAL | Workbench and paper are ready; video and live trace remain | Record, trim, and rehearse Q&A |
| Claim discipline | PASS | Null handling and explicit non-claims | Do not say "predicts immunogenicity" |

### Rubric verdict

The project is scientifically presentable now if the claim remains narrow. The two largest remaining risks are demo completeness and biological overstatement. A short live trace plus explicit limitation language moves the project from a strong model report to a credible end-to-end scientific agent.

## 12. Final submission checklist

### Must complete

- [ ] Record and trim the demo video.
- [ ] Run one explicitly approved low-budget LangGraph campaign.
- [ ] Retain the LangSmith trace and local hashed manifests.
- [ ] Confirm the final workbench loads on desktop and mobile.
- [ ] Confirm `main.pdf` compiles with the mechanistic validation section.
- [ ] Use only deployed v4 metrics in the video and rubric.
- [ ] Keep score wording consistent: inspectable composite is N-cleavage, C-cleavage, BA; campaign rank is EL^0.70 times generation^0.30.
- [ ] Put "not validated for immunogenicity prediction" in the spoken demo.
- [ ] Show at least one citation and one provenance record on screen.
- [ ] Show how a candidate is down-selected, not merely scored.

### If time remains

- [ ] Add the final approved run screenshot to the paper.
- [ ] Package one tiny reproducible example with sequence, checkpoint hash, and JSON output.
- [ ] Add a local-NetMHCpan latency result if a licensed binary is available.
- [ ] Add one slide separating BA, EL, cleavage/TAP, and downstream immune recognition.

## 13. Artifact map

| Artifact | Purpose |
| --- | --- |
| `main.tex` and `main.pdf` | Scientific manuscript |
| `figures/workbench.png` | Workbench architecture and result view |
| `docs/SUNDIAL_WRITEUP.md` | Full deployed-model methods and quantitative results |
| `docs/MHCI_NETMHCPAN_TRAINING.md` | Reproduction workflow for corpus and training |
| `docs/ORCHESTRATION.md` | Agent runtime, gates, and validation policy |
| `models/a0201-netmhcpan-pda-cv5-v4/checkpoint/` | Deployed five-fold student |
| `results/benchmarks/nyeso_a0201/` | Independent peptide cohort |
| `results/benchmarks/screening_calibration/` | High-recall threshold evidence |
| `results/benchmarks/epitope_accessibility/` | Structure and editability analysis |
| `results/benchmarks/screening_latency/` | Operational latency benchmark |
| Proposed measured HLA-A*02:01 affinity benchmark | Direct BA validation |
| Proposed monoallelic HLA-A*02:01 ligand benchmark | Direct EL and processing validation |

## 14. Core references

1. Reynisson et al. NetMHCpan-4.1 and NetMHCIIpan-4.0. *Nucleic Acids Research* (2020). DOI: `10.1093/nar/gkaa379`.
2. Lin et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* (2023). DOI: `10.1126/science.ade2574`.
3. Visani et al. T cell receptor specificity landscape revealed through de novo peptide design. *PNAS* (2025). DOI: `10.1073/pnas.2504783122`.
4. Butcher et al. De novo design of all-atom biomolecular interactions with RFdiffusion3. bioRxiv (2025). DOI: `10.1101/2025.09.18.676967`.
