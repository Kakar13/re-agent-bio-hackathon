# Expansion track — what we build *after* the core works

**Status:** drafted Sat Aug 15, 2026 · 5:05 PM PDT. Submission Sun 10:45 AM (~17.5 h).
**Companion to:** [`docs/STRATEGY.md`](STRATEGY.md) — the core plan. This doc is the *parallel* track.

---

## 0. The governing decision

> **Get something functional first.** Everything in this document is contingent on the core
> immunogenicity risk predictor training and producing a number. Nothing here merges into the demo
> until the gate in §2 is passed.

The reason to write it down anyway: training runs take wall-clock time we cannot compress, and we
will be waiting regardless. This is the list of what to work on *during* that wait, ordered so that
whatever finishes is additive rather than half-finished.

**The trap to avoid:** chasing something that trumps the core before the core exists. A working
mediocre predictor with honest evidence beats an impressive architecture diagram with no numbers.

---

## 1. Where things actually stand (5:05 PM Sat)

| Thing | State |
| --- | --- |
| Branch | `data-exploration` (main is behind) |
| Core `immuno` package | **Not in this repo.** Exists on a teammate's machine. Blocking. |
| `src/re_agent/` | `config.py`, `checks.py` only |
| `.venv` | Created 16:54, `proto-tools[mcp]` being wired |
| Pi harness | Present, class I oriented, **contains fabricated numbers** (see STRATEGY §2) |
| Mock/test data | 6 FASTA fixtures in `data/raw/immuno_tests/`: flu M1 window, SARS-2 spike (full + N200), human insulin, lysozyme, ubiquitin |
| Disk | 17 GiB free on a 99%-full volume; `data/iedb_export/` holds 27 GB of mostly-deletable XML |

The six fixtures are enough to stand up and smoke-test a training/scoring loop end to end. They are
**not** enough to train anything meaningful — they're plumbing tests, not a dataset.

### Judge signal so far

One judge was *semi*-impressed by the semi-supervised angle; the immunology domain framing landed
better. Read that as: **lead with the biological problem and the data asymmetry, not with Mean
Teacher.** The method is the answer to a problem the audience must feel first. See §5.

---

## 2. The gate

Expansion work merges into the demo only once **all** of these hold:

1. Training runs end to end on real labeled data and reports a test metric on a **group split by
   source protein** (no leakage).
2. The scoring path produces the three deliverables for at least one de novo sequence: risk score,
   confidence, per-residue heatmap.
3. Rung 1 of the validation ladder (label-efficiency ablation) has a number, positive or negative.

Until then, expansion work lives on its own branch and stays out of the demo narrative.

---

## 3. Immediate parallel work (do now, while training)

| # | Task | Why it's safe to do in parallel | Est. |
| --- | --- | --- | --- |
| P1 | **Stand up training infra on mock data** — wire the loop against the 6 fixtures so the pipeline is exercised before real data lands | Pure plumbing; de-risks the real run | 1 h |
| P2 | **Free disk** — delete `data/iedb_export/` (27 GB); only `1036766.xml` mattered and it's replaceable by a 12.7 MB Atlas download | Unblocks embedding caches (~3 GB) | 15 m |
| P3 | **PDA extraction** — 1,509 unique sequences, 40–150 aa subset, novelty bins from `struct_max_sim_natural`, tiled to 15-mers | Feeds rung 3; independent of model code | 2 h |
| P4 | **HLA Ligand Atlas fetch** — 5 tables, retry + gzip validation (host throws intermittent 502s, can take minutes) | Feeds tolerance + E1 | 45 m |
| P5 | **E1 manifold figure** (§4) once embeddings exist | Highest visual payoff per hour | 1.5 h |

---

## 4. Expansion candidates, with verdicts

### E1 — Latent manifold diagnostic ✅ **BUILD**

Embed three cohorts and plot them together: **human self peptides** (HLA Ligand Atlas),
**pathogen-derived epitopes** (IEDB viral + bacterial positives), and the **de novo pool** (PDA).
Color by cohort, reduce with UMAP or PCA.

This is the single best expansion-per-hour on the list because it *visualizes the thesis instead of
asserting it*: where do designed proteins actually sit relative to self and to known-immunogenic
foreign? It pairs directly with rung 3 in STRATEGY §6 — that experiment shows the semi-supervised
gain grows with distance from natural; this figure shows what that distance looks like.

Cost is near zero: embeddings already exist for the labeled and de novo sets, and the Atlas is
12.7 MB. This is also the honest, defensible version of the "human-safe space" intuition — the
latent geometry is a **diagnostic**, not something we steer through (see E4).

### E2 — Cleavage head on elution-derived termini ⚠️ **CONDITIONAL STRETCH**

IEDB elution records carry `StartingPosition`, `EndingPosition`, and a UniProt accession, so real
cleavage boundaries are recoverable — genuine supervision, not a relabeling of the presentation
signal.

Two cautions. It drags class I data back in (MS immunopeptidomics is HLA-I dominated), which cuts
against STRATEGY decision 1. And eluted ligands are *already* the joint product of cleavage,
transport, and binding, so heads trained on the same records don't decompose the process — they
factorize one signal three ways and invite us to believe the factorization.

Only attempt if the core is done early and someone wants a second head.

### E3 — Three-head decomposition (cleavage / chaperone / presentation) ❌ **DEFER**

The proposal is a shared vector feeding cleavage + TAP + MHC heads. That is a modern rewrite of the
**class I** antigen processing chain — the published version is NetChop (cleavage) + a TAP
predictor + NetMHCpan, integrated as NetCTLpan. Coherent, and familiar.

**But TAP is class I only.** TAP transports cytosolic peptides into the ER for MHC-I loading. Class
II peptides never touch it — they are generated in the endolysosome and loaded there, with the
invariant chain / CLIP occupying the groove and **HLA-DM** performing peptide editing. Our diagram
(cathepsin S/L/B, AEP, DRB1 panel) is unambiguously class II, so a TAP head is biologically
incoherent in it.

The correct mapping:

| Head | Class I version | Class II version (ours) |
| --- | --- | --- |
| Cleavage | proteasome / immunoproteasome | cathepsin S, L, B, AEP |
| Chaperone | TAP transport | **HLA-DM editing** (DM-resistant peptides survive to the surface) |
| Presentation | MHC-I binding | MHC-II binding |

HLA-DM susceptibility is arguably *the* determinant of which class II peptides reach the surface and
is genuinely underexploited — but there is no public dataset at IEDB scale to train it. So the
three-head design either pushes us back to class I (wrong axis for anti-drug antibodies) or leaves
the middle head unlabeled.

Cost side: three label sets, three loss weights, and a consistency term defined per head. Too much
new surface for the time remaining.

### E4 — Deimmunization / "human-safe space" steering ❌ **REJECTED for this weekend**

Team call, and I agree: moving a design toward human-safe space and then *rescuing affinity
post-hoc* is not the outcome we are trying to demonstrate. Our job is triage — telling a designer
which candidates to spend wet-lab budget on — not redesign.

The technical reasons it also isn't buildable now, recorded so nobody re-proposes it at 1 AM:

- **It can't be a latent walk.** The architecture is frozen ESM-2 + an attention-pooling head — an
  encoder with no decoder. ESM embeddings aren't invertible; you cannot move a point and read a
  sequence back out.
- **The implementable version is discrete search**, not latent steering: take the residues the
  heatmap flags, propose substitutions (ProteinMPNN or ESM pseudo-likelihood), re-score, keep edits
  that lower predicted immunogenicity. That's the standard deimmunization workflow.
- **Binding loss is the entire difficulty.** Constraints needed: edits restricted to positions off
  the binding interface (requires structure + interface annotation), preference for high solvent
  accessibility, a hard cap on edit count, and re-scoring every candidate with binding/stability
  models (AF2 or Boltz interface metrics, ProteinMPNN log-likelihood, ThermoMPNN ΔΔG). It is a
  Pareto problem over immunogenicity × affinity × stability.
- **None of it is validatable in a weekend.** If shown at all, show it as unvalidated capability:
  five designs, three flagged positions each, proposed substitutions that hold ProteinMPNN
  likelihood and predicted contacts roughly constant.

On "thinking in the latent space" generally: representation-space *diagnostics* (E1) are real, cheap
and defensible. Latent-space *optimization* is real in protein design but requires a generative
decoder we don't have. Latent steering specifically for immunogenicity has no established precedent
— it would be a research programme, not a weekend feature.

### E5 — Agent loop / harness orchestration ⚖️ **DIFFERENTIATE, DON'T COMPETE**

Observation from the room: many teams are building agent loops that "do the whole science." That is
the modal project here, which means it is the hardest place to stand out and the easiest place to
look thin.

Our differentiator is the opposite of an agent loop: **a trained model plus real evidence that it
generalizes to a distribution it was never trained on.** The harness should therefore be a thin,
inspectable wrapper that calls the model and writes a legible trace — not a competing scoring path.

Concretely, per STRATEGY §2: delete the fabricated `mhc.ts` pseudo-ranks and the mislabeled
`tolerance.ts` "self" list (it contains flu, CMV, EBV, NY-ESO-1, and HIV epitopes scored as
tolerated), or drop the harness from the demo entirely. Fake numbers next to real ones discredit
both.

### E6 — Persistence / aggregation axis ⏸️ **DEFERRED** (already in STRATEGY §9)

Real: a protein proteases won't touch persists, accumulates, and can aggregate — and aggregation is
a named driver of immune activation. It conflicts in *sign* with processing liability, so the two
must never be summed. Needs ThermoMPNN-style stability prediction, with MegaScale as its training
set. Post-hackathon.

---

## 5. Narrative split for the demo

The judge signal says the immunology framing lands harder than the method. Structure accordingly.

**Arc (target ~5 min):**

1. **The problem, biologically.** You designed 200 binders. Some fraction will provoke an immune
   response in a patient and you cannot tell which. Existing pre-screens — Rosetta `ddg`, contact
   molecular surface, AF2 pLDDT/PAE — all answer "does it fold and bind." None answers "will the
   immune system see it." *(Gao et al. 2026 makes exactly this argument; Table 3 has the failure
   cases: a DARPin withdrawn over anti-drug antibodies despite excellent biophysics.)*
2. **Why you can't just use an existing predictor.** Every immunogenicity model is trained on
   natural proteins. Show the E1 manifold: here is self, here is pathogen, here is where designed
   proteins actually land.
3. **The data asymmetry.** Natural sequence is finite and labeled. De novo sequence is effectively
   unlimited and unlabeled. *This is the setup that makes the method inevitable rather than clever.*
4. **The method, briefly.** Consistency between a student and its own EMA teacher on unlabeled de
   novo windows. No pseudo-labels, no distilling another predictor. One slide.
5. **The evidence.** Rung 1 (label efficiency), rung 2 (heatmap recovers known epitopes), rung 3
   (gain grows with novelty). Framed as **enrichment at top-k**, not AUROC.
6. **The honest limits.** Consistency regularization smooths a decision function; it does not create
   knowledge about immune outcomes. No de novo protein has a ground-truth immunogenicity label.
   Route, formulation, and patient status are invisible to a sequence model.

**Speaker split:** immunology/problem framing and evidence interpretation from the biology side;
data + method + validation design from the ML side. The handoff should be at step 3, where the data
asymmetry motivates the method.

**What not to do:** open with "we used Mean Teacher." That's step 4 of 6.

---

## 6. Decision log

| Decision | Call | Reason |
| --- | --- | --- |
| Functional core before expansion | **Locked** | A working predictor with honest evidence beats an architecture diagram |
| Deimmunization + affinity rescue | **Rejected** | Not the outcome we're demonstrating; also needs a decoder we don't have |
| TAP / three-head decomposition | **Deferred** | TAP is class I; our axis is class II; the class II chaperone analog (HLA-DM) has no training data |
| Latent-space steering | **Rejected** | Encoder-only architecture; no precedent for immunogenicity |
| Latent-space *diagnostic* (E1) | **Build** | Cheap, visualizes the thesis, defensible |
| FDA label ADA benchmark | **Rejected** | Cross-product ADA rates aren't comparable; confounders dominate (STRATEGY §9) |
| Persistence/aggregation axis | **Deferred** | Conflicts in sign with processing; needs stability model |
| Harness as competing scoring path | **Rejected** | Fabricated numbers next to real ones discredit both |
| Training infra on mock data now | **Do it** | Plumbing de-risk while real data lands |

---

## 7. Open questions

1. Who owns getting the `immuno` branch pushed, and by when? Everything downstream is blocked on it.
2. If rung 1 shows **no** semi-supervised gain, do we present the negative result (recommended) or
   pivot the narrative to the data/manifold contribution?
3. Does the demo include the Pi harness at all, or is the artifact a Python CLI + `results/` JSON?
4. Which DRB1 alleles and population frequencies for the population-weighted score?
