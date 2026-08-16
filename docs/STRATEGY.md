# Immunogenicity risk for de novo binders — strategy

**Status:** drafted Sat Aug 15, 2026 · 5:00 PM PDT. Submission Sun 10:45 AM.
**Supersedes:** ad-hoc plans in chat. Aligns with [`harness/TASK.md`](../harness/TASK.md) — see [Reconciling the two pipelines](#reconciling-the-two-pipelines).

---

## 1. Thesis

> Immunogenicity predictors are trained exclusively on natural proteins and are systematically
> out-of-distribution on de novo binders. Natural sequence is finite; de novo sequence is effectively
> unlimited but unlabeled. We close that gap with **consistency regularization on unlabeled de novo
> windows**, and ship a risk score, a calibrated confidence, and a per-residue heatmap so a designer
> can triage which designs to take into the wet lab.

The defense of the approach, in one line: designed proteins are **near-distribution, not
out-of-distribution** — de novo design works precisely because it exploits the latent structure
already present in models trained on natural sequence. That is exactly the regime where consistency
regularization is known to help.

What consistency regularization **cannot** do: invent knowledge about immune outcomes. It smooths
the decision function over the de novo manifold. Say this out loud in the writeup; overclaiming here
is the fastest way to lose a technically strong judge.

---

## 2. Locked scope decisions

| # | Decision | Rationale |
| --- | --- | --- |
| 1 | **MHC class II / CD4 T-cell axis** | Anti-drug antibodies run through CD4 help and class II presentation. Class I is for intracellular/plasmid delivery only. |
| 2 | **One risk axis this weekend: processing liability** | Epitopes generated → presented → not tolerated. Persistence/aggregation is deferred (§9). |
| 3 | **No pseudo-labels, no distillation from another predictor** | Mean Teacher / Π-model consistency only. The teacher is the student's own EMA. |
| 4 | **Headline metric is enrichment, not AUROC** | The product is a triage filter. See §6. |
| 5 | **One scoring path in the demo** | See [Reconciling the two pipelines](#reconciling-the-two-pipelines). |

### Reconciling the two pipelines

There are currently two overlapping implementations:

- **`harness/tools/*.ts`** — Pi harness, class I / proteasome oriented, with **fabricated numbers**:
  `mhc.ts` returns deterministic pseudo-ranks from character codes, and `tolerance.ts` scores
  peptides against a "self" list that is actually flu M1, CMV pp65, EBV BMLF1, NY-ESO-1, and HIV pol
  epitopes — foreign peptides labeled as tolerated.
- **The Python `immuno` package** (teammate's branch, not yet in this repo) — real trained model.

**Decision:** the Python model is the science. The harness either (a) becomes a thin inspectable
agent wrapper that calls it, or (b) is dropped from the demo. Do **not** ship both scoring paths — a
judge who reads the harness will find fake numbers next to real ones and discount both.

`harness/TASK.md` prioritizes class I; that is now superseded by decision 1.

---

## 3. Biological framing

```
de novo sequence (+ predicted structure)
        │
        ▼  endolysosomal processing — cathepsin S, L, B, AEP     ← accessibility: RSA, disorder, ΔG
        ▼  peptide pool — ragged nested sets (15-mer core + flanks)
        ▼  MHC-II presentation — DRB1 supertype panel
        ▼  tolerance check — self immunopeptidome, TCR-facing residues
        ▼  risk score — population-weighted
```

Notes that change implementation:

- Class II peptides are **ragged nested sets** around a 9-mer binding core, not clean 8–11mers.
  15-mer windows with flanks are the right unit.
- Processing is **endolysosomal**, not proteasomal. The ten simplified P1 motifs in
  `harness/tools/catalytic-sites.ts` (trypsin, chymotrypsin, caspase…) are the wrong protease set.
  Replace with cathepsin S/L/B and AEP specificity matrices (MEROPS + published proteomic
  specificity profiling) if the processing stage stays in scope.
- **Tolerance ≠ non-presentation.** Most self peptides are presented and tolerated. The self
  reference is the HLA Ligand Atlas, not a hand-typed list.

---

## 4. Label semantics — read before touching IEDB

This distinction has already caused one error and will cause more if not written down.

| IEDB table | `Positive` means | `Negative` means | Use |
| --- | --- | --- | --- |
| `mhc_search` | peptide was eluted off / bound that MHC | tested, did not bind | **Presentation.** Not immunogenicity. |
| `tcell_search` | T cells from a host responded (IFNγ/IL-2 release, proliferation) | tested, no response observed | **Immunogenicity.** This is our label. |

Verified counts (IEDB Query API, Aug 15 2026):

| Slice | Rows |
| --- | --- |
| `mhc_search` class I, Positive, 8–11mer | 3,241,985 |
| `mhc_search` class I, **Negative** | 124,932 |
| `tcell_search` class II, Positive (all variants) | **100,059** |
| `tcell_search` class II, **Negative** | **150,225** |
| `tcell_search` class II, human host | 172,256 (64,427 pos / 107,829 neg) |
| `mhc_search` class II, Positive | 1,419,132 |

**Why this matters.** The elution table is 26:1 positive-skewed and its negatives are peptides
someone had a reason to suspect — a biased hard-negative set, not a background sample. Training on
it would force us to invent decoys, and the student would inherit whatever bias that choice carried.
The T-cell class II table is 40/60 positive/negative with *real* measured negatives. This is the
single strongest argument for the label choice and belongs in the writeup.

Caveat to carry: a T-cell negative means "no response in these donors, with their HLA types, at this
assay's sensitivity." It is not proof of non-immunogenicity in the population.

`qualitative_measure` has variants — `Positive`, `Positive-High`, `Positive-Intermediate`,
`Positive-Low`, `Negative`. Filter with `like.Positive*`, not `eq.Positive`, or you undercount by ~14%.

---

## 5. Data inventory (verified)

| Source | Role | Scale | Access |
| --- | --- | --- | --- |
| **IEDB T-cell class II** | Teacher's labeled set | 172,256 human rows → ~160k 15-mer windows, ~46% pos, 10,933 protein groups (teammate's run) | Query API (`query-api.iedb.org/tcell_search`) or weekly `tcell_full_v3.zip` (1.3 GB) |
| **Protein Design Archive** | Unlabeled de novo pool **+ novelty stratification** | 1,963 designs → **1,509 unique sequences** → 272,280 15-mers (53,311 in 40–150 aa band) | Single 58 MB JSON: `backend/scripts/data.json` in [wells-wood-research/protein-design-archive](https://github.com/wells-wood-research/protein-design-archive) |
| **HLA Ligand Atlas** | Self / tolerance reference | ~12.7 MB gzipped, all five tables | `hla-ligand-atlas.org/rel/2020.12/{peptides,aggregated,sample_hits,protein_map,donors}.tsv.gz` |
| **RCSB "de novo protein"** | Alternative unlabeled pool | 2,408 entities (teammate's run) | RCSB search API |
| MegaScale (HF mirror) | **Stability only — deferred** | 17,093 sites / 365 domains / 117 designed | `LiteFold/MegaScale-Tsuboyama2023` |

### PDA is the upgrade, not just another pool

PDA and the RCSB query cover largely the same PDB universe, but PDA adds curation and **precomputed
novelty**: `struct_max_sim_natural` on 1,413 designs (range 0.324–1.000, median 0.839), plus
DE-STRESS metrics (packing density, hydrophobic fitness, solubility, charge, isoelectric point) on
~1,880 of 1,963.

That lets us bin the unlabeled pool by structural distance from anything natural and ask the
question that actually tests the thesis:

> **Does the consistency gain grow as designs get more novel?**

That is a far better plot than a single ablation bar, and no one else will have it.

Other verified PDA facts: lengths min 3 / p10 46 / median 151 / p90 340 / max 2,375; 620 unique
sequences (41%) in the 40–150 aa binder band; recency is good (254 designs from 2024, 270 from 2025,
157 from 2026), so the RFdiffusion-era distribution is represented.

### Why MegaScale is deferred

The HF mirror carries **no de novo sequences** — only 1,292 of 299,271 rows have a sequence at all,
spanning 10 *natural* PDB domains. Its trypsin/chymotrypsin ΔK50 columns correlate with ΔΔG at
**r = 0.97**, i.e. they are a stability readout, not a cleavage readout, so they cannot supervise the
processing stage. Real designed sequences require Zenodo record `7992926`
(`Processed_K50_dG_datasets.zip`, 1 GB; `AlphaFold_model_PDBs.zip`, 14 MB).

One finding from it is still worth a slide: designed domains are **10.3 pp more helical and 9.9 pp
less loop** than natural domains under consistent DSSP annotation. Loops are where proteases get
purchase, so that is a mechanistic reason designs may generate fewer epitopes — a quantitative
version of the qualitative claim in Gao et al.

---

## 6. Validation ladder

The mentor said it twice: the biggest risk is benchmark data — enough to train, enough held out.

| Rung | Test | Metric | Priority |
| --- | --- | --- | --- |
| **1** | Label-efficiency ablation on naturals: hide 90% of labels, does Mean Teacher beat supervised-only? | ΔAUROC / ΔAUPRC at matched label budget, **group split by source protein** | Must ship |
| **2** | Epitope recovery on held-out antigens: does the heatmap light up known positive windows? | Per-residue saliency AUC vs known epitope mask | Must ship |
| **3** | **Novelty-stratified gain**: bin PDA pool by `struct_max_sim_natural`, does the gain grow with novelty? | Δ metric per novelty bin | Must ship — this is the thesis |
| **4** | OOD behavior on de novo: stability and calibration of scores | Prediction variance under jitter, student–teacher agreement, reliability diagram | Should ship |
| **5** | De novo outcome ordering (Gao Table 3) | Qualitative rank check only | Nice to have |

### Headline framing

Report **enrichment / precision at top-k** and a **selective-prediction curve** (accuracy vs
coverage, thresholded on confidence). "If you can only express 8 of 200 designs, this picks 8 with
N-fold fewer predicted epitopes than random" is the same model as "test AUROC 0.695", described
honestly and usefully. The selective-prediction curve *is* the product.

### Rung 5 detail

From Gao et al. Table 3, roughly five informative points per direction:

- **Immunogenic:** IL-21 mimics (5/30 mice made anti-binder antibodies), SARS-CoV-2 RBD binder
  (1/10), anti-scaffold responses to RSV and SARS-2 nanoparticle vaccines, DARPin withdrawn from the
  clinic over ADA + intraocular inflammation despite excellent biophysics.
- **Clean:** Neo-2/15 (only human data point), IL-2/IL-15 mimics, D-protein binders.

Caveats to state: tiny n, heterogeneous assays, mostly murine (mouse MHC ≠ the human HLA we score).
Use for ordering, never as a metric.

---

## 7. Deliverables

1. **Risk score** — calibrated, reported as a percentile against a natural-protein reference cohort
   ("riskier than 80% of natural proteins"), not a raw logit.
2. **Confidence score** — from MC-dropout variance, student–teacher disagreement, and distance to
   the training manifold. Must honestly flag unfamiliar territory.
3. **Per-residue heatmap** — attention + gradient saliency, validated quantitatively by rung 2.

Outputs land in `results/immuno_risk/` as JSON + markdown so judges can re-run from inputs.

---

## 8. Workstreams and timeline

**Step zero, blocking everything:** the `immuno` implementation is **not in this repo**. Local clone
is still the scaffold plus harness (`src/re_agent/` has only `config.py` and `checks.py`, no
`.venv`, nothing since `ae799fd`). Get that branch pushed.

**Disk is a hard constraint.** 17 GiB free on a 99%-full volume; `~/.cache/huggingface` already
holds 59 GB; the teammate's embedding caches alone were 2.9 GB. `data/iedb_export/` is 27 GB and the
only file in it we need is the HLA Ligand Atlas submission (`1036766.xml`, 2.65 GB), which is
replaceable by a 12.7 MB download. **Delete `data/iedb_export/` before anything else.**

| When | Workstream | Owner |
| --- | --- | --- |
| Now → +1h | Push teammate branch; free disk; `./scripts/setup.sh` + ML extra; reproduce baseline run | Mihir + teammate |
| Now → +3h | PDA extraction: dedup to 1,509 sequences, 40–150 aa subset, novelty bins from `struct_max_sim_natural`, tile to 15-mers | Mihir |
| Now → +2h | HLA Ligand Atlas fetch (retry + gzip validation — host throws intermittent 502s and can take minutes) | Mihir |
| +1h → +6h | Rungs 1, 2, 4 on existing code | teammate |
| +3h → +7h | Rung 3 (novelty-stratified) once PDA bins land | joint |
| +7h → +9h | Demo artifact, README, figures, writeup | joint |
| Sleep | — | non-negotiable, ~4h |

---

## 9. Deferred, with reasons

| Item | Why not now |
| --- | --- |
| **Persistence / aggregation axis** | Real (a protein proteases won't touch persists, accumulates, aggregates) and it conflicts in sign with processing liability — Gao et al. argue high stability *lowers* immunogenicity by limiting antigen processing, while `risk.ts` adds +15 for `hyperstable_low_rsa`. Two axes, never summed. Needs ThermoMPNN-style stability prediction (MegaScale is its training set). **Next step, not this weekend.** |
| **FDA label ADA benchmark** | Overreaching. Labels themselves warn cross-product ADA comparison is misleading; assays differ in sensitivity and drug tolerance; route (SC vs IV), dose frequency, immunosuppression, and formulation-driven aggregation dominate and are invisible to a sequence model; most approved biologics are 150 kDa mAbs (wrong size class, germline-humanized); filtering to comparable non-mAb products leaves a handful; and it needs an unvalidated peptide→protein aggregation step, so a null result would be uninterpretable. Rung 2 gives the protein-level anchor with matched label semantics. |
| **Proto protease-interface check** | Expect it to be permissive. Proteases have broad specificity; structure predictors on an OOD interface can be confidently wrong; and cathepsin cleavage is substrate-threading-through-a-cleft, not rigid-body docking, so a predicted interface is not a predicted scissile bond. Cheaper substitute with more signal: cathepsin S/L/B + AEP specificity matrices from MEROPS. |
| **Class I / CD8 scoring** | Only relevant for intracellular or plasmid delivery. |

---

## 10. Risks

| Risk | Mitigation |
| --- | --- |
| Teammate's branch doesn't land | Everything downstream stalls. Chase it first, before any other work. |
| Disk fills mid-run | Delete `data/iedb_export/`, set `HF_HOME` off the full volume. |
| Semi-supervised arm shows no gain | Still shippable: report it honestly, with the novelty-stratified breakdown showing *where* it does and doesn't help. A clean negative result beats a vague positive one. |
| Heatmap is decorative | Rung 2 makes it a measurement. If rung 2 fails, present the heatmap as exploratory. |
| Scope creep into the persistence axis | It is explicitly deferred in §9. Point here. |

---

## 11. Open questions

1. Which DRB1 alleles for the population-weighted panel, and what population frequencies?
2. Does the processing/cleavage stage stay in the weekend scope at all, or does the model score
   sliding 15-mers directly? (Sliding windows keep teacher and student inputs identically
   distributed — a real argument for dropping the cleavage stage from the modeling path and keeping
   it as interpretation only.)
3. Reference cohort for percentile calibration: which natural proteome subset, and does it need
   whole-protein tilings rather than random windows?

---

## 12. References

- Gao et al. *De novo protein design: a transformative frontier in clinical protein applications.*
  J Transl Med 2026;24:319 — [`docs/de_novo_protein_design.pdf`](de_novo_protein_design.pdf).
  Immunogenicity §pp. 10–12; developability/CMC §p. 12; Table 3 §p. 11.
- Reynisson et al. *NetMHCpan-4.1 / NetMHCIIpan-4.0.* — [`docs/netmhc4-1.pdf`](netmhc4-1.pdf)
- Nilsson et al. *NetMHCpan-4.2.* — [`docs/netmhc4-2.pdf`](netmhc4-2.pdf)
- Chronowska, Stam, Wood. *The Protein Design Archive.* — [`docs/protein_design_archive.pdf`](protein_design_archive.pdf)
- Tsuboyama et al. *Mega-scale experimental analysis of protein folding stability.* Nature
  2023;620:434–444. Zenodo `10.5281/zenodo.7992926`.
