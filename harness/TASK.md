# De novo immuno-risk task (harness = **last stage**)

**Status:** drafted from team sync (Aug 15) — refine numbers/tools as Mark locks ground-truth sources.

The Pi harness in this folder owns the **late pipeline**, not backbone diffusion / MPNN sequence design. Upstream (elsewhere) produces candidate sequences + structures (e.g. after backbone diffusion → ProteinMPNN / “NTNN” filtering). Downstream here: **cleavage → peptide pool → MHC presentation → tolerance → risk score**, with an inspectable trace for judges.

Biochemically, natural and de novo binders are the same composition; the difference is immune familiarity. Labeled cleavage/tolerance data exists mainly for naturals; de novo has little/no labels — validate with structure predictors and natural holdouts, not by pretending de novo has ground truth.

## Spec

| Field | Value |
| --- | --- |
| **Goal** | Given candidate binder sequence(s) ± structure, emit an **immunogenicity / processing risk score** and the supporting peptide → MHC → tolerance evidence. |
| **Success criteria** | (1) End-to-end run on ≥1 natural control with known processing behavior and ≥1 de novo candidate. (2) Every flagged peptide has cleavage rationale + MHC1 (and optionally MHC2) call + tolerance note. (3) Judges can re-run from inputs under `../results/`. |
| **Constraints** | Prefer **~10 well-characterized protease / catalytic sites** first (not a full proteome of proteases). Focus **MHC class I** if only one path fits the weekend (intracellular / plasmid delivery); keep MHC II as a second score if time. Human HLA context unless TASK is updated. |
| **Out of scope (harness)** | Backbone diffusion, large-scale MPNN redesign, training new foundation models, full semi-supervised teacher–student training (ideas welcome in notes, not the demo bar). |
| **Primary tools** | **Harness custom tools** (`run_immuno_pipeline`, cleavage/MHC/tolerance/risk), Paperclip (lit + HLA refs), Proto/Boltz optional for structure checks |
| **Demo artifact** | Table or JSON: input id → cleavage peptides → MHC scores → tolerance → **risk score** + short markdown/findings under `../results/immuno_risk/`. |

## Pipeline (what the agent implements)

```text
sequence (+ structure / MPNN scorefile)
        │
        ▼
accessibility / disorder / ΔG-proxy features   ← RSA, loop length, unfolding proxy
        │
        ▼
cleavage vs known catalytic sites (start with ~10)  →  peptide pool
        │
        ▼
MHC I presentation (priority)  [+ MHC II if time]
        │
        ▼
tolerance / self check (e.g. HLA Ligand Atlas–style healthy refs)
        │
        ▼
risk score  (+ later: aggregation / “too protease-resistant” flags)
```

### Stage notes

1. **Structure features** — From structure ± MPNN score files: relative solvent accessibility (RSA), local disorder / loop length, unfolding ΔG proxy, secondary structure. Sequence-only specificity is comparatively “solved”; accessibility and dynamics are the gap the team cares about.
2. **Cleavage** — Combinatorial match of the chain against selected catalytic sites → predicted cut products (final peptide sequences). Keep cleavage prediction **last among feature builders**, then feed the peptide pool forward.
3. **MHC** — Class I for intracellular expression (plasmid / genetic delivery). Class II if administered as extracellular drug. Prefer **two scores** (I and II) when feasible.
4. **Tolerance** — Compare predicted ligands to healthy-tissue / self immunopeptidome references (team pointed at **HLA Ligand Atlas** and similar). Native-like → lower immune-flag risk; foreign-like → flag.
5. **Risk score** — Combine presentation + foreignness (+ optional aggregation / hyperstable “won’t clear” flags). Polarize clearly for the demo (high vs low risk), not a black-box logit.
6. **Optional structure validation** — For peptides the model flags, run Boltz/AF/ESM-style complex prediction on peptide + recognition complex as a **consistency check**, not as primary ground truth.

## Evaluation / ground truth

| Cohort | How we evaluate |
| --- | --- |
| **Natural proteins / known splicing or cleavage events** | Hold out experimentally described events. Model must not be told the answer; ask whether it recovers the known products / sites. Decouple any fine-tuning from the test set. |
| **De novo candidates** | No reliable labels. Use natural-trained behavior + structure-based MHC/peptide checks as **proxy validation**; report uncertainty. |
| **Semi-supervised / teacher–student** | Future enhancement (student aggressive on unlabeled de novo space; teacher EMA). **Not required** for Sunday demo unless a thin proof-of-concept fits. |

## Inputs (fill paths as upstream lands)

- Candidate FASTA / sequences: _TBD → `../data/…`_
- Structures or MPNN scorefiles (RSA / energy terms): _TBD_
- Seed list of ~10 protease / catalytic sites: _TBD (Mark / lit)_
- HLA / tolerance reference: HLA Ligand Atlas (and Paperclip-cited methods for tumor vs healthy immunopeptidomes)

## Agent working rules

- Orchestrate this **last stage** with an inspectable agent loop (tools + files), not a silent one-shot script dump.
- Cite literature via Paperclip line URLs; measure MHC/structure outputs with versions, seeds, and paths beside each row.
- Prefer CPU / small runs first; GPU/Modal only when required and configured.
- If upstream inputs are missing, scaffold the pipeline on a **tiny natural control** and document the expected de novo handoff — do not invent a full design campaign here.
