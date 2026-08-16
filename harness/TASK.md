# De novo immuno-risk task (harness = **last stage**)

**Status:** drafted from team sync (Aug 15) — refine numbers/tools as Mark locks ground-truth sources.

The Pi harness owns orchestration. Upstream Proto produces candidate sequences and structures; this stage runs **response-model adapters → MHC-II EL/BA evidence → processing/tolerance proxies → transparent candidate ranking**, with an inspectable trace for judges.

Biochemically, natural and de novo binders are the same composition; the difference is immune familiarity. Labeled cleavage/tolerance data exists mainly for naturals; de novo has little/no labels — validate with structure predictors and natural holdouts, not by pretending de novo has ground truth.

## Spec

| Field | Value |
| --- | --- |
| **Goal** | Given candidate binder sequences ± structure, emit a computational triage rank and separate response, MHC-II EL, MHC-II BA, processing, and tolerance evidence. |
| **Success criteria** | (1) End-to-end run on ≥1 natural control and ≥1 de novo candidate. (2) Every value carries provider/version/input provenance. (3) EL and BA remain separate. (4) Missing providers remain missing—never synthetic. (5) Judges can re-run from inputs under `../results/`. |
| **Constraints** | MHC class II / CD4 axis; fixed versioned DR/DP/DQ panel; computational risk proxy only. |
| **Out of scope (harness)** | Training models, clinical prediction, and silent substitution of heuristics for unavailable providers. |
| **Primary tools** | **Harness custom tools** (`immuno_architecture_status`, `run_immuno_pipeline`), Paperclip for evidence, Proto for upstream design/structure |
| **Demo artifact** | Versioned JSON: input id → all adapter outputs → EL/BA disagreement → processing/tolerance → transparent rank under `../results/immuno_risk/`. |

## Pipeline (what the agent implements)

```text
candidate sequence (+ optional structure)
        ├── response-model adapter registry (teammate / ESM / non-ESM)
        ├── NetMHCIIpan EL + BA on fixed HLA-II panel
        ├── MixMHC2pred default second opinion (+ optional Graph-pMHC)
        ├── processing/accessibility evidence
        └── HLA-gated TCR-face self similarity
                         │
                         ▼
            separate proxy outputs + versioned late-fusion rank
```

### Stage notes

1. **Response models** — Plug checkpoints or immutable prediction artifacts into one contract. The harness does not train them.
2. **MHC evidence** — NetMHCIIpan is primary. Preserve EL presentation rank and BA binding rank/IC50 independently. MixMHC2pred is the default second opinion; Graph-pMHC is optional.
3. **Processing** — EL supports presentation but is not relabeled as direct cleavage. Add explicit cleavage/accessibility outputs only when their provider artifacts exist.
4. **Tolerance** — TCR-facing self similarity supports tolerance only when matched-self and query peptides share predicted HLA presentation.
5. **Ranking** — Apply the predeclared versioned fusion rule. Withhold the combined rank when required evidence is missing.
6. **Structure validation** — Proto/AlphaFold outputs are supporting designability evidence, not MHC ground truth.

## Evaluation / ground truth

| Cohort | How we evaluate |
| --- | --- |
| **Natural proteins / known splicing or cleavage events** | Hold out experimentally described events. Model must not be told the answer; ask whether it recovers the known products / sites. Decouple any fine-tuning from the test set. |
| **De novo candidates** | No reliable labels. Use natural-trained behavior + structure-based MHC/peptide checks as **proxy validation**; report uncertainty. |
| **Model adapters** | Compare already-produced checkpoints on the frozen species/study benchmark. The harness performs no fitting. |

## Inputs (fill paths as upstream lands)

- Candidate FASTA / sequences from the Proto handoff
- One or more versioned response-model prediction artifacts
- NetMHCIIpan 4.3 cache with separate EL and BA columns
- Optional MixMHC2pred / Graph-pMHC challenger caches
- Self-proteome feature table and optional shared-HLA tolerance artifact

## Agent working rules

- Orchestrate this stage with an inspectable agent loop (tools + files), not a silent one-shot script dump.
- Cite literature via Paperclip line URLs; measure MHC/structure outputs with versions, seeds, and paths beside each row.
- Prefer CPU / small runs first; GPU/Modal only when required and configured.
- If an input/provider is missing, emit a structured unavailable status and stop that evidence lane. Never manufacture a score.
