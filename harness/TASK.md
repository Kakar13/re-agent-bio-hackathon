# De novo immuno-risk task (harness = **last stage**)

**Status:** MHC-I existing-predictor MVP locked — see [`docs/IMMUNO_RISK_DESIGN.md`](../docs/IMMUNO_RISK_DESIGN.md).

The Pi harness owns the **late pipeline**, not backbone diffusion / MPNN. Upstream produces candidate sequences ± structures. Downstream here: **presentation → tolerance evidence → dual-arm risk + aggregation report**, with an inspectable trace for judges.

## Spec

| Field | Value |
| --- | --- |
| **Goal** | Given candidate binder sequence(s) ± structure and delivery mode, emit **MHC-I risk**, thin **MHC-II presentation**, **aggregation** report, and supporting evidence. |
| **Success criteria** | (1) End-to-end on ≥1 natural control + ≥1 de novo candidate. (2) Every flagged peptide has MHC call + tolerance/Atlas note + residue projection. (3) Judges can re-run from `../results/immuno_risk/<run-id>/`. (4) Benchling pull/publish when credentials exist. |
| **Constraints** | **MHC-I full** (MHCflurry required; NetMHCpan-4.2e optional). **MHC-II thin** (HLAIIPred required; NetMHCIIpan-4.3k optional). Human HLA. Intracellular/plasmid delivery is the default story. |
| **Out of scope** | Backbone diffusion, retraining MHCflurry/NetMHCpan, clinical ADA claims, vector innate sensing. |
| **Primary tools** | Python `re_agent.immuno_risk` backend + harness Pi tools + Benchling + Paperclip + LangSmith |
| **Demo artifact** | `manifest.json`, `peptides.csv`, `summary.json`, `residue-risk.json`, `aggregation.json`, `report.md` under `../results/immuno_risk/<run-id>/` |

## Pipeline

```text
sequence (+ optional structure)  [± Benchling AA Sequence]
        │
        ├─► MHC-I scan (MHCflurry ± NetMHCpan-4.2e) ──► IEDB-calibrated risk head
        ├─► MHC-II presentation (HLAIIPred ± NetMHCIIpan-4.3k)   [thin]
        ├─► HLA Ligand Atlas / self evidence join
        ├─► Aggregation / persistence report (separate score)
        └─► confidence + residue heatmap + artifacts (± Benchling publish)
```

Cleavage vs ~10 catalytic sites remains a **diagnostic** layer for interpretability; MHC-I peptide funnel is driven by MHCflurry processing/presentation, not the general protease motif count.

## Score contract

- Screening / ranking only — not clinical probability.
- Expose `mhcflurry_baseline_score` and `iedb_risk_score` side by side (never silent blend).
- MHC-II is presentation-only in this MVP.
- Aggregation is a separate report; do not infer it from protease accessibility.
- Confidence from calibration, allele coverage, predictor support, reference coverage, agreement — not LLM confidence.

## Evaluation

| Cohort | How |
| --- | --- |
| Known MHC-I epitopes | Sanity ordering (may overlap predictor training) |
| IEDB T-cell assays | Grouped holdout; explicit negatives only |
| Healthy self / Atlas | Foreignness reference |
| De novo | Proxy validation + uncertainty; never claim labels |

## Inputs

- Candidate FASTA / Benchling AA Sequence IDs
- Optional structure for aggregation / RSA proxies
- HLA panel (default common HLA-A/B/C + DRB1 panel)
- `.env` credentials: Anthropic, LangSmith, optional Benchling / NetMHC\* paths

## Agent working rules

- Orchestrate with inspectable tools + files; prefer structured JSON/CSV/MD.
- Fail closed when licensed predictors or Benchling are misconfigured — still finish local analysis.
- Cite Paperclip line URLs; record predictor versions beside each row.
- Ask before committing; never invent API keys.
