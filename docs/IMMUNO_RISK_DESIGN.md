# Immuno-risk design (Track C — late stage)

Screening / ranking scores for de novo binder candidates under **intracellular / plasmid delivery**.
Not a clinical immunogenicity probability.

## User problem

Upstream design produces candidate sequences (± structures). This stage answers:

> Given this sequence and delivery mode, which peptides are most likely presented,
> how foreign do they look relative to healthy-tissue ligands, and what is a
> polarized risk ranking we can defend to judges?

## Delivery route and MHC arms

| Arm | Why it matters for plasmid / cytosolic expression | MVP depth |
| --- | --- | --- |
| **MHC-I** | Primary: cytosolic antigen processing → CD8 surveillance of transduced cells | Full path + calibrated IEDB risk head |
| **MHC-II** | Secondary: APC expression / autophagy / uptake of released protein → CD4 help | Thin presentation baseline only |
| **Aggregation** | Separate physical-risk channel; can adjuvant responses without epitopes | Lightweight sequence/structure report |

Vector innate sensing (cGAS-STING, CpG, LNP/AAV) is **out of scope** — documented so the score is not over-read.

## Biological pathway (MHC-I)

1. Translation (± misfolded / rapidly degraded products as *contextual* risk, not assumed dominant).
2. Proteasome / immunoproteasome C-terminal generation.
3. TAP + ERAP trimming (captured indirectly by MHCflurry processing / presentation).
4. Peptide–MHC-I binding and surface presentation.
5. Tolerance vs foreignness vs T-cell recognition (distinct evidence layers).

## Biological pathway (MHC-II, thin)

Endolysosomal processing (cathepsins, AEP) → nested ~13–25mer ligand sets → MHC-II binding.
For v1 we **slide windows / predictor-native digestion** rather than insisting on recursive cathepsin simulation.
Cleavage tools remain diagnostic / interpretability, not the MHC-II score.

## Score semantics

| Score | Meaning | Not |
| --- | --- | --- |
| `mhcflurry_baseline_score` | Presentation propensity from MHCflurry | Immunogenicity probability |
| `iedb_risk_score` | Calibrated MHC-I head on explicit IEDB T-cell assay labels | De novo ground truth |
| MHC-II presentation ranks | Likelihood of class-II ligand presentation | ADA prediction |
| Aggregation report | Self-association / solubility / cytosolic cysteine flags | Epitope content |
| Confidence | Coverage, calibration, predictor agreement | LLM confidence |

Raw neural-net scores are never treated as probabilities. Prefer **percentile ranks** for cross-allele comparison.

## Evidence layers (keep separate)

1. **Presentation** — MHCflurry (± NetMHCpan-4.2e); HLAIIPred (± NetMHCIIpan-4.3k).
2. **T-cell assay labels** — IEDB positives/negatives only; never invent negatives from untested peptides.
3. **Self / benign ligands** — HLA Ligand Atlas (observed healthy-tissue ligands). Presence supports presentation of self, not absolute tolerance.
4. **Aggregation / persistence** — separate artifact; do **not** infer aggregation from low protease accessibility (relationship is non-monotonic).

## NetMHCpan literature (locked findings)

### NetMHCpan-4.1 (Reynisson et al., NAR 2020)

- Joint BA + EL training; NNAlign_MA motif deconvolution for multi-allelic MS data.
- Gains strongest for **ligand presentation**; epitope benchmark roughly at par with 4.0.
- Default SB/WB: %Rank < 0.5 / < 2. Prefer %Rank over raw score.

### NetMHCpan-4.2 (Nilsson et al., Front Immunol 2025; package **4.2e**)

- Extended BA/EL coverage; deletion composition + MHC interaction-frequency features for 10–14mers.
- Optional peptide **context** encoding helps EL prediction; does **not** clearly help CD8 epitope prediction.
- Fine-tuned heads on IEDB pathogen epitopes and CEDAR neoepitopes: **modest, domain-specific** gains.
  Pathogen-tuned models do not transfer to neoepitopes and vice versa.
- For de novo proteins, expose pathogen/neo heads only as **domain-mismatched sensitivity analyses**.
- Pin local binary to **4.2e** (context retrain + binding-core voting fix). Academic license; user-provided path only.

### NetMHCIIpan-4.3 (Nilsson et al., Sci Adv 2023; package **4.3k**)

- DR/DQ/DP coverage with inverted binders; context encoding option.
- Optional licensed comparator for the thin MHC-II arm; open baseline is **HLAIIPred**.

## Assumptions and limits

- Screening score for weekend demo / design filtering — not clinical go/no-go.
- MHC binding ≠ immunogenicity; Atlas ligand ≠ tolerance proven.
- ESM / foundation models may have seen UniProt; evaluate on de novo sequences and time-held-out assays where possible.
- Do not claim DRiPs are always the main MHC-I source; treat poor foldability as contextual evidence only.
- Do not collapse aggregation into epitope risk.

## Evaluation design

- **Sanity**: known epitopes (e.g. GILGFVFTL / HLA-A\*02:01) must rank above decoys.
- **Leakage-safe**: group by publication + source protein; exclude 8-mer overlap with NetMHCpan / MHCflurry training when provenance allows.
- **Metrics**: AUROC, AUPRC, Brier / calibration, top-k recall, runtime, coverage — reported **separately** for MHC-I, MHC-II, aggregation.
- **De novo**: no reliable labels; report uncertainty and structure proxies.

## Benchling role

System of record for candidates and run summaries (pull AA sequences → publish idempotent run records with provenance + artifact checksums + LangSmith URL). Detailed peptide tables stay local unless a tenant schema exists.

## Future work

- Full cathepsin recursive processing model with RSA / ΔG gates.
- Calibrated MHC-II T-cell head; population-weighted DRB3/4/5 + DQ/DP panels.
- MAPPs / lysosomal LC-MS labels on our own designs.
- Differentiable risk for design-loop deimmunization.
- Aggregation module with AGGRESCAN3D / CamSol / SEC-MALS validation.
