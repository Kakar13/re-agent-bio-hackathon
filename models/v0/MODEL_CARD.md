# Team E2E-PLS checkpoint v0

Status: archived. The active LangGraph and workbench screening profile uses chao1.

## Artifact

- File: `models/v0/cv5_heads.pkl`
- SHA-256: `c69b509a2a3e557dec04f304336e0c2108579231c92a44c0afa0862b90a4667c`
- Bundle format: `e2e_pls-heads-pickle-v1`
- Model version: `e2e_pls-heads-v1`
- Dataset version hash: `7d335d75e690ce2c`
- Encoder: `esm2_t33_650M_UR50D` (1,280-dimensional final-layer residue embeddings)
- Runtime adapter ID: `team-e2e-pls-v0`

The adapter statically audits the pickle globals, loads it through a restricted unpickler, and
validates every expected tensor and array shape. The only permitted constructors reconstruct Torch
tensors, NumPy arrays, and ordered dictionaries. This does not make arbitrary third-party pickle
files safe. Only the repository-owned team checkpoint should be supplied.

The separately named chao1 artifact currently loads to exactly the same parameters and calibration
values as v0. See `models/chao1/MODEL_CARD.md`.

## Scientific scope

This model is a single-allele MHC class I antigen-processing surrogate for `HLA-A*02:01`. It scores
overlapping 9-mers with separate outputs for:

- N- and C-terminal cleavage probability,
- relative TAP log-IC50 and bootstrap uncertainty,
- MHC-I presentation propensity,
- a geometric processing composite and decomposed confidence.

It is not an MHC class II model, CD4 response model, therapeutic immunogenicity probability, or
experimental measurement. The workbench preserves it as an independent evidence lane and never
adds it to NetMHCIIpan EL/BA late fusion.

## Validation boundary

The pickle contains weights and calibration breakpoints but no held-out metrics artifact, dataset
manifest, or training code snapshot. Results from this bundle are therefore smoke inferences, not
reproduced benchmark performance.

The default PDA smoke sequence is `pda:9s14:0` (PDB `9s14`, 99 residues, PDA novelty bin `novel`,
nearest-natural Foldseek TM-score 0.4728). Its exact parent hash and all 91 overlapping 9-mers are
absent from the supplied 5,000,613-row teammate dataset. That reduces direct sequence leakage for
this smoke case but does not establish biological generalization.

## Run

```bash
uv run python scripts/smoke_team_model_pda.py
```

The command runs the team surrogate and the existing MHC-II screening pipeline, writes a versioned
workbench artifact, and prints a compact result summary. The first run downloads the ESM-2 650M
weights if they are not already cached.
