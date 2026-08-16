# Team E2E-PLS checkpoint: chao1

## Artifact

- File: `models/chao1/cv5_heads.pkl 2`
- SHA-256: `a6a0e265e8466b08efb30f5a508eeeba954e205da66222efabc246c5a1dca41f`
- Bundle format: `e2e_pls-heads-pickle-v1`
- Model version: `e2e_pls-heads-v1`
- Dataset version hash: `7d335d75e690ce2c`
- Encoder: `esm2_t33_650M_UR50D`
- Runtime adapter ID: `team-e2e-pls-chao1`

## Runtime role

Chao1 is the primary custom model selected by the LangGraph screening profile. The restricted
loader identifies this exact artifact by its SHA-256 and rejects incompatible checkpoint shapes
before inference.

## Scientific scope

Like v0, chao1 is an HLA-A*02:01 MHC class I antigen-processing surrogate over 9-mers. It is not an
MHC class II model, CD4 response model, measured immune-response probability, or replacement for
NetMHCIIpan EL/BA. The pipeline displays it as a separate evidence lane and excludes it from MHC-II
late fusion.
