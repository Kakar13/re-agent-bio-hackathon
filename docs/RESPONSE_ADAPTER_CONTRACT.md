# Response-model adapter handoff

The agent does **not** train response models. Each teammate hands off either:

1. a callable implementing `ResponseModelAdapter`, or
2. an immutable JSON prediction artifact consumed by `ArtifactResponseAdapter`.

## Prediction artifact

```json
{
  "model_card": {
    "adapter_id": "team-model",
    "version": "checkpoint-sha-or-semver",
    "score_scale": "calibrated_probability",
    "calibration_id": "team-model-study-heldout-v1",
    "parameters": {}
  },
  "predictions_by_sha256": {
    "<sha256 of parent sequence>": [
      {
        "start": 0,
        "end": 15,
        "sequence": "ABCDEFGHIJKLMNO",
        "score": 0.73,
        "confidence": 0.81
      }
    ]
  }
}
```

If scores are not calibrated, set `score_scale` to `uncalibrated`. The agent will
show the adapter output but withhold late fusion.

## Optional frozen calibration artifact

```json
{
  "calibration_id": "team-model-study-heldout-v1",
  "adapter_id": "team-model",
  "adapter_version": "checkpoint-sha-or-semver",
  "raw_score_knots": [0.0, 0.3, 0.7, 1.0],
  "calibrated_probability_knots": [0.02, 0.18, 0.63, 0.91]
}
```

Calibration is fit outside the agent on the frozen validation split. The agent
only applies this versioned mapping. Model selection uses the frozen benchmark
under `data/processed/benchmarks/`; it does not retrain any adapter.

## Independent MHC-I processing surrogate

`models/chao1/cv5_heads.pkl 2` is not a response-model artifact. It enters through
`TeamE2EPLSAdapter` as an optional `mhc_i_processing_surrogate` lane. Its
HLA-A*02:01 9-mer cleavage, TAP, and MHC-I outputs are displayed and projected
onto residues, but they are never fused with NetMHCIIpan EL/BA or used to fill a
missing CD4 response score. See `models/chao1/MODEL_CARD.md` for the artifact hash,
scope, and smoke-test boundary.
