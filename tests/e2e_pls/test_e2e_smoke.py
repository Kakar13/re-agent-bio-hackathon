"""Deterministic end-to-end smoke test: score -> steer -> trace, all in mock mode."""

import time

from re_agent.e2e_pls import fixtures
from re_agent.e2e_pls import score as score_mod
from re_agent.e2e_pls.encoder import ProteinEncoder
from re_agent.e2e_pls.steer import SteeringConfig, steer_to_safety

DEMO_SEQ = fixtures.DEMO_SEQUENCES["synthetic_designed_binder_demo"]["sequence"]


def test_scoring_is_deterministic_across_fresh_clients(trained_heads):
    client_a = ProteinEncoder(mode="mock")
    client_b = ProteinEncoder(mode="mock")
    scores_a = score_mod.score_sequence(DEMO_SEQ, "HLA-A*02:01", client_a, trained_heads)
    scores_b = score_mod.score_sequence(DEMO_SEQ, "HLA-A*02:01", client_b, trained_heads)
    assert [s.to_dict() for s in scores_a] == [s.to_dict() for s in scores_b]


def test_full_pipeline_smoke_and_runtime_budget(trained_heads, mock_client):
    t0 = time.monotonic()
    scores = score_mod.score_sequence(DEMO_SEQ, "HLA-A*02:01", mock_client, trained_heads)
    risk = score_mod.aggregate_protein_risk(scores)
    w = risk.max_risk_window
    trace = steer_to_safety(
        DEMO_SEQ,
        w.start,
        w.end,
        mock_client,
        trained_heads,
        SteeringConfig(hla_allele="HLA-A*02:01"),
    )
    elapsed = time.monotonic() - t0

    assert len(scores) == len(DEMO_SEQ) - 9 + 1
    assert all(0 <= s.composite_risk <= 1 for s in scores)
    assert len(trace.mutations) <= 3
    assert (
        trace.output_sequence[: trace.target_window_start]
        == trace.input_sequence[: trace.target_window_start]
    )
    assert elapsed < 15, "mock-mode pipeline should stay fast enough for interactive dashboard use"
