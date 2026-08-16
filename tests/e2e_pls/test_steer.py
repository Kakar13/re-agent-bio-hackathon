import json

from re_agent.e2e_pls import score as score_mod
from re_agent.e2e_pls.steer import SteeringConfig, steer_to_safety

SEQ = "GSHMSEEELKEAVKLLKKAEELVKKGD"


def _highest_risk_window(heads, client):
    scores = score_mod.score_sequence(SEQ, "HLA-A*02:01", client, heads)
    return score_mod.aggregate_protein_risk(scores).max_risk_window


def test_steer_respects_max_mutations(trained_heads, mock_client):
    w = _highest_risk_window(trained_heads, mock_client)
    config = SteeringConfig(
        hla_allele="HLA-A*02:01", max_mutations=2, max_steps=6, mask_candidates_per_step=6
    )
    trace = steer_to_safety(SEQ, w.start, w.end, mock_client, trained_heads, config)
    assert len(trace.mutations) <= 2


def test_steer_never_touches_protected_positions(trained_heads, mock_client):
    w = _highest_risk_window(trained_heads, mock_client)
    protected = frozenset(range(w.start, w.end))  # protect the entire target window
    config = SteeringConfig(
        hla_allele="HLA-A*02:01", max_mutations=3, protected_positions=protected
    )
    trace = steer_to_safety(SEQ, w.start, w.end, mock_client, trained_heads, config)
    assert trace.mutations == []
    assert trace.output_sequence == trace.input_sequence


def test_steer_partial_protection_avoids_protected_positions(trained_heads, mock_client):
    w = _highest_risk_window(trained_heads, mock_client)
    protected = frozenset([w.start, w.start + 1, w.start + 2])
    config = SteeringConfig(
        hla_allele="HLA-A*02:01", max_mutations=3, protected_positions=protected
    )
    trace = steer_to_safety(SEQ, w.start, w.end, mock_client, trained_heads, config)
    mutated_positions = {m["position"] for m in trace.mutations}
    assert mutated_positions.isdisjoint(protected)


def test_steer_only_accepts_improving_mutations(trained_heads, mock_client):
    w = _highest_risk_window(trained_heads, mock_client)
    config = SteeringConfig(hla_allele="HLA-A*02:01", max_mutations=3)
    trace = steer_to_safety(SEQ, w.start, w.end, mock_client, trained_heads, config)
    assert trace.final_score["composite_risk"] <= trace.initial_score["composite_risk"]
    for step in trace.steps:
        for candidate in step.candidates:
            if candidate.accepted:
                assert candidate.composite_risk < trace.initial_score["composite_risk"]


def test_steer_trace_is_json_serializable(trained_heads, mock_client):
    w = _highest_risk_window(trained_heads, mock_client)
    config = SteeringConfig(hla_allele="HLA-A*02:01", max_mutations=3)
    trace = steer_to_safety(SEQ, w.start, w.end, mock_client, trained_heads, config)
    payload = json.dumps(trace.to_dict())
    assert "claim_disclaimer" in payload
    assert trace.claim_disclaimer in payload


def test_steer_output_sequence_length_unchanged(trained_heads, mock_client):
    w = _highest_risk_window(trained_heads, mock_client)
    trace = steer_to_safety(
        SEQ, w.start, w.end, mock_client, trained_heads, SteeringConfig(hla_allele="HLA-A*02:01")
    )
    assert len(trace.output_sequence) == len(SEQ)
