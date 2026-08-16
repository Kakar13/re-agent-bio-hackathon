from __future__ import annotations

from langchain_core.messages import HumanMessage

from re_agent.agent.graph import _agent_node, _compute_approval_summary, _review_artifact


def test_reviewer_enforces_placeholder_response_gating() -> None:
    artifact = {
        "id": "test-run",
        "payload": {
            "citations": [{"url": "https://example.test"}],
            "assessment": {
                "mhc_results": [
                    {
                        "provider_id": "netmhciipan",
                        "status": "ok",
                        "supported_alleles": [f"HLA-{index}" for index in range(18)],
                    }
                ],
                "response_results": [{"status": "unavailable"}],
                "combined_rank_score": None,
            },
        },
    }

    review = _review_artifact(artifact)

    assert review["status"] == "pass"


def test_keyless_team_model_screen_attaches_separate_checkpoint(
    monkeypatch,
) -> None:
    monkeypatch.setenv("RE_AGENT_FORCE_KEYLESS", "true")
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    update = _agent_node(
        {
            "messages": [HumanMessage(content=f"Screen candidate {sequence}")],
            "artifacts": [],
            "reviews": [],
            "screening_profile": "mhc_ii_plus_chao1",
        }
    )

    tool_call = update["messages"][0].tool_calls[0]
    assert tool_call["name"] == "screen_candidate"
    assert tool_call["args"]["sequence"] == sequence
    assert tool_call["args"]["mhci_surrogate_checkpoints"] == [
        "models/chao1/cv5_heads.pkl 2"
    ]


def test_keyless_standard_profile_does_not_attach_custom_model(monkeypatch) -> None:
    monkeypatch.setenv("RE_AGENT_FORCE_KEYLESS", "true")

    update = _agent_node(
        {
            "messages": [HumanMessage(content="Screen candidate ACDEFGHIKLMNPQRSTVWY")],
            "artifacts": [],
            "reviews": [],
            "screening_profile": "mhc_ii_standard",
        }
    )

    tool_call = update["messages"][0].tool_calls[0]
    assert "mhci_surrogate_checkpoints" not in tool_call["args"]


def test_direct_chao1_screen_bypasses_model_tool_selection(monkeypatch) -> None:
    monkeypatch.setenv("ANTHROPIC_API_KEY", "not-used-by-direct-route")
    monkeypatch.delenv("RE_AGENT_FORCE_KEYLESS", raising=False)
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    update = _agent_node(
        {
            "messages": [HumanMessage(content="Run visual screening")],
            "artifacts": [],
            "reviews": [],
            "screening_profile": "mhc_ii_standard",
            "direct_screen_request": {
                "sequence": sequence,
                "candidate_id": "pda:9s14:0",
            },
        }
    )

    assert update["direct_screen_request"] is None
    tool_call = update["messages"][0].tool_calls[0]
    assert tool_call["name"] == "screen_candidate"
    assert tool_call["args"]["sequence"] == sequence
    assert tool_call["args"]["mhci_surrogate_checkpoints"] == [
        "models/chao1/cv5_heads.pkl 2"
    ]
    assert tool_call["args"]["source_metadata"]["ui_action"] == "direct_chao1_screen"
    assert tool_call["args"]["structure_pdb_id"] == "9S14"
    assert tool_call["args"]["structure_chain_id"] == "A"


def test_reviewer_rejects_rank_when_response_is_unavailable() -> None:
    artifact = {
        "id": "test-run",
        "payload": {
            "citations": [{"url": "https://example.test"}],
            "assessment": {
                "mhc_results": [
                    {
                        "provider_id": "netmhciipan",
                        "status": "ok",
                        "supported_alleles": [f"HLA-{index}" for index in range(18)],
                    }
                ],
                "response_results": [{"status": "unavailable"}],
                "combined_rank_score": 0.7,
            },
        },
    }

    review = _review_artifact(artifact)

    assert review["status"] == "fail"
    assert any(
        row["name"] == "placeholder_response_gating" and not row["passed"]
        for row in review["checks"]
    )


def test_reviewer_keeps_mhci_surrogate_outside_netmhciipan_fusion() -> None:
    artifact = {
        "id": "test-run",
        "payload": {
            "citations": [{"url": "https://example.test"}],
            "assessment": {
                "sequence": "ACDEFGHIK",
                "mhc_results": [
                    {
                        "provider_id": "netmhciipan",
                        "status": "ok",
                        "supported_alleles": [f"HLA-{index}" for index in range(18)],
                    }
                ],
                "response_results": [{"status": "unavailable"}],
                "combined_rank_score": None,
                "component_summary": {
                    "proxy_values": {
                        "netmhciipan_el": 0.8,
                        "netmhciipan_ba": 0.7,
                    }
                },
                "mhc_i_surrogate_results": [
                    {
                        "provenance": {"capability": "mhc_i_processing_surrogate"},
                        "spatial_tracks": {"mhci_processing_risk_max": [0.1] * 9},
                    }
                ],
            },
        },
    }

    review = _review_artifact(artifact)

    assert review["status"] == "pass"
    assert any(
        row["name"] == "mhci_surrogate_separate" and row["passed"]
        for row in review["checks"]
    )


def test_reviewer_rejects_structurally_blocked_candidate_marked_screened() -> None:
    artifact = {
        "id": "campaign-run",
        "payload": {
            "citations": [{"url": "https://example.test"}],
            "manifest": {
                "candidate_counts": {"screened": 1},
                "candidates": [
                    {
                        "candidate_id": "candidate-1",
                        "validation_status": "fail",
                        "screening_status": "screened",
                        "validation_checks": [{"name": "plddt", "passed": False}],
                    }
                ],
            },
            "assessments": [{}],
        },
    }

    review = _review_artifact(artifact)

    assert review["status"] == "fail"
    assert any(
        row["name"] == "structural_gate_enforced" and not row["passed"]
        for row in review["checks"]
    )


def test_compute_approval_discloses_campaign_budget() -> None:
    summary = _compute_approval_summary({"spec_path": "docs/design_spec.json"})

    assert "8 RFdiffusion3 backbone(s)" in summary
    assert "32 ProteinMPNN candidate(s)" in summary
    assert "64 AlphaFold2" in summary
    assert "Modal" in summary
