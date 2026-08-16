from __future__ import annotations

import pytest

from re_agent.agent import run_store
from re_agent.agent.tools import replay_latest_design_campaign


def test_replay_uses_only_hash_valid_local_campaign(monkeypatch, tmp_path) -> None:
    runs_dir = tmp_path / "results" / "workbench" / "runs"
    structure = tmp_path / "results" / "design_campaigns" / "design-1" / "candidate.pdb"
    structure.parent.mkdir(parents=True)
    structure.write_text("ATOM\n")
    monkeypatch.setattr(run_store, "ROOT", tmp_path)
    monkeypatch.setattr(run_store, "RUNS_DIR", runs_dir)

    source = run_store.write_run_artifact(
        run_id="design-live",
        kind="design_to_screen_campaign",
        payload={
            "title": "Live campaign",
            "manifest": {
                "status": "completed",
                "candidates": [
                    {
                        "refolded_structure_path": (
                            "results/design_campaigns/design-1/candidate.pdb"
                        )
                    }
                ],
            },
            "citations": [{"url": "https://paperclip.gxl.ai/citations/papers/test#L1"}],
        },
    )

    replay = replay_latest_design_campaign.invoke({})

    assert replay["kind"] == "design_to_screen_campaign_replay"
    assert replay["payload"]["execution_mode"] == "replay"
    assert replay["payload"]["replay_of"] == source["id"]
    assert replay["payload"]["replay_source_sha256"] == source["sha256"]


def test_latest_campaign_rejects_tampered_artifact(monkeypatch, tmp_path) -> None:
    runs_dir = tmp_path / "results" / "workbench" / "runs"
    monkeypatch.setattr(run_store, "ROOT", tmp_path)
    monkeypatch.setattr(run_store, "RUNS_DIR", runs_dir)
    source = run_store.write_run_artifact(
        run_id="design-live",
        kind="design_to_screen_campaign",
        payload={"title": "Live campaign", "manifest": {"status": "completed"}},
    )
    (tmp_path / source["path"]).write_text('{"tampered": true}\n')

    with pytest.raises(FileNotFoundError, match="no valid design_to_screen_campaign"):
        run_store.latest_run_artifact("design_to_screen_campaign")
