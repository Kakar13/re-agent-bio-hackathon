from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import pytest

import re_agent.immuno.e2e_pls_pickle as team_model
from re_agent.e2e_pls.netmhcpan_student import (
    MHCINetMHCpanStudent,
    save_student_cv_checkpoints,
)

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = ROOT / "models" / "chao1" / "cv5_heads.pkl 2"


class _ForbiddenPayload:
    def __reduce__(self):
        return eval, ("1 + 1",)


class _DeterministicEncoder:
    def __init__(self, *, batch_size: int = 32) -> None:
        self.batch_size = batch_size

    def embed(self, sequences: list[str]) -> list[np.ndarray]:
        outputs = []
        for sequence in sequences:
            matrix = np.zeros((len(sequence), team_model.EMBEDDING_DIM), dtype=np.float32)
            for index, residue in enumerate(sequence):
                matrix[index, (ord(residue) + index) % team_model.EMBEDDING_DIM] = 1.0
            outputs.append(matrix)
        return outputs


def test_pickle_audit_rejects_unapproved_global(tmp_path: Path) -> None:
    checkpoint = tmp_path / "unsafe.pkl"
    checkpoint.write_bytes(pickle.dumps(_ForbiddenPayload()))

    with pytest.raises(team_model.UnsafeCheckpointError):
        team_model._audit_pickle(checkpoint)


@pytest.mark.skipif(not CHECKPOINT.exists(), reason="team checkpoint was not supplied")
def test_team_checkpoint_contract_and_inference(monkeypatch: pytest.MonkeyPatch) -> None:
    bundle, digest = team_model.load_validated_bundle(CHECKPOINT)
    monkeypatch.setattr(team_model, "_ESM2SegmentEncoder", _DeterministicEncoder)

    result = team_model.TeamE2EPLSAdapter(CHECKPOINT).predict("ACDEFGHIKLMNPQRST")

    assert len(digest) == 64
    assert bundle["manifest"]["encoder_model_id"] == team_model.ENCODER_MODEL_ID
    assert result.status == "ok"
    assert result.adapter_id == "team-e2e-pls-chao1"
    assert len(result.predictions) == 9
    assert result.provenance.capability == "mhc_i_processing_surrogate"
    assert set(result.spatial_tracks) == {
        "mhci_processing_risk_max",
        "mhci_presentation_propensity_mean",
    }
    assert all(len(track) == 17 for track in result.spatial_tracks.values())
    assert any("do not interpret as MHC-II" in warning for warning in result.warnings)


@pytest.mark.skipif(not CHECKPOINT.exists(), reason="team checkpoint was not supplied")
def test_team_adapter_combines_processing_with_netmhcpan_profile(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    student_dir = tmp_path / "student"
    models = {fold: MHCINetMHCpanStudent(hidden_dim=8, dropout=0.0) for fold in range(5)}
    metrics = {
        "folds": [
            {
                "test_fold": fold,
                "validation_fold": (fold + 1) % 5,
                "training_folds": [
                    value for value in range(5) if value not in (fold, (fold + 1) % 5)
                ],
                "metrics": {"fold": fold},
            }
            for fold in range(5)
        ]
    }
    save_student_cv_checkpoints(
        student_dir,
        models,
        metrics,
        corpus_sha256="0" * 64,
    )
    monkeypatch.setattr(team_model, "_ESM2SegmentEncoder", _DeterministicEncoder)

    result = team_model.TeamE2EPLSAdapter(
        CHECKPOINT,
        netmhcpan_checkpoint_dir=student_dir,
    ).predict("ACDEFGHIKLMNPQRST")

    assert result.status == "ok"
    assert result.adapter_id.endswith("-netmhcpan")
    assert set(result.spatial_tracks) == {
        "mhci_processing_risk_max",
        "mhci_presentation_propensity_mean",
        "mhci_binding_propensity_mean",
    }
    assert result.protein_summary["risk_profile"]["overall"]["score"] >= 0.0
    assert all(row.mhc_i_binding_propensity is not None for row in result.predictions)
    assert all(
        row.mhc_i_binder_class in {"strong", "weak", "nonbinder"}
        for row in result.predictions
    )
