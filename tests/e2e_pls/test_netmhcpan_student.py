from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from re_agent.e2e_pls.netmhcpan_student import (
    MAX_IC50_NM,
    StudentTrainConfig,
    affinity_to_ic50_nm,
    format_mhci_profile,
    ic50_nm_to_affinity,
    rank_to_propensity,
    save_student_checkpoint,
    save_student_cv_checkpoints,
    screening_thresholds,
    train_student,
    train_student_cv,
)
from re_agent.e2e_pls.schema import EMBEDDING_DIM


def test_student_uses_grouped_splits_and_writes_checkpoint(tmp_path: Path) -> None:
    rng = np.random.default_rng(7)
    embeddings = rng.normal(size=(32, EMBEDDING_DIM)).astype(np.float32)
    signal = 1 / (1 + np.exp(-embeddings[:, 0]))
    frame = pd.DataFrame(
        {
            "split": ["train"] * 20 + ["val"] * 6 + ["test"] * 4 + ["challenge"] * 2,
            "netmhcpan_el_score": signal,
            "netmhcpan_ba_score": np.clip(0.8 * signal + 0.1, 0, 1),
            "netmhcpan_el_rank": 100 * (1 - signal),
            "netmhcpan_ba_rank": 100 * (1 - signal),
        }
    )

    model, metrics = train_student(
        embeddings,
        frame,
        StudentTrainConfig(hidden_dim=8, batch_size=8, epochs=2, patience=2, seed=3),
        device="cpu",
    )
    save_student_checkpoint(tmp_path, model, metrics, corpus_sha256="a" * 64)

    assert metrics["splits"]["test"]["el"]["n_rows"] == 4
    assert (tmp_path / "weights.pt").exists()
    assert (tmp_path / "manifest.json").exists()
    assert (tmp_path / "metrics.json").exists()


def test_student_cv_holds_out_each_fold_and_writes_deployment_manifest(
    tmp_path: Path,
) -> None:
    rng = np.random.default_rng(11)
    embeddings = rng.normal(size=(27, EMBEDDING_DIM)).astype(np.float32)
    signal = 1 / (1 + np.exp(-embeddings[:, 0]))
    frame = pd.DataFrame(
        {
            "split": ["train"] * 25 + ["challenge"] * 2,
            "cv_fold": list(range(5)) * 5 + [np.nan, np.nan],
            "netmhcpan_el_score": signal,
            "netmhcpan_ba_score": np.clip(0.8 * signal + 0.1, 0, 1),
            "netmhcpan_el_rank": 100 * (1 - signal),
            "netmhcpan_ba_rank": 100 * (1 - signal),
        }
    )

    models, metrics = train_student_cv(
        embeddings,
        frame,
        StudentTrainConfig(hidden_dim=8, batch_size=8, epochs=1, patience=1, seed=3),
        device="cpu",
    )
    save_student_cv_checkpoints(
        tmp_path,
        models,
        metrics,
        corpus_sha256="b" * 64,
    )

    assert set(models) == set(range(5))
    assert len(metrics["folds"]) == 5
    for fold_report in metrics["folds"]:
        assert fold_report["validation_fold"] == (fold_report["test_fold"] + 1) % 5
        assert len(fold_report["training_folds"]) == 3
        assert fold_report["metrics"]["splits"]["test"]["el"]["n_rows"] == 5
        assert fold_report["metrics"]["splits"]["challenge"]["el"]["n_rows"] == 2
    assert "mean" in metrics["aggregate"]["held_out_test"]["el"]["propensity_mae"]
    assert "std" in metrics["aggregate"]["held_out_test"]["ba"]["propensity_rmse"]
    assert metrics["aggregate"]["pooled_out_of_fold"]["el"]["n_rows"] == 25

    for fold in range(5):
        assert (tmp_path / f"fold_{fold}" / "weights.pt").exists()
    deployment = json.loads((tmp_path / "deployment_manifest.json").read_text())
    assert deployment["artifact_type"] == "five_fold_cv_ensemble"
    assert len(deployment["checkpoints"]) == 5
    assert deployment["not_ground_truth"] is True
    assert (tmp_path / "metrics.json").exists()


@pytest.mark.parametrize(
    ("rank", "binder_class", "risk_band"),
    [
        (0.5, "strong", "high"),
        (2.0, "weak", "moderate"),
        (2.01, "nonbinder", "low"),
    ],
)
def test_mhci_profile_formats_propensities_without_double_counting(
    rank: float,
    binder_class: str,
    risk_band: str,
) -> None:
    el_propensity = float(rank_to_propensity(np.asarray(rank)))
    profile = format_mhci_profile(el_propensity, 0.1)

    assert profile["predicted_el_rank"] == pytest.approx(rank, abs=1e-4)
    assert profile["binder_class"] == binder_class
    assert profile["binder_class_basis"] == "predicted_el_rank"
    assert profile["risk_band"] == risk_band
    assert profile["overall_mhci_risk"] == profile["el_presentation_propensity"]
    assert profile["overall_mhci_risk"] != pytest.approx(
        (profile["el_presentation_propensity"] + profile["ba_binding_propensity"]) / 2
    )
    assert profile["caveats"]
    assert "predicted_ba_ic50_nm" not in profile


def test_affinity_transform_round_trips_to_nanomolar() -> None:
    ic50 = np.array([1.0, 50.0, 500.0, 5_000.0, MAX_IC50_NM])

    recovered = affinity_to_ic50_nm(ic50_nm_to_affinity(ic50))

    assert recovered == pytest.approx(ic50, rel=1e-9)
    # A strong binder must map to a higher score than a weak one.
    assert ic50_nm_to_affinity(np.asarray(50.0)) > ic50_nm_to_affinity(np.asarray(5_000.0))


def test_profile_reports_affinity_when_the_channel_is_present() -> None:
    score = float(ic50_nm_to_affinity(np.asarray(120.0)))

    profile = format_mhci_profile(0.9, 0.8, score)

    assert profile["predicted_ba_ic50_nm"] == pytest.approx(120.0, rel=1e-6)
    assert profile["ba_affinity_score"] == pytest.approx(score)


def test_screening_thresholds_trade_precision_for_recall() -> None:
    # A student that ranks binders correctly but compresses them toward the
    # middle: the default rank cutoff misses some, a lower cutoff recovers them.
    teacher_ranks = np.concatenate([np.full(20, 0.1), np.linspace(3.0, 90.0, 180)])
    propensities = np.concatenate(
        [np.linspace(0.95, 0.60, 20), np.linspace(0.59, 0.01, 180)]
    )

    report = screening_thresholds(
        propensities, teacher_ranks, binder_rank=0.5, target_recalls=(0.5, 1.0)
    )

    assert report["n_teacher_binders"] == 20
    assert report["default_cutoff"]["recall"] < 1.0
    points = {point["target_recall"]: point for point in report["operating_points"]}
    assert points[1.0]["recall"] == pytest.approx(1.0)
    assert points[1.0]["n_flagged"] >= points[0.5]["n_flagged"]
    assert points[1.0]["precision"] <= points[0.5]["precision"]


def test_screen_flags_peptides_the_binder_class_calls_nonbinders() -> None:
    # A peptide just past the weak-binder rank: NetMHCpan convention says
    # nonbinder, but the recall-tuned screen should still surface it.
    borderline = float(rank_to_propensity(np.asarray(1.4)))

    profile = format_mhci_profile(borderline, 0.5)

    assert profile["binder_class"] == "weak"
    assert profile["screening_flag"] is True

    clear_nonbinder = format_mhci_profile(float(rank_to_propensity(np.asarray(40.0))), 0.1)
    assert clear_nonbinder["binder_class"] == "nonbinder"
    assert clear_nonbinder["screening_flag"] is False


def test_screening_threshold_is_configurable() -> None:
    propensity = 0.80

    assert format_mhci_profile(propensity, 0.5, screening_threshold=0.70)["screening_flag"]
    assert not format_mhci_profile(
        propensity, 0.5, screening_threshold=0.95
    )["screening_flag"]


def test_screening_thresholds_require_both_classes() -> None:
    with pytest.raises(ValueError, match="both binders and non-binders"):
        screening_thresholds(
            np.array([0.9, 0.8, 0.7]),
            np.array([10.0, 20.0, 30.0]),
            binder_rank=0.5,
        )
