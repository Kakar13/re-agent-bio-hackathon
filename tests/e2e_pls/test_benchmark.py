from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from re_agent.e2e_pls.benchmark import (
    assign_pda_novelty_splits,
    binder_agreement,
    evaluate_external_scores,
    imitation_metrics,
    teacher_agreement,
)


def test_novelty_split_uses_highest_parent_novelty_per_component() -> None:
    frame = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD", "EEEEEEEEE"],
            "parent_component_id": ["shared", "train", "val", "excluded"],
            "occurrence_parent_ids": [
                json.dumps(["p-natural", "p-novel"]),
                json.dumps(["p-similar"]),
                json.dumps(["p-distant"]),
                json.dumps(["p-unscored"]),
            ],
        }
    )
    novelty = {
        "p-natural": "natural_like",
        "p-novel": "novel",
        "p-similar": "similar",
        "p-distant": "distant",
        "p-unscored": "unscored",
    }

    split, manifest = assign_pda_novelty_splits(frame, novelty)

    assert split.set_index("peptide")["split"].to_dict() == {
        "AAAAAAAAA": "test",
        "CCCCCCCCC": "train",
        "DDDDDDDDD": "val",
        "EEEEEEEEE": "exclude",
    }
    assert manifest["row_counts"] == {"test": 1, "train": 1, "val": 1, "exclude": 1}


def test_external_metrics_reward_correct_score_direction() -> None:
    frame = pd.DataFrame(
        {
            "resp": [False, False, True, True],
            "r_mean": [0.0, 0.1, 0.9, 1.0],
            "risk": [0.1, 0.2, 0.8, 0.9],
            "reversed": [0.9, 0.8, 0.2, 0.1],
        }
    )

    report = evaluate_external_scores(
        frame,
        {"risk": "risk", "reversed": "reversed"},
        bootstrap_samples=20,
        seed=4,
    )

    assert report["metrics"]["risk"]["roc_auc"] == 1.0
    assert report["metrics"]["risk"]["average_precision"] == 1.0
    assert report["metrics"]["reversed"]["roc_auc"] == 0.0
    assert report["metrics"]["risk"]["bootstrap_95_ci"]["roc_auc"] == [1.0, 1.0]


def test_imitation_metrics_are_zero_error_on_an_exact_copy() -> None:
    teacher = np.array([0.05, 0.2, 0.55, 0.9])

    report = imitation_metrics(teacher.copy(), teacher)

    assert report["spearman_vs_teacher_propensity"] == pytest.approx(1.0)
    assert report["propensity_mae"] == pytest.approx(0.0)
    assert report["propensity_rmse"] == pytest.approx(0.0)


def test_imitation_metrics_reject_misaligned_inputs() -> None:
    with pytest.raises(ValueError, match="align row for row"):
        imitation_metrics(np.array([0.1, 0.2]), np.array([0.1]))


def test_binder_agreement_counts_each_side_separately() -> None:
    # The student misses one teacher binder and invents none, the conservative
    # direction of error seen in the real cohorts.
    report = binder_agreement(
        np.array([0.1, 0.4, 3.0, 8.0]),
        np.array([0.1, 0.4, 0.3, 8.0]),
        0.5,
    )

    assert report["teacher_binders"] == 3
    assert report["student_binders"] == 2
    assert report["agreement"] == pytest.approx(0.75)


def test_teacher_agreement_separates_auc_from_average_precision() -> None:
    # One strong binder among nine peptides: a score that ranks it last still
    # earns a respectable AUC, so average precision has to carry the signal.
    teacher_rank = np.array([0.2, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    teacher_propensity = 1.0 - np.log1p(teacher_rank) / np.log(101.0)
    faithful = teacher_propensity.copy()
    blunt = np.array([0.4, 0.9, 0.3, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1])

    good = teacher_agreement(
        faithful, teacher_rank, teacher_propensity, strong_rank=0.5, weak_rank=2.0
    )
    poor = teacher_agreement(
        blunt, teacher_rank, teacher_propensity, strong_rank=0.5, weak_rank=2.0
    )

    assert good["strong_binder_average_precision"] == pytest.approx(1.0)
    assert good["strong_binder_roc_auc"] == pytest.approx(1.0)
    assert poor["strong_binder_average_precision"] < good["strong_binder_average_precision"]
    assert poor["strong_binder_roc_auc"] > 0.5


def test_teacher_agreement_requires_both_classes() -> None:
    ranks = np.array([5.0, 6.0, 7.0])
    with pytest.raises(ValueError, match="both strong binders and non-binders"):
        teacher_agreement(
            np.array([0.1, 0.2, 0.3]),
            ranks,
            1.0 - np.log1p(ranks) / np.log(101.0),
            strong_rank=0.5,
            weak_rank=2.0,
        )
