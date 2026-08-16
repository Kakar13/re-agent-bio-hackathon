"""Leakage-safe PDA novelty splits and external binary benchmark metrics."""

from __future__ import annotations

import json
from collections.abc import Mapping

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

NOVELTY_PRIORITY = {
    "unscored": 0,
    "natural_like": 1,
    "similar": 2,
    "distant": 3,
    "novel": 4,
}
NOVELTY_SPLIT = {
    0: "exclude",
    1: "train",
    2: "train",
    3: "val",
    4: "test",
}


def assign_pda_novelty_splits(
    frame: pd.DataFrame,
    parent_novelty: Mapping[str, str],
) -> tuple[pd.DataFrame, dict]:
    """Assign whole shared-peptide components to train/val/novel-test cohorts."""

    required = {"parent_component_id", "occurrence_parent_ids"}
    missing = required - set(frame)
    if missing:
        raise ValueError(f"PDA frame is missing novelty split columns: {sorted(missing)}")

    component_parents: dict[str, set[str]] = {}
    for component, encoded_parents in zip(
        frame["parent_component_id"].astype(str),
        frame["occurrence_parent_ids"].astype(str),
        strict=True,
    ):
        try:
            parents = json.loads(encoded_parents)
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid occurrence_parent_ids JSON: {encoded_parents!r}") from exc
        if not isinstance(parents, list) or not all(isinstance(value, str) for value in parents):
            raise ValueError("occurrence_parent_ids must encode a list of strings")
        component_parents.setdefault(component, set()).update(parents)

    missing_parents = sorted(
        {
            parent
            for parents in component_parents.values()
            for parent in parents
            if parent not in parent_novelty
        }
    )
    if missing_parents:
        raise ValueError(f"PDA novelty metadata is missing parents: {missing_parents[:5]}")

    component_level = {
        component: max(NOVELTY_PRIORITY.get(str(parent_novelty[parent]), 0) for parent in parents)
        for component, parents in component_parents.items()
    }
    out = frame.copy()
    out["component_novelty_level"] = (
        out["parent_component_id"].astype(str).map(component_level).astype("int8")
    )
    out["split"] = out["component_novelty_level"].map(NOVELTY_SPLIT)
    if out["split"].isna().any():
        raise RuntimeError("failed to assign every PDA row to a novelty cohort")

    component_summary = pd.Series(component_level).map(NOVELTY_SPLIT).value_counts()
    row_summary = out["split"].value_counts()
    manifest = {
        "definition": {
            "train": ["natural_like", "similar"],
            "val": ["distant"],
            "test": ["novel"],
            "exclude": ["unscored"],
            "component_rule": "highest novelty among all parents in a shared-9-mer component",
        },
        "row_counts": {str(key): int(value) for key, value in row_summary.items()},
        "component_counts": {str(key): int(value) for key, value in component_summary.items()},
    }
    return out, manifest


def imitation_metrics(student: np.ndarray, teacher: np.ndarray) -> dict:
    """Score one student channel against its teacher on matched rows."""

    student = np.asarray(student, dtype=np.float64)
    teacher = np.asarray(teacher, dtype=np.float64)
    if student.shape != teacher.shape:
        raise ValueError("student and teacher propensities must align row for row")
    if not (np.isfinite(student).all() and np.isfinite(teacher).all()):
        raise ValueError("imitation inputs must be finite")
    return {
        "spearman_vs_teacher_propensity": float(
            pd.Series(student).corr(pd.Series(teacher), method="spearman")
        ),
        "pearson_vs_teacher_propensity": float(
            pd.Series(student).corr(pd.Series(teacher), method="pearson")
        ),
        "propensity_mae": float(np.mean(np.abs(student - teacher))),
        "propensity_rmse": float(np.sqrt(np.mean((student - teacher) ** 2))),
    }


def binder_agreement(
    predicted_rank: np.ndarray,
    teacher_rank: np.ndarray,
    cutoff: float,
) -> dict:
    """Compare student and teacher binder calls at one percentile-rank cutoff."""

    predicted = np.asarray(predicted_rank, dtype=np.float64) <= cutoff
    actual = np.asarray(teacher_rank, dtype=np.float64) <= cutoff
    if predicted.shape != actual.shape:
        raise ValueError("predicted and teacher ranks must align row for row")
    return {
        "rank_cutoff": float(cutoff),
        "teacher_binders": int(actual.sum()),
        "student_binders": int(predicted.sum()),
        "agreement": float(np.mean(predicted == actual)),
    }


def teacher_agreement(
    values: np.ndarray,
    teacher_rank: np.ndarray,
    teacher_propensity: np.ndarray,
    *,
    strong_rank: float,
    weak_rank: float,
) -> dict:
    """Rank a single score against NetMHCpan ranks on identical rows.

    Average precision is reported alongside ROC AUC because strong binders are a
    small minority; AUC alone flatters a score that cannot surface them.
    """

    values = np.asarray(values, dtype=np.float64)
    teacher_rank = np.asarray(teacher_rank, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("scores must be finite")
    strong = teacher_rank <= strong_rank
    weak = teacher_rank <= weak_rank
    if not strong.any() or strong.all():
        raise ValueError("teacher ranks must contain both strong binders and non-binders")
    return {
        "spearman_vs_teacher_el_rank": float(
            pd.Series(values).corr(pd.Series(teacher_propensity), method="spearman")
        ),
        "strong_binder_roc_auc": float(roc_auc_score(strong, values)),
        "strong_binder_average_precision": float(average_precision_score(strong, values)),
        "weak_binder_roc_auc": float(roc_auc_score(weak, values)),
        "weak_binder_average_precision": float(average_precision_score(weak, values)),
    }


def _score_metrics(labels: np.ndarray, values: np.ndarray, continuous: np.ndarray) -> dict:
    return {
        "roc_auc": float(roc_auc_score(labels, values)),
        "average_precision": float(average_precision_score(labels, values)),
        "spearman_vs_activity": float(
            pd.Series(values).corr(pd.Series(continuous), method="spearman")
        ),
    }


def evaluate_external_scores(
    frame: pd.DataFrame,
    score_columns: Mapping[str, str],
    *,
    label_column: str = "resp",
    activity_column: str = "r_mean",
    bootstrap_samples: int = 1_000,
    seed: int = 0,
) -> dict:
    """Evaluate model scores against experimental activation with bootstrap CIs."""

    required = {label_column, activity_column, *score_columns.values()}
    missing = required - set(frame)
    if missing:
        raise ValueError(f"external benchmark is missing columns: {sorted(missing)}")

    labels = frame[label_column].astype(bool).to_numpy()
    activity = frame[activity_column].to_numpy(dtype=np.float64)
    if len(labels) < 2 or labels.all() or (~labels).all():
        raise ValueError("external benchmark requires both positive and negative outcomes")
    if not np.isfinite(activity).all():
        raise ValueError("external activity values must be finite")

    rng = np.random.default_rng(seed)
    report = {}
    for name, column in score_columns.items():
        values = frame[column].to_numpy(dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError(f"benchmark score {name!r} contains non-finite values")
        metrics = _score_metrics(labels, values, activity)
        bootstraps = {metric: [] for metric in metrics}
        completed = 0
        attempts = 0
        while completed < bootstrap_samples and attempts < bootstrap_samples * 10:
            attempts += 1
            indices = rng.integers(0, len(frame), size=len(frame))
            sampled_labels = labels[indices]
            if sampled_labels.all() or (~sampled_labels).all():
                continue
            sampled = _score_metrics(
                sampled_labels,
                values[indices],
                activity[indices],
            )
            for metric, value in sampled.items():
                if np.isfinite(value):
                    bootstraps[metric].append(value)
            completed += 1
        metrics["bootstrap_95_ci"] = {
            metric: [
                float(np.quantile(samples, 0.025)),
                float(np.quantile(samples, 0.975)),
            ]
            for metric, samples in bootstraps.items()
            if samples
        }
        report[name] = metrics

    return {
        "n_rows": len(frame),
        "n_positive": int(labels.sum()),
        "n_negative": int((~labels).sum()),
        "metrics": report,
        "interpretation": (
            "Experimental T-cell activation is downstream of MHC presentation and "
            "TCR recognition; it is not a direct NetMHCpan teacher target."
        ),
    }
