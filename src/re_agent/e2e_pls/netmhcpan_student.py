"""Frozen-ESM student head for separate NetMHCpan EL and BA targets.

Local HLA-A*02:01 teacher surrogate for campaign downselection. Measured MHC
labels on de novo proteins are almost unavailable, so this head distills
NetMHCpan onto ESM embeddings of Protein Design Archive 9-mers. Primary
intended use is cytosolic MHC-I presentation when a designed protein is
expressed inside a nucleated cell from a genetic payload. It does not
predict TCR recognition, danger, MHC-II/ADA, or clinical immunogenicity.
"""

from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score
from torch import nn

from re_agent.e2e_pls.schema import EMBEDDING_DIM, ENCODER_MODEL_ID

TARGET_COLUMNS = ("netmhcpan_el_score", "netmhcpan_ba_score")
RANK_COLUMNS = ("netmhcpan_el_rank", "netmhcpan_ba_rank")
AFFINITY_COLUMN = "netmhcpan_ba_score"
MODEL_VERSION = "mhci-netmhcpan-student-v3"
# v2 heads predict the two rank channels only. They still load so that a
# deployment can keep serving while a v3 ensemble is being retrained.
SUPPORTED_MODEL_VERSIONS = ("mhci-netmhcpan-student-v2", "mhci-netmhcpan-student-v3")
RANK_OUTPUTS = 2
AFFINITY_INDEX = 2
DEFAULT_OUTPUTS = 3
# NetMHCpan reports BA on the field-standard 1 - log(IC50)/log(50000) scale,
# which is exactly invertible back to nanomolar.
MAX_IC50_NM = 50_000.0
CV_FOLDS = tuple(range(5))
# Screening cutoff on EL propensity, chosen for 95% recall of the teacher's
# strong binders over 159,038 out-of-fold PDA peptides. The binder_class rule
# below keeps NetMHCpan's conventional rank thresholds, which recover only 52%
# of them; a screen needs the opposite error direction, so it gets its own
# cutoff instead of redefining what "strong binder" means.
# Source: results/benchmarks/screening_calibration/metrics.json
SCREENING_PROPENSITY_THRESHOLD = 0.7892
SCREENING_RECALL_TARGET = 0.95
SCREENING_PRECISION_AT_TARGET = 0.439
STRONG_BINDER_RANK = 0.5
WEAK_BINDER_RANK = 2.0
PROFILE_CAVEATS = (
    "Student outputs estimate agreement with the NetMHCpan teacher, not biological ground truth.",
    "EL presentation and BA binding are related teacher channels; overall MHC-I risk uses EL "
    "presentation propensity alone to avoid double-counting BA.",
    "Predicted percentile ranks are deterministic inversions of student propensities, not direct "
    "NetMHCpan calculations.",
)


def rank_to_propensity(ranks: np.ndarray) -> np.ndarray:
    """Map NetMHCpan percentile ranks to a smooth higher-is-riskier target."""

    clipped = np.clip(np.asarray(ranks, dtype=np.float32), 0.0, 100.0)
    return 1.0 - np.log1p(clipped) / np.log(101.0)


def propensity_to_rank(propensities: np.ndarray) -> np.ndarray:
    """Invert :func:`rank_to_propensity` into predicted percentile ranks."""

    clipped = np.clip(np.asarray(propensities, dtype=np.float64), 0.0, 1.0)
    return np.expm1((1.0 - clipped) * np.log(101.0))


def affinity_to_ic50_nm(scores: np.ndarray) -> np.ndarray:
    """Convert NetMHCpan BA scores back to nanomolar IC50."""

    clipped = np.clip(np.asarray(scores, dtype=np.float64), 0.0, 1.0)
    return np.power(MAX_IC50_NM, 1.0 - clipped)


def ic50_nm_to_affinity(ic50_nm: np.ndarray) -> np.ndarray:
    """Convert nanomolar IC50 to the bounded NetMHCpan BA score."""

    clipped = np.clip(np.asarray(ic50_nm, dtype=np.float64), 1.0, MAX_IC50_NM)
    return 1.0 - np.log(clipped) / np.log(MAX_IC50_NM)


def format_mhci_profile(
    el_presentation_propensity: float,
    ba_binding_propensity: float,
    ba_affinity_score: float | None = None,
    screening_threshold: float = SCREENING_PROPENSITY_THRESHOLD,
) -> dict[str, Any]:
    """Format student outputs as a deterministic, agent-facing MHC-I profile.

    ``ba_affinity_score`` is optional so that a v2 rank-only ensemble still
    formats; when present it adds the nanomolar affinity an immunologist reads.

    Two different calls come out of this. ``binder_class`` answers "what would
    NetMHCpan call this", on NetMHCpan's own rank thresholds. ``screening_flag``
    answers "should a designer look at this", at a cutoff tuned for recall.
    They disagree on purpose: the second deliberately over-calls.
    """

    propensities = np.asarray(
        [el_presentation_propensity, ba_binding_propensity], dtype=np.float64
    )
    if not np.isfinite(propensities).all():
        raise ValueError("MHC-I propensities must be finite")
    if ((propensities < 0.0) | (propensities > 1.0)).any():
        raise ValueError("MHC-I propensities must be between 0 and 1")

    el_rank, ba_rank = (float(value) for value in propensity_to_rank(propensities))
    if el_rank <= STRONG_BINDER_RANK or np.isclose(el_rank, STRONG_BINDER_RANK):
        binder_class = "strong"
        risk_band = "high"
    elif el_rank <= WEAK_BINDER_RANK or np.isclose(el_rank, WEAK_BINDER_RANK):
        binder_class = "weak"
        risk_band = "moderate"
    else:
        binder_class = "nonbinder"
        risk_band = "low"

    profile = {
        "el_presentation_propensity": float(propensities[0]),
        "ba_binding_propensity": float(propensities[1]),
        "predicted_el_rank": el_rank,
        "predicted_ba_rank": ba_rank,
        "binder_class": binder_class,
        "binder_class_basis": "predicted_el_rank",
        "screening_flag": bool(propensities[0] >= screening_threshold),
        "screening_basis": (
            f"EL propensity >= {screening_threshold:.4f}, calibrated for "
            f"{SCREENING_RECALL_TARGET:.0%} recall of teacher strong binders "
            f"at {SCREENING_PRECISION_AT_TARGET:.0%} precision"
        ),
        "overall_mhci_risk": float(propensities[0]),
        "risk_band": risk_band,
        "caveats": list(PROFILE_CAVEATS),
    }
    if ba_affinity_score is not None:
        score = float(ba_affinity_score)
        if not np.isfinite(score) or not 0.0 <= score <= 1.0:
            raise ValueError("BA affinity score must be finite and between 0 and 1")
        profile["ba_affinity_score"] = score
        profile["predicted_ba_ic50_nm"] = float(affinity_to_ic50_nm(np.asarray(score)))
    return profile


class MHCINetMHCpanStudent(nn.Module):
    """Small multi-task head over a mean-pooled frozen ESM-2 9-mer.

    Channel 0 is EL rank propensity, channel 1 is BA rank propensity, and
    channel 2 (when present) is the BA affinity score that inverts to IC50.
    """

    def __init__(
        self,
        hidden_dim: int = 256,
        dropout: float = 0.1,
        n_outputs: int = DEFAULT_OUTPUTS,
    ) -> None:
        super().__init__()
        if n_outputs not in (RANK_OUTPUTS, DEFAULT_OUTPUTS):
            raise ValueError(f"student supports {RANK_OUTPUTS} or {DEFAULT_OUTPUTS} outputs")
        self.hidden_dim = hidden_dim
        self.dropout = dropout
        self.n_outputs = n_outputs
        self.net = nn.Sequential(
            nn.LayerNorm(EMBEDDING_DIM),
            nn.Linear(EMBEDDING_DIM, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_outputs),
        )

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        return torch.sigmoid(self.net(embeddings))


def load_student_ensemble(
    checkpoint_dir: Path,
) -> tuple[list[MHCINetMHCpanStudent], str, dict[str, Any]]:
    """Load a validated five-fold deployment ensemble without unsafe pickle globals."""

    checkpoint_dir = Path(checkpoint_dir).resolve()
    manifest_path = checkpoint_dir / "deployment_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("artifact_type") != "five_fold_cv_ensemble":
        raise ValueError(f"not a five-fold student ensemble: {manifest_path}")
    if manifest.get("model_version") not in SUPPORTED_MODEL_VERSIONS:
        raise ValueError(
            f"unsupported student version {manifest.get('model_version')!r}; "
            f"expected one of {list(SUPPORTED_MODEL_VERSIONS)}"
        )
    if manifest.get("encoder_model_id") != ENCODER_MODEL_ID:
        raise ValueError("student ensemble encoder does not match the runtime encoder")

    models: list[MHCINetMHCpanStudent] = []
    digest = hashlib.sha256(manifest_path.read_bytes())
    checkpoints = manifest.get("checkpoints", [])
    if len(checkpoints) != len(CV_FOLDS):
        raise ValueError("student deployment manifest must reference five fold checkpoints")
    for expected_fold, checkpoint in zip(CV_FOLDS, checkpoints, strict=True):
        if checkpoint.get("test_fold") != expected_fold:
            raise ValueError("student fold checkpoints must be ordered 0..4")
        weights_path = (checkpoint_dir / checkpoint["checkpoint"]).resolve()
        fold_manifest_path = (checkpoint_dir / checkpoint["manifest"]).resolve()
        if (
            checkpoint_dir not in weights_path.parents
            or checkpoint_dir not in fold_manifest_path.parents
        ):
            raise ValueError("student checkpoint path escapes its deployment directory")
        fold_manifest = json.loads(fold_manifest_path.read_text())
        model = MHCINetMHCpanStudent(
            hidden_dim=int(fold_manifest["hidden_dim"]),
            dropout=float(fold_manifest["dropout"]),
            n_outputs=int(fold_manifest.get("n_outputs", RANK_OUTPUTS)),
        )
        state = torch.load(weights_path, map_location="cpu", weights_only=True)
        model.load_state_dict(state, strict=True)
        model.eval()
        models.append(model)
        digest.update(fold_manifest_path.read_bytes())
        digest.update(weights_path.read_bytes())
    return models, digest.hexdigest(), manifest


def predict_student_ensemble(
    models: list[MHCINetMHCpanStudent],
    embeddings: np.ndarray,
    *,
    batch_size: int = 512,
) -> np.ndarray:
    """Average the five fold heads into EL and BA rank-propensity predictions."""

    if len(models) != len(CV_FOLDS):
        raise ValueError("student ensemble requires exactly five fold heads")
    widths = {model.n_outputs for model in models}
    if len(widths) != 1:
        raise ValueError(f"fold heads disagree on output width: {sorted(widths)}")
    values = np.asarray(embeddings, dtype=np.float32)
    if values.ndim != 2 or values.shape[1] != EMBEDDING_DIM:
        raise ValueError(f"expected embeddings with shape (n, {EMBEDDING_DIM})")
    predictions = np.empty((len(values), widths.pop()), dtype=np.float32)
    tensor = torch.from_numpy(values)
    with torch.inference_mode():
        for offset in range(0, len(values), batch_size):
            batch = tensor[offset : offset + batch_size]
            fold_outputs = torch.stack([model(batch) for model in models])
            predictions[offset : offset + len(batch)] = fold_outputs.mean(dim=0).numpy()
    return predictions


@dataclass(frozen=True)
class StudentTrainConfig:
    hidden_dim: int = 256
    dropout: float = 0.1
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    batch_size: int = 512
    epochs: int = 60
    patience: int = 8
    binder_rank_threshold: float = 2.0
    binder_weight: float = 4.0
    n_outputs: int = DEFAULT_OUTPUTS
    seed: int = 0


def _metrics(
    frame: pd.DataFrame,
    predictions: np.ndarray,
    mask: np.ndarray,
    binder_rank_threshold: float,
) -> dict:
    if not mask.any():
        return {"warning": "split has no rows"}
    metrics: dict[str, dict] = {}
    for index, channel in enumerate(("el", "ba")):
        rank = frame[RANK_COLUMNS[index]].to_numpy(dtype=np.float32)[mask]
        raw_score = frame[TARGET_COLUMNS[index]].to_numpy(dtype=np.float32)[mask]
        target = rank_to_propensity(rank)
        predicted = predictions[mask, index]
        binders = rank <= binder_rank_threshold
        channel_metrics = {
            "n_rows": int(mask.sum()),
            "spearman_vs_rank": float(
                pd.Series(predicted).corr(pd.Series(-rank), method="spearman")
            ),
            "spearman_vs_raw_score": float(
                pd.Series(predicted).corr(pd.Series(raw_score), method="spearman")
            ),
            "propensity_mae": float(np.mean(np.abs(predicted - target))),
            "propensity_rmse": float(np.sqrt(np.mean((predicted - target) ** 2))),
            "binder_prevalence": float(binders.mean()),
        }
        if binders.any() and (~binders).any():
            channel_metrics["binder_auprc"] = float(
                average_precision_score(binders.astype(np.int8), predicted)
            )
        metrics[channel] = channel_metrics

    if predictions.shape[1] > AFFINITY_INDEX:
        actual_score = frame[AFFINITY_COLUMN].to_numpy(dtype=np.float64)[mask]
        predicted_score = predictions[mask, AFFINITY_INDEX].astype(np.float64)
        actual_ic50 = affinity_to_ic50_nm(actual_score)
        predicted_ic50 = affinity_to_ic50_nm(predicted_score)
        # Affinity error is meaningful multiplicatively, not additively: being
        # 100 nM off matters at 50 nM and is noise at 30,000 nM.
        fold_error = np.exp(np.abs(np.log(predicted_ic50) - np.log(actual_ic50)))
        binders = actual_ic50 <= 500.0
        affinity_metrics = {
            "n_rows": int(mask.sum()),
            "score_mae": float(np.mean(np.abs(predicted_score - actual_score))),
            "spearman_vs_ic50": float(
                pd.Series(predicted_score).corr(pd.Series(-actual_ic50), method="spearman")
            ),
            "median_ic50_fold_error": float(np.median(fold_error)),
            "binder_prevalence_ic50_500nm": float(binders.mean()),
        }
        if binders.any():
            affinity_metrics["median_ic50_fold_error_binders"] = float(
                np.median(fold_error[binders])
            )
        if binders.any() and (~binders).any():
            affinity_metrics["binder_auprc_ic50_500nm"] = float(
                average_precision_score(binders.astype(np.int8), predicted_score)
            )
        metrics["affinity"] = affinity_metrics
    return metrics


def train_student(
    embeddings: np.ndarray,
    frame: pd.DataFrame,
    config: StudentTrainConfig | None = None,
    device: str | torch.device | None = None,
) -> tuple[MHCINetMHCpanStudent, dict]:
    """Train on ``train``, select on ``val``, and report untouched test/challenge."""

    config = config or StudentTrainConfig()
    if embeddings.shape != (len(frame), EMBEDDING_DIM):
        raise ValueError(
            f"expected embeddings shape {(len(frame), EMBEDDING_DIM)}, got {embeddings.shape}"
        )
    missing = set(TARGET_COLUMNS + RANK_COLUMNS + ("split",)) - set(frame)
    if missing:
        raise ValueError(f"training frame is missing columns: {sorted(missing)}")

    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    ranks = np.array(frame[list(RANK_COLUMNS)], dtype=np.float32, copy=True)
    raw_scores = np.array(frame[list(TARGET_COLUMNS)], dtype=np.float32, copy=True)
    target = rank_to_propensity(ranks)
    if not np.isfinite(raw_scores).all() or not np.isfinite(ranks).all():
        raise ValueError("teacher targets and ranks must be finite")

    if config.n_outputs > RANK_OUTPUTS:
        affinity = np.array(frame[AFFINITY_COLUMN], dtype=np.float32, copy=True)
        if not np.isfinite(affinity).all():
            raise ValueError("teacher affinity scores must be finite")
        if ((affinity < 0.0) | (affinity > 1.0)).any():
            raise ValueError(f"{AFFINITY_COLUMN} must lie on the bounded 0-1 BA scale")
        target = np.concatenate([target, affinity[:, None]], axis=1)

    split = frame["split"].astype(str).to_numpy()
    train_indices = np.flatnonzero(split == "train")
    val_indices = np.flatnonzero(split == "val")
    if not len(train_indices) or not len(val_indices):
        raise ValueError("grouped corpus must contain non-empty train and val splits")

    resolved_device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
    model = MHCINetMHCpanStudent(
        config.hidden_dim,
        config.dropout,
        n_outputs=config.n_outputs,
    ).to(resolved_device)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.learning_rate,
        weight_decay=config.weight_decay,
    )
    x = torch.from_numpy(embeddings.astype(np.float32, copy=False))
    y = torch.from_numpy(target)
    sample_weight = 1.0 + config.binder_weight * (
        ranks <= config.binder_rank_threshold
    ).mean(axis=1)
    weights = torch.from_numpy(sample_weight.astype(np.float32))

    rng = np.random.default_rng(config.seed)
    best_state = copy.deepcopy(model.state_dict())
    best_val_loss = float("inf")
    best_epoch = -1
    history = []
    stale_epochs = 0

    for epoch in range(config.epochs):
        model.train()
        order = rng.permutation(train_indices)
        train_loss_sum = 0.0
        for offset in range(0, len(order), config.batch_size):
            indices = order[offset : offset + config.batch_size]
            xb = x[indices].to(resolved_device)
            yb = y[indices].to(resolved_device)
            wb = weights[indices].to(resolved_device)
            optimizer.zero_grad()
            prediction = model(xb)
            per_row = torch.nn.functional.smooth_l1_loss(
                prediction, yb, reduction="none"
            ).mean(dim=1)
            loss = (per_row * wb).sum() / wb.sum()
            loss.backward()
            optimizer.step()
            train_loss_sum += float(loss.item()) * len(indices)

        model.eval()
        with torch.no_grad():
            val_prediction = model(x[val_indices].to(resolved_device))
            val_loss = float(
                torch.nn.functional.smooth_l1_loss(
                    val_prediction,
                    y[val_indices].to(resolved_device),
                ).item()
            )
        history.append(
            {
                "epoch": epoch,
                "train_loss": train_loss_sum / max(1, len(order)),
                "val_loss": val_loss,
            }
        )
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            stale_epochs = 0
        else:
            stale_epochs += 1
            if stale_epochs >= config.patience:
                break

    model.load_state_dict(best_state)
    model.eval()
    predictions = np.empty((len(frame), config.n_outputs), dtype=np.float32)
    with torch.no_grad():
        for offset in range(0, len(frame), config.batch_size):
            batch = x[offset : offset + config.batch_size].to(resolved_device)
            predictions[offset : offset + len(batch)] = model(batch).cpu().numpy()

    metrics = {
        "model_version": MODEL_VERSION,
        "encoder_model_id": ENCODER_MODEL_ID,
        "interpretation": (
            "Teacher-imitation metrics measure NetMHCpan fidelity, not independent "
            "biological accuracy."
        ),
        "best_epoch": best_epoch,
        "best_val_loss": best_val_loss,
        "history": history,
        "splits": {
            name: _metrics(
                frame,
                predictions,
                split == name,
                config.binder_rank_threshold,
            )
            for name in ("train", "val", "test", "challenge")
        },
    }
    return model.cpu(), metrics


def _validated_cv_folds(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    if "cv_fold" not in frame:
        raise ValueError("cross-validation frame is missing column: cv_fold")

    challenge = (
        frame["split"].astype(str).eq("challenge").to_numpy()
        if "split" in frame
        else np.zeros(len(frame), dtype=bool)
    )
    raw_folds = frame["cv_fold"]
    numeric_folds = pd.to_numeric(raw_folds, errors="coerce")
    eligible = ~challenge
    invalid_numeric = eligible & raw_folds.notna().to_numpy() & numeric_folds.isna().to_numpy()
    if invalid_numeric.any():
        raise ValueError("cv_fold must contain integer fold IDs 0..4")
    if numeric_folds[eligible].isna().any():
        raise ValueError("all non-challenge rows must have a cv_fold")

    fold_array = numeric_folds.fillna(-1).to_numpy(dtype=np.float64)
    eligible_values = fold_array[eligible]
    if not np.equal(eligible_values, np.floor(eligible_values)).all():
        raise ValueError("cv_fold must contain integer fold IDs 0..4")
    observed = set(eligible_values.astype(np.int64).tolist())
    if observed != set(CV_FOLDS):
        raise ValueError(f"cv_fold must contain every fold 0..4; observed {sorted(observed)}")
    return fold_array.astype(np.int64), challenge


def _aggregate_held_out_metrics(fold_reports: list[dict]) -> dict:
    aggregate: dict[str, dict] = {}
    for channel in ("el", "ba"):
        channel_reports = [report["metrics"]["splits"]["test"][channel] for report in fold_reports]
        metric_names = sorted(
            set.intersection(*(set(channel_report) for channel_report in channel_reports))
        )
        aggregate[channel] = {}
        for metric_name in metric_names:
            values = np.asarray(
                [channel_report[metric_name] for channel_report in channel_reports],
                dtype=np.float64,
            )
            if np.isfinite(values).all():
                aggregate[channel][metric_name] = {
                    "mean": float(values.mean()),
                    "std": float(values.std(ddof=0)),
                }
    return aggregate


def train_student_cv(
    embeddings: np.ndarray,
    frame: pd.DataFrame,
    config: StudentTrainConfig | None = None,
    device: str | torch.device | None = None,
) -> tuple[dict[int, MHCINetMHCpanStudent], dict]:
    """Train five heads, holding out each ``cv_fold`` once.

    Validation uses the next fold modulo five. The other three folds train the
    head. Rows marked ``split == "challenge"`` remain external to all fitting
    folds and are reported separately by each head.
    """

    config = config or StudentTrainConfig()
    folds, challenge = _validated_cv_folds(frame)
    models: dict[int, MHCINetMHCpanStudent] = {}
    fold_reports: list[dict] = []
    out_of_fold_predictions = np.full((len(frame), config.n_outputs), np.nan, dtype=np.float32)
    embedding_tensor = torch.from_numpy(np.asarray(embeddings, dtype=np.float32))

    for test_fold in CV_FOLDS:
        validation_fold = (test_fold + 1) % len(CV_FOLDS)
        training_folds = [
            fold for fold in CV_FOLDS if fold not in (test_fold, validation_fold)
        ]
        cv_frame = frame.copy()
        split = np.full(len(frame), "challenge", dtype=object)
        split[(~challenge) & np.isin(folds, training_folds)] = "train"
        split[(~challenge) & (folds == validation_fold)] = "val"
        split[(~challenge) & (folds == test_fold)] = "test"
        cv_frame["split"] = split

        fold_config = replace(config, seed=config.seed + test_fold)
        model, fold_metrics = train_student(
            embeddings,
            cv_frame,
            config=fold_config,
            device=device,
        )
        models[test_fold] = model
        test_indices = np.flatnonzero((~challenge) & (folds == test_fold))
        with torch.inference_mode():
            for offset in range(0, len(test_indices), config.batch_size):
                indices = test_indices[offset : offset + config.batch_size]
                out_of_fold_predictions[indices] = model(embedding_tensor[indices]).numpy()
        fold_reports.append(
            {
                "test_fold": test_fold,
                "validation_fold": validation_fold,
                "training_folds": training_folds,
                "metrics": fold_metrics,
            }
        )

    metrics = {
        "model_version": MODEL_VERSION,
        "encoder_model_id": ENCODER_MODEL_ID,
        "validation_scheme": {
            "n_folds": len(CV_FOLDS),
            "validation_fold_rule": "(test_fold + 1) modulo 5",
            "training_folds_per_head": 3,
        },
        "interpretation": (
            "Cross-validation metrics measure held-out NetMHCpan teacher imitation, not "
            "independent biological accuracy."
        ),
        "folds": fold_reports,
        "aggregate": {
            "held_out_test": _aggregate_held_out_metrics(fold_reports),
            "pooled_out_of_fold": _metrics(
                frame,
                out_of_fold_predictions,
                ~challenge,
                config.binder_rank_threshold,
            ),
            "std_definition": "population standard deviation across five held-out folds (ddof=0)",
        },
    }
    return models, metrics


def screening_thresholds(
    propensities: np.ndarray,
    teacher_ranks: np.ndarray,
    *,
    binder_rank: float,
    target_recalls: tuple[float, ...] = (0.80, 0.90, 0.95, 0.99),
) -> dict[str, Any]:
    """Find student cutoffs that hit target recall of the teacher's binders.

    The default class thresholds invert a predicted rank and therefore inherit
    whatever bias the student has at that rank. A screen wants recall instead:
    missing a real binder is worse than flagging an extra peptide for review.
    """

    propensities = np.asarray(propensities, dtype=np.float64)
    teacher_ranks = np.asarray(teacher_ranks, dtype=np.float64)
    if propensities.shape != teacher_ranks.shape:
        raise ValueError("propensities and teacher ranks must align row for row")
    if not np.isfinite(propensities).all():
        raise ValueError("student propensities must be finite")

    positives = teacher_ranks <= binder_rank
    n_positive = int(positives.sum())
    if not n_positive or positives.all():
        raise ValueError("teacher ranks must contain both binders and non-binders")

    order = np.argsort(-propensities, kind="stable")
    sorted_positive = positives[order]
    true_positives = np.cumsum(sorted_positive)
    flagged = np.arange(1, len(order) + 1)
    recall = true_positives / n_positive

    # The default cutoff, for comparison against every recall-driven option.
    default_flagged = propensity_to_rank(propensities) <= binder_rank
    default_true_positives = int((default_flagged & positives).sum())
    operating_points = []
    for target in target_recalls:
        reached = np.flatnonzero(recall >= target)
        if not len(reached):
            continue
        index = int(reached[0])
        operating_points.append(
            {
                "target_recall": float(target),
                "propensity_threshold": float(propensities[order][index]),
                "implied_predicted_rank": float(
                    propensity_to_rank(np.asarray(propensities[order][index]))
                ),
                "recall": float(recall[index]),
                "precision": float(true_positives[index] / flagged[index]),
                "n_flagged": int(flagged[index]),
            }
        )

    return {
        "teacher_binder_rank": float(binder_rank),
        "n_rows": int(len(propensities)),
        "n_teacher_binders": n_positive,
        "default_cutoff": {
            "rule": f"predicted rank <= {binder_rank}",
            "n_flagged": int(default_flagged.sum()),
            "recall": float(default_true_positives / n_positive),
            "precision": float(
                default_true_positives / max(1, int(default_flagged.sum()))
            ),
        },
        "operating_points": operating_points,
    }


def _output_names(n_outputs: int) -> list[str]:
    names = ["netmhcpan_el_rank_propensity", "netmhcpan_ba_rank_propensity"]
    if n_outputs > AFFINITY_INDEX:
        names.append("netmhcpan_ba_affinity_score")
    return names


def save_student_checkpoint(
    output_dir: Path,
    model: MHCINetMHCpanStudent,
    metrics: dict,
    *,
    corpus_sha256: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), output_dir / "weights.pt")
    manifest = {
        "model_version": MODEL_VERSION,
        "encoder_model_id": ENCODER_MODEL_ID,
        "embedding_dim": EMBEDDING_DIM,
        "pooling": "mean over 9 peptide residues; BOS/EOS excluded",
        "hidden_dim": model.hidden_dim,
        "dropout": model.dropout,
        "n_outputs": model.n_outputs,
        "outputs": _output_names(model.n_outputs),
        "profile_formatter": {
            "function": "re_agent.e2e_pls.netmhcpan_student.format_mhci_profile",
            "overall_mhci_risk": "EL presentation propensity only",
            "binder_class_basis": "predicted EL rank",
            "strong_binder_rank_max": STRONG_BINDER_RANK,
            "weak_binder_rank_max": WEAK_BINDER_RANK,
        },
        "corpus_sha256": corpus_sha256,
        "intended_use": "fast HLA-A*02:01 NetMHCpan surrogate and separate MHC-I risk lane",
        "not_ground_truth": True,
        "caveats": list(PROFILE_CAVEATS),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")


def save_student_cv_checkpoints(
    output_dir: Path,
    models: dict[int, MHCINetMHCpanStudent],
    metrics: dict,
    *,
    corpus_sha256: str,
) -> None:
    """Save all fold heads and a deployment manifest for their ensemble."""

    if set(models) != set(CV_FOLDS):
        raise ValueError(f"expected models for folds 0..4, got {sorted(models)}")
    output_dir.mkdir(parents=True, exist_ok=True)

    checkpoints = []
    fold_reports = {report["test_fold"]: report for report in metrics["folds"]}
    for fold in CV_FOLDS:
        relative_path = Path(f"fold_{fold}")
        save_student_checkpoint(
            output_dir / relative_path,
            models[fold],
            fold_reports[fold]["metrics"],
            corpus_sha256=corpus_sha256,
        )
        checkpoints.append(
            {
                "test_fold": fold,
                "validation_fold": fold_reports[fold]["validation_fold"],
                "training_folds": fold_reports[fold]["training_folds"],
                "checkpoint": str(relative_path / "weights.pt"),
                "manifest": str(relative_path / "manifest.json"),
            }
        )

    deployment_manifest = {
        "model_version": MODEL_VERSION,
        "artifact_type": "five_fold_cv_ensemble",
        "encoder_model_id": ENCODER_MODEL_ID,
        "embedding_dim": EMBEDDING_DIM,
        "pooling": "mean over 9 peptide residues; BOS/EOS excluded",
        "n_outputs": models[CV_FOLDS[0]].n_outputs,
        "outputs": _output_names(models[CV_FOLDS[0]].n_outputs),
        "affinity_inversion": (
            f"IC50 nM = {MAX_IC50_NM:.0f} ** (1 - netmhcpan_ba_affinity_score)"
        ),
        "ensemble_method": "arithmetic mean of each fold head for each output channel",
        "profile_formatter": {
            "function": "re_agent.e2e_pls.netmhcpan_student.format_mhci_profile",
            "overall_mhci_risk": "EL presentation propensity only",
            "binder_class_basis": "predicted EL rank",
            "strong_binder_rank_max": STRONG_BINDER_RANK,
            "weak_binder_rank_max": WEAK_BINDER_RANK,
        },
        "checkpoints": checkpoints,
        "corpus_sha256": corpus_sha256,
        "intended_use": "fast HLA-A*02:01 NetMHCpan surrogate and separate MHC-I risk lane",
        "not_ground_truth": True,
        "caveats": list(PROFILE_CAVEATS),
    }
    (output_dir / "deployment_manifest.json").write_text(
        json.dumps(deployment_manifest, indent=2) + "\n"
    )
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
