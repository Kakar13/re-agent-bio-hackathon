"""Inference adapter for the team's frozen E2E-PLS pickle bundle.

This is an MHC-I antigen-processing surrogate, not a CD4 response model. Its
outputs are retained as a separate lane and are never fused into MHC-II risk.
"""

from __future__ import annotations

import hashlib
import math
import pickle
import pickletools
import re
import time
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import torch
from torch import nn

from re_agent.e2e_pls.netmhcpan_student import (
    AFFINITY_INDEX,
    format_mhci_profile,
    load_student_ensemble,
    predict_student_ensemble,
)
from re_agent.e2e_pls.score import compute_pathway_rank_score
from re_agent.immuno.adapters import sequence_sha256
from re_agent.immuno.contracts import (
    MHCISurrogatePrediction,
    MHCISurrogateResult,
    Provenance,
)

PICKLE_FORMAT = "e2e_pls-heads-pickle-v1"
ENCODER_MODEL_ID = "esm2_t33_650M_UR50D"
EMBEDDING_DIM = 1280
PEPTIDE_LENGTH = 9
FLANK_LENGTH = 4
ALLOWED_GLOBALS = {
    ("torch._utils", "_rebuild_tensor_v2"),
    ("torch.storage", "_load_from_bytes"),
    ("collections", "OrderedDict"),
    ("numpy._core.numeric", "_frombuffer"),
    ("numpy", "dtype"),
}


class UnsafeCheckpointError(ValueError):
    """Raised when a pickle references constructors outside the frozen allowlist."""


class _RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> Any:
        if (module, name) not in ALLOWED_GLOBALS:
            raise UnsafeCheckpointError(f"checkpoint references forbidden global: {module}.{name}")
        return super().find_class(module, name)


def _checkpoint_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _audit_pickle(path: Path) -> None:
    if path.stat().st_size > 50 * 1024 * 1024:
        raise UnsafeCheckpointError("checkpoint exceeds the 50 MiB inference limit")
    operations = list(pickletools.genops(path.read_bytes()))
    for index, (operation, argument, _position) in enumerate(operations):
        if operation.name == "GLOBAL":
            module, name = str(argument).split(" ", maxsplit=1)
            if (module, name) not in ALLOWED_GLOBALS:
                raise UnsafeCheckpointError(f"checkpoint references {module}.{name}")
        if operation.name == "STACK_GLOBAL":
            strings = [
                previous_argument
                for previous_operation, previous_argument, _ in operations[
                    max(0, index - 8) : index
                ]
                if previous_operation.name in {"SHORT_BINUNICODE", "BINUNICODE", "UNICODE"}
            ]
            if len(strings) < 2 or tuple(strings[-2:]) not in ALLOWED_GLOBALS:
                raise UnsafeCheckpointError("checkpoint contains an unresolved STACK_GLOBAL")


def _as_finite_array(value: Any, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.shape != shape:
        raise ValueError(f"{name} expected shape {shape}, got {array.shape}")
    if not np.issubdtype(array.dtype, np.number) or not np.isfinite(array).all():
        raise ValueError(f"{name} must contain finite numeric values")
    return array


def _validate_calibrator(payload: dict[str, Any], name: str) -> None:
    x_values = np.asarray(payload["x_thresholds"], dtype=np.float64)
    y_values = np.asarray(payload["y_thresholds"], dtype=np.float64)
    if x_values.ndim != 1 or len(x_values) < 2 or x_values.shape != y_values.shape:
        raise ValueError(f"{name} calibrator thresholds are malformed")
    if not np.isfinite(x_values).all() or not np.isfinite(y_values).all():
        raise ValueError(f"{name} calibrator contains non-finite values")
    if np.any(np.diff(x_values) < 0) or np.any(np.diff(y_values) < 0):
        raise ValueError(f"{name} calibrator must be monotonic")
    if np.any((y_values < 0) | (y_values > 1)):
        raise ValueError(f"{name} calibrated probabilities must lie in [0, 1]")


def load_validated_bundle(path: Path) -> tuple[dict[str, Any], str]:
    """Audit, deserialize, and shape-check the trusted team checkpoint."""

    _audit_pickle(path)
    # This repository-owned team checkpoint is intentionally pickle-based. The
    # restricted unpickler and exact shape checks run before any inference use.
    with path.open("rb") as handle:
        bundle = _RestrictedUnpickler(handle).load()
    if not isinstance(bundle, dict) or bundle.get("format") != PICKLE_FORMAT:
        raise ValueError(f"checkpoint format must be {PICKLE_FORMAT!r}")

    manifest = bundle["manifest"]
    if manifest.get("encoder_model_id") != ENCODER_MODEL_ID:
        raise ValueError("checkpoint encoder does not match ESM-2 t33 650M")

    cleavage = bundle["cleavage"]
    if int(cleavage["hidden_dim"]) != 128:
        raise ValueError("only the validated 128-unit cleavage head is supported")
    expected_state_shapes = {
        "net.0.weight": (128, EMBEDDING_DIM * 2),
        "net.0.bias": (128,),
        "net.2.weight": (64, 128),
        "net.2.bias": (64,),
        "net.4.weight": (2, 64),
        "net.4.bias": (2,),
    }
    state_dict = cleavage["state_dict"]
    if set(state_dict) != set(expected_state_shapes):
        raise ValueError("cleavage state_dict keys do not match the frozen architecture")
    for key, shape in expected_state_shapes.items():
        _as_finite_array(state_dict[key], shape, f"cleavage.{key}")
    _validate_calibrator(cleavage["calibrator_n"], "cleavage_n")
    _validate_calibrator(cleavage["calibrator_c"], "cleavage_c")

    tap = bundle["tap"]
    _as_finite_array(tap["coef"], (EMBEDDING_DIM,), "tap.coef")
    _as_finite_array(tap["bootstrap_coef"], (50, EMBEDDING_DIM), "tap.bootstrap_coef")
    _as_finite_array(tap["bootstrap_intercept"], (50,), "tap.bootstrap_intercept")

    mhc = bundle["mhc"]
    _as_finite_array(mhc["projection_weight"], (64, EMBEDDING_DIM), "mhc.projection_weight")
    _as_finite_array(mhc["projection_bias"], (64,), "mhc.projection_bias")
    if set(mhc["centroids"]) != {"HLA-A*02:01"}:
        raise ValueError("checkpoint must expose only its calibrated HLA-A*02:01 centroid")
    _as_finite_array(mhc["centroids"]["HLA-A*02:01"], (64,), "mhc.centroid")
    _validate_calibrator(mhc["calibrators"]["HLA-A*02:01"], "mhc.HLA-A*02:01")
    return bundle, _checkpoint_sha256(path)


class _CleavageMLP(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(EMBEDDING_DIM * 2, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 2),
        )

    def forward(self, n_vectors: torch.Tensor, c_vectors: torch.Tensor) -> torch.Tensor:
        return self.net(torch.cat([n_vectors, c_vectors], dim=-1))


class _ESM2SegmentEncoder:
    _shared_model: ClassVar[Any | None] = None
    _shared_batch_converter: ClassVar[Any | None] = None
    _shared_device: ClassVar[torch.device | None] = None
    _last_sequences: ClassVar[tuple[str, ...] | None] = None
    _last_embeddings: ClassVar[list[np.ndarray] | None] = None

    def __init__(self, *, batch_size: int = 32) -> None:
        if self.__class__._shared_model is None:
            import esm

            model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
            model.eval()
            if torch.cuda.is_available():
                device = torch.device("cuda")
            elif torch.backends.mps.is_available():
                device = torch.device("mps")
            else:
                device = torch.device("cpu")
            self.__class__._shared_model = model.to(device)
            self.__class__._shared_device = device
            self.__class__._shared_batch_converter = alphabet.get_batch_converter()
        self.model = self.__class__._shared_model
        self.device = self.__class__._shared_device
        self.batch_converter = self.__class__._shared_batch_converter
        self.batch_size = batch_size

    def embed(self, sequences: list[str]) -> list[np.ndarray]:
        cache_key = tuple(sequences)
        if (
            cache_key == self.__class__._last_sequences
            and self.__class__._last_embeddings is not None
        ):
            return self.__class__._last_embeddings
        outputs: list[np.ndarray] = []
        for offset in range(0, len(sequences), self.batch_size):
            batch = sequences[offset : offset + self.batch_size]
            _, _, tokens = self.batch_converter(
                [(f"segment-{offset + index}", sequence) for index, sequence in enumerate(batch)]
            )
            tokens = tokens.to(self.device)
            with torch.no_grad():
                result = self.model(tokens, repr_layers=[33], return_contacts=False)
            representations = result["representations"][33].float().cpu()
            outputs.extend(
                representations[index, 1 : 1 + len(sequence)].numpy()
                for index, sequence in enumerate(batch)
            )
        self.__class__._last_sequences = cache_key
        self.__class__._last_embeddings = outputs
        return outputs


def _interp(calibrator: dict[str, Any], values: np.ndarray) -> np.ndarray:
    return np.interp(values, calibrator["x_thresholds"], calibrator["y_thresholds"])


def _confidence(
    tap_uncertainty: float,
    mhc_propensity: float,
    cleavage_n: float,
    cleavage_c: float,
    n_flank_length: int,
    c_flank_length: int,
) -> float:
    tap = 1.0 / (1.0 + max(tap_uncertainty, 0.0))
    mhc = float(np.clip(2.0 * abs(mhc_propensity - 0.5), 0.0, 1.0))
    context = 0.5 * (
        min(n_flank_length, FLANK_LENGTH) / FLANK_LENGTH
        + min(c_flank_length, FLANK_LENGTH) / FLANK_LENGTH
    )
    probabilities = np.asarray([cleavage_n, cleavage_c, mhc_propensity])
    maximum_std = float(np.std([0.0, 0.0, 1.0]))
    agreement = float(np.clip(1.0 - np.std(probabilities) / maximum_std, 0.0, 1.0))
    return float(np.exp(np.mean(np.log(np.clip([tap, mhc, context, agreement], 1e-6, 1.0)))))


class TeamE2EPLSAdapter:
    """Score one parent sequence with the frozen MHC-I processing surrogate."""

    def __init__(
        self,
        checkpoint_path: Path,
        *,
        netmhcpan_checkpoint_dir: Path | None = None,
        adapter_id: str | None = None,
        batch_size: int = 32,
    ) -> None:
        self.checkpoint_path = checkpoint_path
        model_name = re.sub(r"[^a-z0-9-]+", "-", checkpoint_path.parent.name.lower()).strip("-")
        self.netmhcpan_checkpoint_dir = netmhcpan_checkpoint_dir
        self.netmhcpan_students = None
        self.netmhcpan_checkpoint_sha256 = None
        self.netmhcpan_manifest = None
        if netmhcpan_checkpoint_dir is not None:
            (
                self.netmhcpan_students,
                self.netmhcpan_checkpoint_sha256,
                self.netmhcpan_manifest,
            ) = load_student_ensemble(netmhcpan_checkpoint_dir)
        suffix = "-netmhcpan" if self.netmhcpan_students is not None else ""
        self.adapter_id = adapter_id or f"team-e2e-pls-{model_name}{suffix}"
        self.bundle, self.checkpoint_sha256 = load_validated_bundle(checkpoint_path)
        self.version = self.bundle["manifest"]["model_version"]
        self.dataset_version_hash = self.bundle["manifest"]["dataset_version_hash"]
        self.batch_size = batch_size
        self._encoder: _ESM2SegmentEncoder | None = None

        self.cleavage_model = _CleavageMLP()
        self.cleavage_model.load_state_dict(self.bundle["cleavage"]["state_dict"], strict=True)
        self.cleavage_model.eval()

    def predict(self, sequence: str, allele: str = "HLA-A*02:01") -> MHCISurrogateResult:
        started = time.perf_counter()
        sequence = "".join(sequence.split()).upper()
        input_hash = sequence_sha256(sequence)
        try:
            if len(sequence) < PEPTIDE_LENGTH:
                raise ValueError("sequence must contain at least one 9-mer")
            if any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in sequence):
                raise ValueError("sequence contains non-canonical residues")
            if allele != "HLA-A*02:01":
                raise ValueError("checkpoint is calibrated only for HLA-A*02:01")
            predictions, tracks, summary = self._score(sequence, allele)
            status = "ok"
            error = None
        except Exception as exc:
            predictions, tracks, summary = [], {}, {}
            status = "error"
            error = f"{type(exc).__name__}: {exc}"

        return MHCISurrogateResult(
            adapter_id=self.adapter_id,
            status=status,
            allele=allele,
            predictions=predictions,
            protein_summary=summary,
            spatial_tracks=tracks,
            provenance=Provenance(
                provider=self.adapter_id,
                version=self.version,
                capability="mhc_i_processing_surrogate",
                source=str(self.checkpoint_path),
                parameters={
                    "checkpoint_sha256": self.checkpoint_sha256,
                    "dataset_version_hash": self.dataset_version_hash,
                    "encoder_model_id": ENCODER_MODEL_ID,
                    "peptide_length": PEPTIDE_LENGTH,
                    "flank_length": FLANK_LENGTH,
                    "allele": allele,
                    "netmhcpan_student_checkpoint": (
                        str(self.netmhcpan_checkpoint_dir)
                        if self.netmhcpan_checkpoint_dir is not None
                        else None
                    ),
                    "netmhcpan_student_sha256": self.netmhcpan_checkpoint_sha256,
                },
                input_sha256=input_hash,
                runtime_seconds=time.perf_counter() - started,
            ),
            warnings=[
                "MHC-I processing surrogate; do not interpret as MHC-II or CD4 response.",
                "Composite processing risk is not a measured immune-response probability.",
                "PDA smoke inference may overlap the model's de novo training distribution.",
                *(
                    [
                        "Cleavage and TAP use the legacy checkpoint; MHC-I EL/BA use the "
                        "PDA-trained NetMHCpan student ensemble.",
                        "NetMHCpan student outputs measure teacher imitation, not experimental "
                        "binding or immune response.",
                    ]
                    if self.netmhcpan_students is not None
                    else []
                ),
            ],
            error=error,
        )

    def _score(
        self,
        sequence: str,
        allele: str,
    ) -> tuple[list[MHCISurrogatePrediction], dict[str, list[float]], dict[str, Any]]:
        windows = []
        segments = []
        for start in range(len(sequence) - PEPTIDE_LENGTH + 1):
            end = start + PEPTIDE_LENGTH
            n_flank = sequence[max(0, start - FLANK_LENGTH) : start]
            c_flank = sequence[end : end + FLANK_LENGTH]
            windows.append((start, end, sequence[start:end], n_flank, c_flank))
            segments.append(f"{n_flank}{sequence[start:end]}{c_flank}")

        if self._encoder is None:
            self._encoder = _ESM2SegmentEncoder(batch_size=self.batch_size)
        embeddings = self._encoder.embed(segments)
        n_vectors = []
        c_vectors = []
        mer_vectors = []
        for embedding, (_start, _end, peptide, n_flank, _c_flank) in zip(
            embeddings, windows, strict=True
        ):
            n_length = len(n_flank)
            n_vectors.append(
                embedding[max(0, n_length - 3) : min(len(embedding), n_length + 3)].mean(axis=0)
            )
            c_site = n_length + len(peptide)
            c_vectors.append(
                embedding[max(0, c_site - 3) : min(len(embedding), c_site + 3)].mean(axis=0)
            )
            mer_vectors.append(embedding[n_length:c_site].mean(axis=0))

        n_array = np.stack(n_vectors).astype("float32")
        c_array = np.stack(c_vectors).astype("float32")
        mer_array = np.stack(mer_vectors)
        with torch.no_grad():
            logits = self.cleavage_model(
                torch.from_numpy(n_array),
                torch.from_numpy(c_array),
            ).numpy()
        probabilities = 1.0 / (1.0 + np.exp(-logits))
        cleavage = self.bundle["cleavage"]
        cleavage_n = _interp(cleavage["calibrator_n"], probabilities[:, 0])
        cleavage_c = _interp(cleavage["calibrator_c"], probabilities[:, 1])

        tap = self.bundle["tap"]
        tap_mean = mer_array @ tap["coef"] + float(tap["intercept"])
        tap_bootstrap = mer_array @ tap["bootstrap_coef"].T + tap["bootstrap_intercept"]
        tap_uncertainty = tap_bootstrap.std(axis=1)

        student_profiles = None
        if self.netmhcpan_students is not None:
            student_predictions = predict_student_ensemble(self.netmhcpan_students, mer_array)
            has_affinity = student_predictions.shape[1] > AFFINITY_INDEX
            student_profiles = [
                format_mhci_profile(
                    float(row[0]),
                    float(row[1]),
                    float(row[AFFINITY_INDEX]) if has_affinity else None,
                )
                for row in student_predictions
            ]
            mhc_propensity = student_predictions[:, 0]
        else:
            mhc = self.bundle["mhc"]
            projected = mer_array @ mhc["projection_weight"].T + mhc["projection_bias"]
            projected /= np.clip(np.linalg.norm(projected, axis=1, keepdims=True), 1e-6, None)
            cosine = projected @ mhc["centroids"][allele]
            mhc_propensity = _interp(mhc["calibrators"][allele], cosine)

        prediction_rows = []
        for index, (start, end, peptide, n_flank, c_flank) in enumerate(windows):
            profile = student_profiles[index] if student_profiles is not None else None
            binding_propensity = (
                float(profile["ba_binding_propensity"])
                if profile is not None
                else float(mhc_propensity[index])
            )
            composite = float(
                math.exp(
                    np.mean(
                        np.log(
                            np.clip(
                                [cleavage_n[index], cleavage_c[index], binding_propensity],
                                1e-6,
                                1.0,
                            )
                        )
                    )
                )
            )
            el_propensity = (
                float(profile["el_presentation_propensity"])
                if profile is not None
                else float(mhc_propensity[index])
            )
            pathway_rank = compute_pathway_rank_score(
                el_propensity,
                float(cleavage_n[index]),
                float(cleavage_c[index]),
            )
            prediction_rows.append(
                MHCISurrogatePrediction(
                    start=start,
                    end=end,
                    peptide=peptide,
                    cleavage_n_probability=float(cleavage_n[index]),
                    cleavage_c_probability=float(cleavage_c[index]),
                    tap_log_ic50_relative=float(tap_mean[index]),
                    tap_uncertainty=float(tap_uncertainty[index]),
                    mhc_i_presentation_propensity=float(mhc_propensity[index]),
                    mhc_i_binding_propensity=(
                        profile["ba_binding_propensity"] if profile is not None else None
                    ),
                    mhc_i_predicted_el_rank=(
                        profile["predicted_el_rank"] if profile is not None else None
                    ),
                    mhc_i_predicted_ba_rank=(
                        profile["predicted_ba_rank"] if profile is not None else None
                    ),
                    mhc_i_predicted_ba_ic50_nm=(
                        profile.get("predicted_ba_ic50_nm") if profile is not None else None
                    ),
                    mhc_i_screening_flag=(
                        profile.get("screening_flag") if profile is not None else None
                    ),
                    mhc_i_binder_class=(
                        profile["binder_class"] if profile is not None else None
                    ),
                    mhc_i_risk_band=profile["risk_band"] if profile is not None else None,
                    overall_mhci_risk=(
                        profile["overall_mhci_risk"] if profile is not None else None
                    ),
                    composite_processing_risk=composite,
                    pathway_rank_score=pathway_rank,
                    confidence=_confidence(
                        float(tap_uncertainty[index]),
                        binding_propensity,
                        float(cleavage_n[index]),
                        float(cleavage_c[index]),
                        len(n_flank),
                        len(c_flank),
                    ),
                )
            )

        tracks = self._spatial_tracks(len(sequence), prediction_rows)
        ranks = np.asarray(
            [
                row.pathway_rank_score
                if row.pathway_rank_score is not None
                else row.composite_processing_risk
                for row in prediction_rows
            ]
        )
        order = np.argsort(-ranks)
        top_count = min(5, len(ranks))
        summary = {
            "n_windows": len(prediction_rows),
            "top_k": top_count,
            "top_k_mean_risk": float(ranks[order[:top_count]].mean()),
            "max_risk": float(ranks[order[0]]),
            "max_risk_window_start": int(prediction_rows[int(order[0])].start),
            "mean_confidence": float(
                np.mean([row.confidence for row in prediction_rows])
            ),
        }
        if student_profiles is not None:
            summary["risk_profile"] = {
                "processing": {
                    "max_cleavage_n_probability": float(
                        max(row.cleavage_n_probability for row in prediction_rows)
                    ),
                    "max_cleavage_c_probability": float(
                        max(row.cleavage_c_probability for row in prediction_rows)
                    ),
                },
                "transport": {
                    "mean_tap_log_ic50_relative": float(
                        np.mean([row.tap_log_ic50_relative for row in prediction_rows])
                    ),
                    "mean_tap_uncertainty": float(
                        np.mean([row.tap_uncertainty for row in prediction_rows])
                    ),
                },
                "mhc_i": {
                    "max_el_presentation_propensity": float(
                        max(row.mhc_i_presentation_propensity for row in prediction_rows)
                    ),
                    "max_ba_binding_propensity": float(
                        max(
                            row.mhc_i_binding_propensity or 0.0
                            for row in prediction_rows
                        )
                    ),
                    "highest_binder_class": next(
                        (
                            binder_class
                            for binder_class in ("strong", "weak", "nonbinder")
                            if any(
                                row.mhc_i_binder_class == binder_class
                                for row in prediction_rows
                            )
                        ),
                        "nonbinder",
                    ),
                },
                "overall": {
                    "score": float(risks[order[:top_count]].mean()),
                    "definition": (
                        "top-window mean geometric processing score from N-cleavage, "
                        "C-cleavage, and BA binding; TAP and EL are reported separately"
                    ),
                },
            }
        return prediction_rows, tracks, summary

    @staticmethod
    def _spatial_tracks(
        sequence_length: int,
        predictions: list[MHCISurrogatePrediction],
    ) -> dict[str, list[float]]:
        risk_values: list[list[float]] = [[] for _ in range(sequence_length)]
        mhc_values: list[list[float]] = [[] for _ in range(sequence_length)]
        binding_values: list[list[float]] = [[] for _ in range(sequence_length)]
        for prediction in predictions:
            for position in range(prediction.start, prediction.end):
                risk_values[position].append(prediction.composite_processing_risk)
                mhc_values[position].append(prediction.mhc_i_presentation_propensity)
                if prediction.mhc_i_binding_propensity is not None:
                    binding_values[position].append(prediction.mhc_i_binding_propensity)
        tracks = {
            "mhci_processing_risk_max": [
                max(values) if values else 0.0 for values in risk_values
            ],
            "mhci_presentation_propensity_mean": [
                float(np.mean(values)) if values else 0.0 for values in mhc_values
            ],
        }
        if any(binding_values):
            tracks["mhci_binding_propensity_mean"] = [
                float(np.mean(values)) if values else 0.0 for values in binding_values
            ]
        return tracks
