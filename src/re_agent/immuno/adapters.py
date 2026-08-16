"""Adapters for teammate artifacts and externally generated MHC evidence.

These adapters intentionally consume versioned artifacts. They do not train a
model, silently substitute a heuristic, or fabricate values when a provider is
unavailable.
"""

from __future__ import annotations

import hashlib
import json
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

from re_agent.immuno.contracts import (
    MHCHit,
    MHCProviderResult,
    MHCRequest,
    Provenance,
    ResponseModelRequest,
    ResponseModelResult,
    ResponsePrediction,
)


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


class CallableResponseAdapter:
    """Wrap a teammate predictor without constraining its internal model stack."""

    def __init__(
        self,
        adapter_id: str,
        version: str,
        predict_fn: Callable[[ResponseModelRequest], list[ResponsePrediction]],
        source: str,
        parameters: dict[str, Any] | None = None,
        score_scale: str = "unknown",
        calibration_id: str | None = None,
    ) -> None:
        self.adapter_id = adapter_id
        self.version = version
        self._predict_fn = predict_fn
        self.source = source
        self.parameters = parameters or {}
        self.score_scale = score_scale
        self.calibration_id = calibration_id

    def predict(self, request: ResponseModelRequest) -> ResponseModelResult:
        started = time.perf_counter()
        try:
            predictions = self._predict_fn(request)
            status = "ok"
            error = None
        except Exception as exc:  # provider boundary: preserve the rest of the agent run
            predictions = []
            status = "error"
            error = f"{type(exc).__name__}: {exc}"
        provenance = Provenance(
            provider=self.adapter_id,
            version=self.version,
            capability="response_model",
            source=self.source,
            parameters=self.parameters,
            input_sha256=sequence_sha256(request.parent_sequence),
            runtime_seconds=time.perf_counter() - started,
        )
        return ResponseModelResult(
            adapter_id=self.adapter_id,
            status=status,
            predictions=predictions,
            score_scale=self.score_scale,
            calibration_id=self.calibration_id,
            provenance=provenance,
            error=error,
        )


class ArtifactResponseAdapter:
    """Read immutable predictions keyed by sequence hash from a JSON artifact."""

    def __init__(self, path: Path) -> None:
        self.path = path
        payload = json.loads(path.read_text())
        self.adapter_id = payload["model_card"]["adapter_id"]
        self.version = payload["model_card"]["version"]
        self.model_card = payload["model_card"]
        self._predictions = payload["predictions_by_sha256"]

    def predict(self, request: ResponseModelRequest) -> ResponseModelResult:
        key = sequence_sha256(request.parent_sequence)
        rows = self._predictions.get(key)
        status = "ok" if rows is not None else "unsupported"
        predictions = [ResponsePrediction.model_validate(row) for row in rows or []]
        return ResponseModelResult(
            adapter_id=self.adapter_id,
            status=status,
            predictions=predictions,
            score_scale=self.model_card.get("score_scale", "unknown"),
            calibration_id=self.model_card.get("calibration_id"),
            provenance=Provenance(
                provider=self.adapter_id,
                version=self.version,
                capability="response_model",
                source=str(self.path),
                parameters=self.model_card.get("parameters", {}),
                input_sha256=key,
                cached=True,
            ),
            warnings=[] if rows is not None else ["sequence is absent from prediction artifact"],
        )


class UnavailableResponseAdapter:
    """Explicit missing-provider lane used while waiting for a teammate checkpoint."""

    version = "unavailable"

    def __init__(self, adapter_id: str, reason: str) -> None:
        self.adapter_id = adapter_id
        self.reason = reason

    def predict(self, request: ResponseModelRequest) -> ResponseModelResult:
        return ResponseModelResult(
            adapter_id=self.adapter_id,
            status="unavailable",
            provenance=Provenance(
                provider=self.adapter_id,
                version=self.version,
                capability="response_model",
                source="provider handoff",
                parameters={},
                input_sha256=sequence_sha256(request.parent_sequence),
            ),
            warnings=[self.reason],
        )


class NetMHCIIpanArtifactProvider:
    """Load cached NetMHCIIpan rows while preserving EL and BA as separate evidence."""

    provider_id = "netmhciipan"

    def __init__(self, path: Path, version: str | None = None) -> None:
        self.path = path
        payload = json.loads(path.read_text())
        artifact_version = payload["provenance"]["version"]
        if version is not None and artifact_version != version:
            raise ValueError(f"expected NetMHCIIpan {version}, artifact is {artifact_version}")
        self.version = artifact_version
        self._rows = payload["predictions_by_sha256"]
        self._parameters = payload["provenance"].get("parameters", {})

    def predict(self, request: MHCRequest) -> MHCProviderResult:
        key = sequence_sha256(request.parent_sequence)
        rows = self._rows.get(key)
        if rows is None:
            return MHCProviderResult(
                provider_id=self.provider_id,
                status="unsupported",
                provenance=self._provenance(key),
                warnings=["sequence is absent from the NetMHCIIpan cache"],
            )

        requested = set(request.alleles)
        hits = [MHCHit.model_validate(row) for row in rows if row["allele"] in requested]
        supported = sorted({hit.allele for hit in hits})
        return MHCProviderResult(
            provider_id=self.provider_id,
            status="ok" if hits else "unsupported",
            hits=hits,
            supported_alleles=supported,
            missing_alleles=sorted(requested - set(supported)),
            provenance=self._provenance(key),
            warnings=[
                "EL and BA are parallel evidence channels; do not average their raw scales."
            ],
        )

    def _provenance(self, key: str) -> Provenance:
        return Provenance(
            provider=self.provider_id,
            version=self.version,
            capability="mhc_evidence",
            source=str(self.path),
            parameters=self._parameters,
            input_sha256=key,
            cached=True,
        )


class ChallengerArtifactProvider:
    """Generic cached adapter for MixMHC2pred or Graph-pMHC outputs."""

    def __init__(self, provider_id: str, version: str, path: Path) -> None:
        self.provider_id = provider_id
        self.version = version
        self.path = path
        payload = json.loads(path.read_text())
        self._rows = payload["predictions_by_sha256"]
        self._parameters = payload.get("provenance", {}).get("parameters", {})

    def predict(self, request: MHCRequest) -> MHCProviderResult:
        key = sequence_sha256(request.parent_sequence)
        rows = self._rows.get(key)
        if rows is None:
            return MHCProviderResult(
                provider_id=self.provider_id,
                status="unsupported",
                provenance=self._provenance(key),
                warnings=["sequence is absent from challenger cache"],
            )
        requested = set(request.alleles)
        hits = [MHCHit.model_validate(row) for row in rows if row["allele"] in requested]
        supported = sorted({hit.allele for hit in hits})
        return MHCProviderResult(
            provider_id=self.provider_id,
            status="ok" if hits else "unsupported",
            hits=hits,
            supported_alleles=supported,
            missing_alleles=sorted(requested - set(supported)),
            provenance=self._provenance(key),
        )

    def _provenance(self, key: str) -> Provenance:
        return Provenance(
            provider=self.provider_id,
            version=self.version,
            capability="mhc_evidence",
            source=str(self.path),
            parameters=self._parameters,
            input_sha256=key,
            cached=True,
        )
