"""Versioned contracts shared by response models, MHC providers, and the agent."""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Any, Literal

from pydantic import BaseModel, Field, model_validator

Status = Literal["ok", "unavailable", "unsupported", "error", "timeout"]


class Provenance(BaseModel):
    provider: str
    version: str
    capability: Literal[
        "response_model",
        "mhc_evidence",
        "mhc_i_processing_surrogate",
        "processing",
        "tolerance",
    ]
    source: str
    parameters: dict[str, Any] = Field(default_factory=dict)
    input_sha256: str
    runtime_seconds: float = 0.0
    cached: bool = False
    timestamp: datetime = Field(default_factory=lambda: datetime.now(UTC))


class Window(BaseModel):
    sequence: str = Field(min_length=1, max_length=30)
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    source_n_flank: str = ""
    source_c_flank: str = ""

    @model_validator(mode="after")
    def validate_coordinates(self) -> Window:
        if self.end - self.start != len(self.sequence):
            raise ValueError("window coordinates do not match sequence length")
        return self


class ResponseModelRequest(BaseModel):
    request_id: str
    parent_name: str
    parent_sequence: str
    windows: list[Window]


class ResponsePrediction(BaseModel):
    start: int
    end: int
    sequence: str
    score: float = Field(ge=0.0, le=1.0)
    confidence: float | None = Field(default=None, ge=0.0, le=1.0)
    attribution: list[float] | None = None


class ResponseModelResult(BaseModel):
    adapter_id: str
    status: Status
    predictions: list[ResponsePrediction] = Field(default_factory=list)
    score_scale: Literal["calibrated_probability", "uncalibrated", "unknown"] = "unknown"
    calibration_id: str | None = None
    provenance: Provenance
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None


class MHCRequest(BaseModel):
    request_id: str
    parent_name: str
    parent_sequence: str
    alleles: list[str]


class MHCHit(BaseModel):
    allele: str
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    peptide: str
    core: str
    el_score: float | None = None
    el_rank: float | None = None
    ba_ic50_nm: float | None = None
    ba_rank: float | None = None
    provider_score: float | None = None
    provider_rank: float | None = None


class MHCProviderResult(BaseModel):
    provider_id: str
    status: Status
    hits: list[MHCHit] = Field(default_factory=list)
    supported_alleles: list[str] = Field(default_factory=list)
    missing_alleles: list[str] = Field(default_factory=list)
    provenance: Provenance
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None


class ProcessingEvidence(BaseModel):
    start: int
    end: int
    sequence: str
    el_presentation_support: float | None = Field(default=None, ge=0.0, le=1.0)
    ba_binding_support: float | None = Field(default=None, ge=0.0, le=1.0)
    cleavage_model_score: float | None = Field(default=None, ge=0.0, le=1.0)
    structure_accessibility: float | None = Field(default=None, ge=0.0, le=1.0)
    caveat: str


class ToleranceEvidence(BaseModel):
    start: int
    end: int
    sequence: str
    exact_self_9mer_count: int = Field(ge=0)
    tcr_face_match_count: int = Field(ge=0)
    shared_hla_alleles: list[str] = Field(default_factory=list)
    hla_gate: Literal["shared", "query_only", "not_available"]
    tolerance_score: float = Field(ge=0.0, le=1.0)
    caveat: str


class StructureReference(BaseModel):
    path: str
    format: Literal["pdb"]
    chain_id: str = Field(min_length=1, max_length=4)
    residue_ids: list[str]
    unresolved_sequence_positions: list[int] = Field(default_factory=list)
    sequence_sha256: str
    structure_sha256: str
    mapping_status: Literal[
        "verified_exact_sequence",
        "verified_terminal_trim",
    ]


class MHCISurrogatePrediction(BaseModel):
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    peptide: str = Field(min_length=9, max_length=9)
    cleavage_n_probability: float = Field(ge=0.0, le=1.0)
    cleavage_c_probability: float = Field(ge=0.0, le=1.0)
    tap_log_ic50_relative: float
    tap_uncertainty: float = Field(ge=0.0)
    mhc_i_presentation_propensity: float = Field(ge=0.0, le=1.0)
    composite_processing_risk: float = Field(ge=0.0, le=1.0)
    confidence: float = Field(ge=0.0, le=1.0)


class MHCISurrogateResult(BaseModel):
    adapter_id: str
    status: Status
    allele: str
    predictions: list[MHCISurrogatePrediction] = Field(default_factory=list)
    protein_summary: dict[str, Any] = Field(default_factory=dict)
    spatial_tracks: dict[str, list[float]] = Field(default_factory=dict)
    provenance: Provenance
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None


class CandidateAssessment(BaseModel):
    schema_version: str = "1.0.0"
    candidate_id: str
    sequence: str
    response_results: list[ResponseModelResult]
    mhc_results: list[MHCProviderResult]
    processing: list[ProcessingEvidence]
    tolerance: list[ToleranceEvidence]
    component_summary: dict[str, Any]
    mhc_i_surrogate_results: list[MHCISurrogateResult] = Field(default_factory=list)
    spatial_tracks: dict[str, list[float]] = Field(default_factory=dict)
    structure: StructureReference | None = None
    combined_rank_score: float | None = None
    warnings: list[str] = Field(default_factory=list)
