"""Shared schemas for immuno-risk predictors and pipeline outputs."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field


class PeptideHit(BaseModel):
    peptide: str
    allele: str
    mhc_class: Literal["I", "II"]
    start: int | None = None
    end: int | None = None
    length: int
    affinity_nm: float | None = None
    presentation_score: float | None = None
    processing_score: float | None = None
    percentile_rank: float | None = None
    el_score: float | None = None
    ba_score: float | None = None
    pathogen_epitope_score: float | None = None
    neoepitope_score: float | None = None
    binding_core: str | None = None
    distance_to_training: float | None = None
    binder: bool = False
    method: str
    version: str | None = None
    provenance: dict[str, Any] = Field(default_factory=dict)
    caveat: str = ""


class ToleranceEvidence(BaseModel):
    peptide: str
    allele: str | None = None
    status: Literal["self_like", "foreign_like", "unknown"]
    nearest_self: str | None = None
    identity: float = 0.0
    atlas_hit: bool = False
    method: str
    caveat: str = ""


class PeptideEvidence(BaseModel):
    peptide: str
    allele: str
    mhc_class: Literal["I", "II"]
    start: int | None = None
    end: int | None = None
    mhc: PeptideHit
    tolerance: ToleranceEvidence | None = None
    presentation_points: int = 0
    tolerance_points: int = 0
    point_score: int = 0
    contribution: float = 0.0


class ResidueRisk(BaseModel):
    position: int  # 1-based
    residue: str
    risk: float
    peptide_count: int
    peptides: list[str] = Field(default_factory=list)


class AggregationReport(BaseModel):
    sequence_id: str
    overall: Literal["low", "moderate", "high"]
    score0to100: float
    factors: list[dict[str, Any]] = Field(default_factory=list)
    net_charge_ph74: float
    hydrophobic_fraction: float
    free_cysteine_count: int
    beta_edge_proxy: float
    solubility_proxy: float
    persistence_caveat: str
    method: str = "lightweight_aggregation_v1"
    caveat: str = (
        "Separate from epitope risk. Do not infer aggregation from protease accessibility alone."
    )


class ArmSummary(BaseModel):
    mhc_class: Literal["I", "II"]
    binder_count: int
    top_hits: list[PeptideHit] = Field(default_factory=list)
    total_points: int = 0
    max_peptide_points: int = 0
    method: str
    caveat: str = ""


class ConfidenceReport(BaseModel):
    score0to1: float
    factors: list[dict[str, Any]] = Field(default_factory=list)
    method: str = "coverage_agreement_v1"


class RiskBreakdown(BaseModel):
    sequence_id: str
    overall: Literal["low", "moderate", "high"]
    score0to100: float
    mhc_i: ArmSummary | None = None
    mhc_ii: ArmSummary | None = None
    total_points: int = 0
    max_points: int = 85
    factors: list[dict[str, Any]] = Field(default_factory=list)
    peptides_flagged: list[str] = Field(default_factory=list)
    method: str
    caveat: str


class ImmunoRunResult(BaseModel):
    run_id: str
    sequence_id: str
    sequence: str
    delivery_mode: str = "intracellular_plasmid"
    alleles_i: list[str] = Field(default_factory=list)
    alleles_ii: list[str] = Field(default_factory=list)
    peptides: list[PeptideEvidence] = Field(default_factory=list)
    risk: RiskBreakdown
    confidence: ConfidenceReport
    residue_risk: list[ResidueRisk] = Field(default_factory=list)
    aggregation: AggregationReport
    predictor_versions: dict[str, str] = Field(default_factory=dict)
    artifact_dir: str | None = None
    caveats: list[str] = Field(default_factory=list)
