"""End-to-end orchestration for candidate immunogenicity screening."""

from __future__ import annotations

import json
import uuid
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from re_agent.immuno.contracts import (
    CandidateAssessment,
    MHCRequest,
    ResponseModelRequest,
)
from re_agent.immuno.fusion import apply_fusion_rule, summarize_components
from re_agent.immuno.proxies import processing_evidence, tile_sequence, tolerance_evidence
from re_agent.immuno.registry import MHCProviderRegistry, ResponseModelRegistry
from re_agent.immuno.spatial import build_spatial_tracks


class ImmunogenicityScreeningAgent:
    def __init__(
        self,
        response_registry: ResponseModelRegistry,
        mhc_registry: MHCProviderRegistry,
        default_response_adapter: str,
        hla_panel_path: Path,
        fusion_rule_path: Path,
        self_proteome_path: Path,
    ) -> None:
        self.response_registry = response_registry
        self.mhc_registry = mhc_registry
        self.default_response_adapter = default_response_adapter
        self.hla_panel_path = hla_panel_path
        self.fusion_rule_path = fusion_rule_path
        panel = json.loads(hla_panel_path.read_text())
        self.hla_panel_id = panel["panel_id"]
        self.alleles = panel["alleles"]
        self.self_proteome = pd.read_parquet(self_proteome_path)

    def assess(
        self,
        candidate_id: str,
        sequence: str,
        *,
        accessibility_by_residue: Sequence[float] | None = None,
        cleavage_by_window: Mapping[str, float] | None = None,
        shared_hla_by_window: Mapping[str, Sequence[str]] | None = None,
    ) -> CandidateAssessment:
        request_id = str(uuid.uuid4())
        windows = tile_sequence(sequence)
        response_results = self.response_registry.predict_all(
            ResponseModelRequest(
                request_id=request_id,
                parent_name=candidate_id,
                parent_sequence=sequence,
                windows=windows,
            )
        )
        mhc_results = self.mhc_registry.predict_all(
            MHCRequest(
                request_id=request_id,
                parent_name=candidate_id,
                parent_sequence=sequence,
                alleles=self.alleles,
            )
        )
        processing = processing_evidence(
            windows,
            mhc_results,
            accessibility_by_residue=accessibility_by_residue,
            cleavage_by_window=cleavage_by_window,
        )
        query_hla = self._query_hla_by_window(windows, mhc_results)
        tolerance = tolerance_evidence(
            windows,
            self.self_proteome,
            query_hla,
            shared_hla_by_window=shared_hla_by_window,
        )

        default_result = next(
            (
                result
                for result in response_results
                if result.adapter_id == self.default_response_adapter
            ),
            None,
        )
        warnings: list[str] = []
        if default_result is None:
            warnings.append(
                f"default response adapter {self.default_response_adapter!r} did not run"
            )
            components = {}
            combined = None
            fusion_detail = {}
        else:
            components = summarize_components(default_result, mhc_results, processing, tolerance)
            combined, fusion_detail = apply_fusion_rule(components, self.fusion_rule_path)
        if combined is None:
            warnings.append("combined rank withheld because required evidence is missing")
        if any(row.hla_gate != "shared" for row in tolerance):
            warnings.append("ungated self-similarity is displayed but excluded from fusion")
        spatial_tracks = build_spatial_tracks(
            len(sequence),
            default_result,
            mhc_results,
            processing,
            tolerance,
        )

        return CandidateAssessment(
            candidate_id=candidate_id,
            sequence=sequence,
            response_results=response_results,
            mhc_results=mhc_results,
            processing=processing,
            tolerance=tolerance,
            component_summary={
                "hla_panel_id": self.hla_panel_id,
                "default_response_adapter": self.default_response_adapter,
                "proxy_values": components,
                "fusion": fusion_detail,
            },
            spatial_tracks=spatial_tracks,
            combined_rank_score=combined,
            warnings=warnings,
        )

    @staticmethod
    def _query_hla_by_window(windows, mhc_results) -> dict[str, list[str]]:
        mapping: dict[str, set[str]] = defaultdict(set)
        hits = [hit for result in mhc_results if result.status == "ok" for hit in result.hits]
        for window in windows:
            for hit in hits:
                if window.start < hit.end and (hit.start - 1) < window.end:
                    mapping[window.sequence].add(hit.allele)
        return {sequence: sorted(alleles) for sequence, alleles in mapping.items()}


def rank_assessments(assessments: Sequence[CandidateAssessment]) -> list[CandidateAssessment]:
    """Highest computational risk first; missing combined scores sort last."""
    return sorted(
        assessments,
        key=lambda row: (
            row.combined_rank_score is not None,
            row.combined_rank_score if row.combined_rank_score is not None else -1.0,
        ),
        reverse=True,
    )


def write_assessments(assessments: Sequence[CandidateAssessment], output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": "1.0.0",
        "ranking_direction": "higher computational risk first",
        "methodology": {
            "claim_boundary": (
                "Computational candidate triage from response propensity, MHC-II "
                "presentation/binding, processing, and tolerance proxies; not a measured "
                "immune response probability."
            ),
            "citations": [
                {
                    "claim": "NetMHCIIpan integrates binding-affinity and eluted-ligand data.",
                    "url": (
                        "https://paperclip.gxl.ai/citations/papers/"
                        "PMC7319546#L11,L16,L36-L42"
                    ),
                },
                {
                    "claim": "NetMHCIIpan 4.3 EL and BA are called through the IEDB Tools API.",
                    "url": "https://tools.iedb.org/main/tools-api/",
                },
                {
                    "claim": (
                        "TCR-facing self similarity may support tolerance, subject to a "
                        "shared-HLA presentation gate."
                    ),
                    "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4635734/",
                },
            ],
        },
        "candidates": [row.model_dump(mode="json") for row in rank_assessments(assessments)],
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")
    return output
