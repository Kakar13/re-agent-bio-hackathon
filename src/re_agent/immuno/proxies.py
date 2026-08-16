"""Separate, inspectable processing and tolerance proxy channels."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import pandas as pd

from re_agent.immuno.contracts import (
    MHCProviderResult,
    ProcessingEvidence,
    ToleranceEvidence,
    Window,
)


def tile_sequence(sequence: str, width: int = 15, stride: int = 1) -> list[Window]:
    if len(sequence) < width:
        return [Window(sequence=sequence, start=0, end=len(sequence))]
    return [
        Window(sequence=sequence[start : start + width], start=start, end=start + width)
        for start in range(0, len(sequence) - width + 1, stride)
    ]


def _rank_support(rank_percent: float | None) -> float | None:
    if rank_percent is None:
        return None
    return 1.0 - min(max(float(rank_percent), 0.0), 100.0) / 100.0


def processing_evidence(
    windows: Sequence[Window],
    mhc_results: Sequence[MHCProviderResult],
    accessibility_by_residue: Sequence[float] | None = None,
    cleavage_by_window: Mapping[str, float] | None = None,
) -> list[ProcessingEvidence]:
    """Derive processing evidence without collapsing EL and BA outputs.

    NetMHCIIpan EL is used as empirical presentation support. BA is retained as
    binding evidence, not mislabeled as cleavage. A separate cleavage score must
    come from an explicit provider if one is available.
    """
    hits = [hit for result in mhc_results if result.status == "ok" for hit in result.hits]
    output = []
    for window in windows:
        overlapping = [
            hit
            for hit in hits
            if window.start < hit.end and (hit.start - 1) < window.end
        ]
        el_values = [_rank_support(hit.el_rank) for hit in overlapping if hit.el_rank is not None]
        ba_values = [_rank_support(hit.ba_rank) for hit in overlapping if hit.ba_rank is not None]
        accessibility = None
        if accessibility_by_residue is not None:
            values = accessibility_by_residue[window.start : window.end]
            if len(values):
                accessibility = float(sum(values) / len(values))
        cleavage = cleavage_by_window.get(window.sequence) if cleavage_by_window else None
        output.append(
            ProcessingEvidence(
                start=window.start,
                end=window.end,
                sequence=window.sequence,
                el_presentation_support=max(el_values) if el_values else None,
                ba_binding_support=max(ba_values) if ba_values else None,
                cleavage_model_score=cleavage,
                structure_accessibility=accessibility,
                caveat=(
                    "EL supports presentation but is not a direct cleavage measurement; "
                    "BA measures binding. Missing cleavage/accessibility values remain null."
                ),
            )
        )
    return output


def tolerance_evidence(
    windows: Sequence[Window],
    self_proteome: pd.DataFrame,
    query_hla_by_window: Mapping[str, Sequence[str]],
    shared_hla_by_window: Mapping[str, Sequence[str]] | None = None,
) -> list[ToleranceEvidence]:
    """Score TCR-face self similarity and explicitly represent HLA gating state."""
    required = {"seq", "self_tcr_matches", "self_exact_9mer", "tcr_log2_enrichment"}
    missing = required - set(self_proteome.columns)
    if missing:
        raise ValueError(f"self-proteome table is missing columns: {sorted(missing)}")
    lookup = self_proteome.drop_duplicates("seq").set_index("seq")
    output = []
    for window in windows:
        row = lookup.loc[window.sequence] if window.sequence in lookup.index else None
        exact = int(row["self_exact_9mer"]) if row is not None else 0
        tcr_matches = int(row["self_tcr_matches"]) if row is not None else 0
        enrichment = float(row["tcr_log2_enrichment"]) if row is not None else float("-inf")
        query_alleles = list(query_hla_by_window.get(window.sequence, ()))
        if shared_hla_by_window is not None:
            shared = sorted(set(query_alleles) & set(shared_hla_by_window.get(window.sequence, ())))
            gate = "shared"
        elif query_alleles:
            shared = []
            gate = "query_only"
        else:
            shared = []
            gate = "not_available"

        # Bounded display score. It is eligible for fusion only when gate == "shared".
        score = 0.0 if not math.isfinite(enrichment) else 1.0 / (1.0 + math.exp(-enrichment))
        output.append(
            ToleranceEvidence(
                start=window.start,
                end=window.end,
                sequence=window.sequence,
                exact_self_9mer_count=exact,
                tcr_face_match_count=tcr_matches,
                shared_hla_alleles=shared,
                hla_gate=gate,
                tolerance_score=score,
                caveat=(
                    "TCR-face self similarity is treated as tolerance-supporting evidence only "
                    "when query and matched-self peptides share predicted HLA presentation."
                ),
            )
        )
    return output
