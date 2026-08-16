"""Project overlapping evidence windows onto inspectable per-residue tracks."""

from __future__ import annotations

from collections.abc import Iterable

from re_agent.immuno.contracts import (
    MHCProviderResult,
    ProcessingEvidence,
    ResponseModelResult,
    ToleranceEvidence,
)


def _scatter(
    length: int,
    spans: Iterable[tuple[int, int, float]],
) -> list[float]:
    totals = [0.0] * length
    counts = [0] * length
    for start, end, value in spans:
        for index in range(max(0, start), min(length, end)):
            totals[index] += value
            counts[index] += 1
    return [
        round(total / count, 6) if count else 0.0
        for total, count in zip(totals, counts, strict=True)
    ]


def build_spatial_tracks(
    sequence_length: int,
    response: ResponseModelResult | None,
    mhc_results: list[MHCProviderResult],
    processing: list[ProcessingEvidence],
    tolerance: list[ToleranceEvidence],
) -> dict[str, list[float]]:
    response_spans = (
        [
            (prediction.start, prediction.end, prediction.score)
            for prediction in response.predictions
        ]
        if response is not None
        else []
    )
    netmhc = next(
        (result for result in mhc_results if result.provider_id == "netmhciipan"),
        None,
    )
    el_spans = []
    ba_spans = []
    if netmhc is not None:
        for hit in netmhc.hits:
            if hit.el_rank is not None:
                el_spans.append((hit.start - 1, hit.end, 1.0 - min(hit.el_rank, 100.0) / 100.0))
            if hit.ba_rank is not None:
                ba_spans.append((hit.start - 1, hit.end, 1.0 - min(hit.ba_rank, 100.0) / 100.0))
    accessibility_spans = [
        (row.start, row.end, row.structure_accessibility)
        for row in processing
        if row.structure_accessibility is not None
    ]
    tolerance_spans = [
        (row.start, row.end, row.tolerance_score)
        for row in tolerance
        if row.hla_gate == "shared"
    ]
    return {
        "response_propensity": _scatter(sequence_length, response_spans),
        "netmhciipan_el_support": _scatter(sequence_length, el_spans),
        "netmhciipan_ba_support": _scatter(sequence_length, ba_spans),
        "structure_accessibility": _scatter(sequence_length, accessibility_spans),
        "hla_gated_tolerance_support": _scatter(sequence_length, tolerance_spans),
    }
