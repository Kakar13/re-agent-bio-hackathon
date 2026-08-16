"""Transparent late fusion over calibrated, independently inspectable proxies."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from re_agent.immuno.contracts import (
    MHCProviderResult,
    ProcessingEvidence,
    ResponseModelResult,
    ToleranceEvidence,
)


def _max_or_none(values: list[float]) -> float | None:
    return max(values) if values else None


def _rank_support(rank: float) -> float:
    return 1.0 - min(max(rank, 0.0), 100.0) / 100.0


def summarize_components(
    response: ResponseModelResult,
    mhc_results: list[MHCProviderResult],
    processing: list[ProcessingEvidence],
    tolerance: list[ToleranceEvidence],
) -> dict[str, float | None]:
    netmhc = next(
        (result for result in mhc_results if result.provider_id == "netmhciipan"),
        None,
    )
    el = (
        _max_or_none([_rank_support(hit.el_rank) for hit in netmhc.hits if hit.el_rank is not None])
        if netmhc and netmhc.status == "ok"
        else None
    )
    ba = (
        _max_or_none([_rank_support(hit.ba_rank) for hit in netmhc.hits if hit.ba_rank is not None])
        if netmhc and netmhc.status == "ok"
        else None
    )
    cleavage = _max_or_none(
        [
            value
            for row in processing
            for value in (
                row.cleavage_model_score,
                row.structure_accessibility,
            )
            if value is not None
        ]
    )
    gated = [row.tolerance_score for row in tolerance if row.hla_gate == "shared"]
    return {
        "response": _max_or_none(
            [prediction.score for prediction in response.predictions]
        )
        if response.score_scale == "calibrated_probability"
        else None,
        "netmhciipan_el": el,
        "netmhciipan_ba": ba,
        "cleavage_accessibility": cleavage,
        "inverse_hla_gated_tolerance": 1.0 - max(gated) if gated else None,
    }


def apply_fusion_rule(
    components: dict[str, float | None],
    rule_path: Path,
) -> tuple[float | None, dict[str, Any]]:
    rule = json.loads(rule_path.read_text())
    required = set(rule["required_components"])
    missing_required = sorted(key for key in required if components.get(key) is None)
    available = {
        key: float(value)
        for key, value in components.items()
        if value is not None and key in rule["weights"]
    }
    available_weight = sum(float(rule["weights"][key]) for key in available)
    detail: dict[str, Any] = {
        "rule_id": rule["rule_id"],
        "components": components,
        "missing_required": missing_required,
        "available_weight": available_weight,
        "component_contributions": {},
    }
    if missing_required or not available_weight:
        return None, detail

    score = 0.0
    for key, value in available.items():
        normalized_weight = float(rule["weights"][key]) / available_weight
        contribution = normalized_weight * value
        score += contribution
        detail["component_contributions"][key] = {
            "value": value,
            "normalized_weight": normalized_weight,
            "contribution": contribution,
        }
    detail["ablations"] = {
        "response_only": components.get("response"),
        "response_plus_el": _weighted_subset(
            components, rule["weights"], ["response", "netmhciipan_el"]
        ),
        "response_plus_el_ba": _weighted_subset(
            components,
            rule["weights"],
            ["response", "netmhciipan_el", "netmhciipan_ba"],
        ),
        "all_available": score,
    }
    return score, detail


def _weighted_subset(
    components: dict[str, float | None],
    weights: dict[str, float],
    keys: list[str],
) -> float | None:
    if any(components.get(key) is None for key in keys):
        return None
    total_weight = sum(float(weights[key]) for key in keys)
    return sum(
        float(components[key]) * float(weights[key]) / total_weight
        for key in keys
        if components[key] is not None
    )
