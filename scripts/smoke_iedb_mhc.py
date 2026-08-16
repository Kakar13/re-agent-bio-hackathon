#!/usr/bin/env python3
"""Smoke-test real NetMHCIIpan EL+BA and placeholder-response gating."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from re_agent.immuno.adapters import UnavailableResponseAdapter  # noqa: E402
from re_agent.immuno.iedb_mhc import IEDBNetMHCIIpanProvider  # noqa: E402
from re_agent.immuno.pipeline import ImmunogenicityScreeningAgent  # noqa: E402
from re_agent.immuno.registry import MHCProviderRegistry, ResponseModelRegistry  # noqa: E402

CONTROL = "ARFTGIKTAARFTGI"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Require an existing cache entry and disable HTTP calls",
    )
    args = parser.parse_args()

    endpoint = (
        "http://127.0.0.1:9/offline"
        if args.offline
        else "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
    )
    mhc = IEDBNetMHCIIpanProvider(
        ROOT / "data" / "processed" / "mhc_cache" / "iedb",
        endpoint=endpoint,
    )
    agent = ImmunogenicityScreeningAgent(
        response_registry=ResponseModelRegistry(
            [
                UnavailableResponseAdapter(
                    "teammate-model",
                    "smoke test intentionally exercises the placeholder handoff",
                )
            ]
        ),
        mhc_registry=MHCProviderRegistry([mhc]),
        default_response_adapter="teammate-model",
        hla_panel_path=ROOT / "docs" / "hla_class_ii_panel.v1.json",
        fusion_rule_path=ROOT / "docs" / "fusion_rule.v1.json",
        self_proteome_path=ROOT / "data" / "processed" / "self_proteome.parquet",
    )
    result = agent.assess("iedb-smoke-control", CONTROL)
    mhc_result = result.mhc_results[0]
    assert mhc_result.status == "ok", mhc_result.error
    assert len(mhc_result.supported_alleles) == 18
    assert not mhc_result.missing_alleles
    assert all(hit.el_rank is not None and hit.ba_rank is not None for hit in mhc_result.hits)
    assert result.response_results[0].status == "unavailable"
    assert result.combined_rank_score is None
    if args.offline:
        assert mhc_result.provenance.cached
    print(
        json.dumps(
            {
                "status": "ok",
                "alleles": len(mhc_result.supported_alleles),
                "hits": len(mhc_result.hits),
                "cached": mhc_result.provenance.cached,
                "combined_rank_score": result.combined_rank_score,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
