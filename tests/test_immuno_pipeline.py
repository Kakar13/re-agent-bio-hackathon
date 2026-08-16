from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from re_agent.immuno.adapters import sequence_sha256
from re_agent.immuno.contracts import (
    MHCHit,
    MHCProviderResult,
    Provenance,
    ResponseModelResult,
    ResponsePrediction,
)
from re_agent.immuno.pipeline import ImmunogenicityScreeningAgent
from re_agent.immuno.proxies import tile_sequence
from re_agent.immuno.registry import MHCProviderRegistry, ResponseModelRegistry

SEQUENCE = "ACDEFGHIKLMNPQRSTVWY"


class ResponseAdapter:
    adapter_id = "teammate-model"
    version = "checkpoint-1"

    def predict(self, request):
        return ResponseModelResult(
            adapter_id=self.adapter_id,
            status="ok",
            score_scale="calibrated_probability",
            calibration_id="unit-test-calibration",
            predictions=[
                ResponsePrediction(
                    start=window.start,
                    end=window.end,
                    sequence=window.sequence,
                    score=0.7,
                )
                for window in request.windows
            ],
            provenance=Provenance(
                provider=self.adapter_id,
                version=self.version,
                capability="response_model",
                source="unit-test",
                input_sha256=sequence_sha256(request.parent_sequence),
            ),
        )


class MHCAdapter:
    provider_id = "netmhciipan"
    version = "4.3"

    def predict(self, request):
        return MHCProviderResult(
            provider_id=self.provider_id,
            status="ok",
            hits=[
                MHCHit(
                    allele=request.alleles[0],
                    start=1,
                    end=15,
                    peptide=request.parent_sequence[:15],
                    core=request.parent_sequence[:9],
                    el_rank=1.5,
                    ba_rank=3.0,
                    ba_ic50_nm=42.0,
                )
            ],
            supported_alleles=[request.alleles[0]],
            missing_alleles=request.alleles[1:],
            provenance=Provenance(
                provider=self.provider_id,
                version=self.version,
                capability="mhc_evidence",
                source="unit-test",
                input_sha256=sequence_sha256(request.parent_sequence),
            ),
        )


class PipelineTests(unittest.TestCase):
    def test_tile_coordinates_match_sequence(self):
        windows = tile_sequence(SEQUENCE)
        self.assertEqual(len(windows), 6)
        for window in windows:
            self.assertEqual(SEQUENCE[window.start : window.end], window.sequence)

    def test_pipeline_keeps_el_ba_and_adapter_output_separate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            panel = root / "panel.json"
            panel.write_text(
                json.dumps(
                    {
                        "panel_id": "test-panel",
                        "alleles": ["DRB1*01:01", "DRB1*15:01"],
                    }
                )
            )
            rule = root / "rule.json"
            rule.write_text(
                json.dumps(
                    {
                        "rule_id": "test-rule",
                        "required_components": [
                            "response",
                            "netmhciipan_el",
                            "netmhciipan_ba",
                        ],
                        "weights": {
                            "response": 0.5,
                            "netmhciipan_el": 0.3,
                            "netmhciipan_ba": 0.2,
                        },
                    }
                )
            )
            self_table = root / "self.parquet"
            pd.DataFrame(
                {
                    "seq": [window.sequence for window in tile_sequence(SEQUENCE)],
                    "self_tcr_matches": [2] * 6,
                    "self_exact_9mer": [0] * 6,
                    "tcr_log2_enrichment": [0.25] * 6,
                }
            ).to_parquet(self_table, index=False)

            agent = ImmunogenicityScreeningAgent(
                ResponseModelRegistry([ResponseAdapter()]),
                MHCProviderRegistry([MHCAdapter()]),
                "teammate-model",
                panel,
                rule,
                self_table,
            )
            result = agent.assess("candidate-1", SEQUENCE)

        self.assertIsNotNone(result.combined_rank_score)
        values = result.component_summary["proxy_values"]
        self.assertNotEqual(values["netmhciipan_el"], values["netmhciipan_ba"])
        ablations = result.component_summary["fusion"]["ablations"]
        self.assertIsNotNone(ablations["response_only"])
        self.assertIsNotNone(ablations["response_plus_el_ba"])
        self.assertEqual(result.response_results[0].adapter_id, "teammate-model")
        self.assertTrue(all(row.hla_gate == "query_only" for row in result.tolerance))
        self.assertTrue(
            any("excluded from fusion" in warning for warning in result.warnings)
        )


if __name__ == "__main__":
    unittest.main()
