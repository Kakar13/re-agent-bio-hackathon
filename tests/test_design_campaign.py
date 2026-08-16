from __future__ import annotations

import json
import shutil
from pathlib import Path

from re_agent.design.campaign import (
    CandidateRecord,
    build_campaign_plan,
    evaluate_candidate_validation,
    run_campaign,
)

ROOT = Path(__file__).resolve().parents[1]


class FailIfCalledRunner:
    def run(self, tool_key, inputs, config):
        raise AssertionError(f"compute should not run: {tool_key}")


def _pdb(chains: list[tuple[str, str]], *, chain_offset: float = 200.0) -> str:
    lines = []
    serial = 1
    three_letter = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR",
    }
    for chain_index, (chain_id, sequence) in enumerate(chains):
        for residue_number, residue in enumerate(sequence, start=1):
            x = chain_index * chain_offset + residue_number * 3.8
            lines.append(
                f"ATOM  {serial:5d}  CA  {three_letter[residue]:>3} {chain_id}"
                f"{residue_number:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}"
                "  1.00 90.00           C"
            )
            serial += 1
    return "\n".join(lines) + "\nEND\n"


class PassingCampaignRunner:
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    def run(self, tool_key, inputs, config):
        if tool_key == "rfdiffusion3-design":
            return {
                "designed_structures": [
                    {
                        "structures": [
                            {
                                "structure": {
                                    "structure": _pdb(
                                        [("A", self.sequence), ("B", "A")],
                                    )
                                }
                            }
                        ]
                    }
                ]
            }
        if tool_key == "proteinmpnn-sample":
            return {
                "design_sets": [
                    {
                        "complexes": [
                            {
                                "chains": [
                                    {"id": "A", "sequence": self.sequence},
                                    {"id": "B", "sequence": "A"},
                                ],
                                "designed": [True, False],
                                "metrics": {"perplexity": 1.2},
                            }
                        ]
                    }
                ]
            }
        if tool_key == "alphafold2-prediction":
            chains = inputs["complexes"][0]["chains"]
            if len(chains) == 1:
                return {
                    "structures": [
                        {
                            "structure": _pdb([("A", self.sequence)]),
                            "metrics": {"avg_plddt": 91.0},
                        }
                    ]
                }
            size = len(self.sequence) + 1
            pae = [[5.0 for _ in range(size)] for _ in range(size)]
            return {
                "structures": [
                    {
                        "structure": _pdb([("A", self.sequence), ("B", "A")]),
                        "metrics": {"pae": pae},
                    }
                ]
            }
        if tool_key == "pymol-rmsd-alignment":
            return {"metrics": {"rmsd": 0.5}}
        raise AssertionError(tool_key)


def test_campaign_plan_matches_frozen_tool_contracts() -> None:
    plan = build_campaign_plan(ROOT / "docs" / "design_spec.json")

    assert [row.tool_key for row in plan.executions] == [
        "rfdiffusion3-design",
        "proteinmpnn-sample",
        "alphafold2-prediction",
        "pymol-rmsd-alignment",
    ]
    rfd3 = plan.executions[0]
    assert rfd3.inputs["design_specs"][0]["select_hotspots"] == "B62,B84,B143"
    assert rfd3.config["diffusion_batch_size"] == 8
    assert [execution.config["device"] for execution in plan.executions] == [
        "cuda",
        "cuda",
        "cuda",
        "cpu",
    ]


def test_campaign_requires_explicit_compute_approval(tmp_path: Path) -> None:
    manifest = run_campaign(
        ROOT / "docs" / "design_spec.json",
        tmp_path,
        FailIfCalledRunner(),
        approved=False,
    )

    assert manifest.status == "awaiting_compute_approval"
    assert (tmp_path / "campaign_manifest.json").exists()


def test_structural_validation_fails_closed_on_missing_metrics() -> None:
    candidate = CandidateRecord(
        candidate_id="candidate-1",
        backbone_id="bb-1",
        sequence="ACDEFGHIKLMNPQRSTVWY",
        designed_structure_path="results/design.pdb",
        validation_metrics={"mean_plddt": 91.0},
    )
    spec = build_campaign_plan(ROOT / "docs" / "design_spec.json")
    validation = json.loads((ROOT / "docs" / "design_spec.json").read_text())["validation"]

    evaluated = evaluate_candidate_validation([candidate], validation)

    assert spec.status == "planned"
    assert evaluated[0].validation_status == "fail"
    assert evaluated[0].screening_status == "blocked"
    assert any(
        check["name"] == "interface_confidence_recorded" and not check["recorded"]
        for check in evaluated[0].validation_checks
    )


def test_structural_validation_allows_only_complete_passing_candidate() -> None:
    candidate = CandidateRecord(
        candidate_id="candidate-1",
        backbone_id="bb-1",
        sequence="ACDEFGHIKLMNPQRSTVWY",
        proteinmpnn_perplexity=1.4,
        designed_structure_path="results/design.pdb",
        refolded_structure_path="results/refold.pdb",
        validation_metrics={
            "mean_plddt": 91.0,
            "interface_pae": 7.0,
            "binder_rmsd": 0.6,
            "clash_count": 0,
        },
    )
    validation = json.loads((ROOT / "docs" / "design_spec.json").read_text())["validation"]

    evaluated = evaluate_candidate_validation([candidate], validation)

    assert evaluated[0].validation_status == "pass"
    assert evaluated[0].screening_status == "eligible"
    assert all(check["passed"] for check in evaluated[0].validation_checks)


def test_campaign_hands_only_structural_passes_to_screening_fasta(tmp_path: Path) -> None:
    output_dir = ROOT / "results" / "test-runs" / tmp_path.name
    try:
        manifest = run_campaign(
            ROOT / "docs" / "design_spec.json",
            output_dir,
            PassingCampaignRunner(),
            approved=True,
        )

        assert manifest.status == "ready_for_screening"
        assert manifest.candidate_counts == {
            "generated": 1,
            "validation_passed": 1,
            "validation_failed": 0,
        }
        assert manifest.candidates[0].validation_status == "pass"
        screening_fasta = output_dir / "screening_candidates.fasta"
        assert PassingCampaignRunner.sequence in screening_fasta.read_text()
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)
