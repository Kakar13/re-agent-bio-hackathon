"""RFdiffusion3 -> ProteinMPNN -> AlphaFold2 -> immunogenicity handoff."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Protocol

from pydantic import BaseModel, Field

ROOT = Path(__file__).resolve().parents[3]


class ToolExecution(BaseModel):
    tool_key: str
    inputs: dict[str, Any]
    config: dict[str, Any]
    output_path: str | None = None
    status: str = "planned"


class CandidateRecord(BaseModel):
    candidate_id: str
    backbone_id: str
    sequence: str
    proteinmpnn_perplexity: float | None = None
    designed_structure_path: str
    refolded_structure_path: str | None = None
    complex_structure_path: str | None = None
    validation_metrics: dict[str, Any] = Field(default_factory=dict)
    validation_checks: list[dict[str, Any]] = Field(default_factory=list)
    validation_status: str = "pending"
    screening_status: str = "blocked"
    assessment_artifact_id: str | None = None


class CampaignManifest(BaseModel):
    schema_version: str = "1.0.0"
    campaign_id: str
    spec_path: str
    status: str
    target: dict[str, Any]
    claim_boundary: str
    executions: list[ToolExecution]
    candidates: list[CandidateRecord] = Field(default_factory=list)
    candidate_counts: dict[str, int] = Field(default_factory=dict)
    screening_fasta: str | None = None
    warnings: list[str] = Field(default_factory=list)


class ProtoRunner(Protocol):
    def run(
        self,
        tool_key: str,
        inputs: Mapping[str, Any],
        config: Mapping[str, Any],
    ) -> dict[str, Any]: ...


class ProtoLocalRunner:
    """Execute exact Proto tool keys through their typed Python APIs."""

    def run(
        self,
        tool_key: str,
        inputs: Mapping[str, Any],
        config: Mapping[str, Any],
    ) -> dict[str, Any]:
        if tool_key == "rfdiffusion3-design":
            from proto_tools.tools.structure_design.rfdiffusion3.rfdiffusion3_sample import (
                RFdiffusion3Config,
                RFdiffusion3Input,
                run_rfdiffusion3,
            )

            result = run_rfdiffusion3(
                RFdiffusion3Input.model_validate(inputs),
                RFdiffusion3Config.model_validate(config),
            )
        elif tool_key == "proteinmpnn-sample":
            from proto_tools.tools.inverse_folding.proteinmpnn.proteinmpnn_sample import (
                ProteinMPNNSampleConfig,
                run_proteinmpnn_sample,
            )
            from proto_tools.tools.inverse_folding.shared_data_models import (
                InverseFoldingInput,
            )

            result = run_proteinmpnn_sample(
                InverseFoldingInput.model_validate(inputs),
                ProteinMPNNSampleConfig.model_validate(config),
            )
        elif tool_key == "alphafold2-prediction":
            from proto_tools.tools.structure_prediction.alphafold2.alphafold2 import (
                AlphaFold2Config,
                AlphaFold2Input,
                run_alphafold2,
            )

            result = run_alphafold2(
                AlphaFold2Input.model_validate(inputs),
                AlphaFold2Config.model_validate(config),
            )
        else:
            raise ValueError(f"unsupported campaign tool: {tool_key}")
        return result.model_dump(mode="json")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _structure_text(payload: Any) -> str:
    value = payload.get("structure") if isinstance(payload, Mapping) else payload
    if not isinstance(value, str):
        raise ValueError("Proto structure output did not contain PDB text or a file path")
    if "\n" in value or value.startswith(("ATOM", "HEADER", "data_")):
        return value
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = ROOT / path
    if not path.is_file():
        raise FileNotFoundError(f"Proto structure artifact is missing: {path}")
    return path.read_text()


def _interchain_clash_count(
    path: Path,
    *,
    binder_chain: str = "A",
    cutoff_angstrom: float = 2.0,
) -> int:
    binder_atoms: list[tuple[float, float, float]] = []
    context_atoms: list[tuple[float, float, float]] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        element = line[76:78].strip() if len(line) >= 78 else line[12:14].strip()
        if element.upper().startswith("H"):
            continue
        try:
            coordinates = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue
        if line[21].strip() == binder_chain:
            binder_atoms.append(coordinates)
        else:
            context_atoms.append(coordinates)
    cutoff_squared = cutoff_angstrom**2
    return sum(
        1
        for x1, y1, z1 in binder_atoms
        for x2, y2, z2 in context_atoms
        if (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2 < cutoff_squared
    )


def _interface_pae_from_matrix(
    metrics: Mapping[str, Any],
    chains: list[Mapping[str, Any]],
    designed: list[bool],
) -> float | None:
    matrix = (
        metrics.get("pae")
        or metrics.get("pae_matrix")
        or metrics.get("predicted_aligned_error")
    )
    if not isinstance(matrix, list) or not matrix or not isinstance(matrix[0], list):
        return None
    binder_indices: list[int] = []
    context_indices: list[int] = []
    offset = 0
    for chain, redesigned in zip(chains, designed, strict=False):
        length = len(chain.get("sequence", ""))
        destination = binder_indices if redesigned else context_indices
        destination.extend(range(offset, offset + length))
        offset += length
    values: list[float] = []
    for binder_index in binder_indices:
        for context_index in context_indices:
            try:
                values.append(float(matrix[binder_index][context_index]))
                values.append(float(matrix[context_index][binder_index]))
            except (IndexError, TypeError, ValueError):
                continue
    return sum(values) / len(values) if values else None


def _nested_value(payload: Mapping[str, Any], dotted_paths: tuple[str, ...]) -> Any:
    for dotted_path in dotted_paths:
        value: Any = payload
        for part in dotted_path.split("."):
            if not isinstance(value, Mapping) or part not in value:
                value = None
                break
            value = value[part]
        if value is not None:
            return value
    return None


def _numeric_metric(payload: Mapping[str, Any], dotted_paths: tuple[str, ...]) -> float | None:
    value = _nested_value(payload, dotted_paths)
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, list):
        numeric = [float(item) for item in value if isinstance(item, (int, float))]
        return sum(numeric) / len(numeric) if numeric else None
    return None


def evaluate_candidate_validation(
    candidates: list[CandidateRecord],
    validation: Mapping[str, Any],
) -> list[CandidateRecord]:
    """Apply the frozen structural gates and block incomplete candidates."""

    thresholds = validation.get("literature_prior_thresholds", {})
    sequence_counts = Counter(row.sequence for row in candidates)
    plddt_min = float(thresholds.get("monomer_plddt_min", 80.0))
    interface_pae_max = float(thresholds.get("interface_pae_max", 10.0))
    binder_rmsd_max = float(thresholds.get("binder_rmsd_angstrom_max", 1.0))

    for candidate in candidates:
        metrics = candidate.validation_metrics
        plddt = _numeric_metric(
            metrics,
            (
                "monomer_plddt",
                "mean_plddt",
                "avg_plddt",
                "plddt",
                "confidence.plddt",
                "confidence.mean_plddt",
            ),
        )
        interface_pae = _numeric_metric(
            metrics,
            (
                "interface_pae",
                "mean_interface_pae",
                "ipae",
                "confidence.interface_pae",
            ),
        )
        binder_rmsd = _numeric_metric(
            metrics,
            (
                "binder_rmsd",
                "binder_rmsd_angstrom",
                "designed_vs_refolded_binder_rmsd",
                "rmsd.binder",
            ),
        )
        clash_count = _numeric_metric(
            metrics,
            ("clash_count", "num_clashes", "clashes.count"),
        )
        clash_passed = _nested_value(
            metrics,
            ("clash_check_passed", "clashes.passed"),
        )
        if clash_passed is None and clash_count is not None:
            clash_passed = clash_count == 0

        checks = [
            {
                "name": "unique_sequence",
                "recorded": True,
                "value": sequence_counts[candidate.sequence],
                "passed": sequence_counts[candidate.sequence] == 1,
            },
            {
                "name": "proteinmpnn_perplexity_recorded",
                "recorded": candidate.proteinmpnn_perplexity is not None,
                "value": candidate.proteinmpnn_perplexity,
                "passed": candidate.proteinmpnn_perplexity is not None,
            },
            {
                "name": "monomer_structure_confidence_recorded",
                "recorded": plddt is not None,
                "value": plddt,
                "threshold": {"operator": ">=", "value": plddt_min},
                "passed": plddt is not None and plddt >= plddt_min,
            },
            {
                "name": "designed_vs_refolded_binder_rmsd_recorded",
                "recorded": binder_rmsd is not None,
                "value": binder_rmsd,
                "threshold": {"operator": "<=", "value": binder_rmsd_max},
                "passed": binder_rmsd is not None and binder_rmsd <= binder_rmsd_max,
            },
            {
                "name": "interface_confidence_recorded",
                "recorded": interface_pae is not None,
                "value": interface_pae,
                "threshold": {"operator": "<=", "value": interface_pae_max},
                "passed": interface_pae is not None and interface_pae <= interface_pae_max,
            },
            {
                "name": "clash_check_recorded",
                "recorded": isinstance(clash_passed, bool),
                "value": clash_count if clash_count is not None else clash_passed,
                "passed": clash_passed is True,
            },
        ]
        candidate.validation_checks = checks
        candidate.validation_status = (
            "pass" if all(bool(check["passed"]) for check in checks) else "fail"
        )
        candidate.screening_status = (
            "eligible" if candidate.validation_status == "pass" else "blocked"
        )
    return candidates


def build_campaign_plan(spec_path: Path, *, device: str = "cuda") -> CampaignManifest:
    spec = json.loads(spec_path.read_text())
    target_path = ROOT / spec["target"]["structure_path"]
    if not target_path.exists():
        raise FileNotFoundError(f"pinned target structure is missing: {target_path}")
    expected_hash = spec["target"]["structure_sha256"]
    actual_hash = _sha256(target_path)
    if actual_hash != expected_hash:
        raise ValueError(
            f"target structure hash mismatch: expected {expected_hash}, got {actual_hash}"
        )

    backbone = spec["generation"]["backbone_provider"]
    sequence = spec["generation"]["sequence_provider"]
    executions = [
        ToolExecution(
            tool_key=backbone["tool_key"],
            inputs={
                "design_specs": [
                    {
                        "input_structure": str(target_path),
                        "length": backbone["binder_length"],
                        "select_hotspots": backbone["select_hotspots"],
                        "infer_ori_strategy": "hotspots",
                    }
                ]
            },
            config={
                "n_batches": 1,
                "diffusion_batch_size": backbone["n_backbones"],
                "seed": backbone["seed"],
                "device": device,
            },
        ),
        ToolExecution(
            tool_key=sequence["tool_key"],
            inputs={"inputs": "<one designed structure per execution>"},
            config={
                "num_sequences_per_structure": sequence["sequences_per_backbone"],
                "temperature": sequence["temperature"],
                "model_choice": sequence["model_choice"],
                "seed": sequence["seed"],
                "device": device,
            },
        ),
        ToolExecution(
            tool_key="alphafold2-prediction",
            inputs={
                "complexes": [
                    "<one binder monomer per candidate>",
                    "<one binder-target complex per candidate>",
                ]
            },
            config={
                "use_msa": False,
                "num_recycles": 3,
                "include_pae_matrix": True,
                "device": device,
            },
        ),
        ToolExecution(
            tool_key="pymol-rmsd-alignment",
            inputs={
                "target_structure": "<designed binder complex>",
                "mobile_structure": "<AlphaFold2-refolded complex>",
            },
            config={
                "method": "align",
                "target_selection": "target and chain A and name CA",
                "mobile_selection": "mobile and chain A and name CA",
                "include_superimposed_pdb": False,
                "device": "cpu",
            },
        ),
    ]
    return CampaignManifest(
        campaign_id=spec["spec_id"],
        spec_path=_relative(spec_path),
        status="planned",
        target=spec["target"],
        claim_boundary=spec["objective"]["claim_boundary"],
        executions=executions,
        warnings=[
            "AlphaFold2 confidence is a computational validation proxy, not binding proof.",
            "Candidate screening starts only after required structural metrics are recorded.",
        ],
    )


def run_campaign(
    spec_path: Path,
    output_dir: Path,
    runner: ProtoRunner,
    *,
    approved: bool,
    device: str = "cuda",
) -> CampaignManifest:
    """Run a bounded campaign; callers must pass an explicit compute approval."""

    manifest = build_campaign_plan(spec_path, device=device)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "campaign_manifest.json"
    if not approved:
        manifest.status = "awaiting_compute_approval"
        manifest_path.write_text(manifest.model_dump_json(indent=2) + "\n")
        return manifest

    spec = json.loads(spec_path.read_text())
    backbone_execution = manifest.executions[0]
    backbone_output = runner.run(
        backbone_execution.tool_key,
        backbone_execution.inputs,
        backbone_execution.config,
    )
    backbone_path = output_dir / "rfdiffusion3_output.json"
    backbone_path.write_text(json.dumps(backbone_output, indent=2) + "\n")
    backbone_execution.output_path = _relative(backbone_path)
    backbone_execution.status = "completed"

    designed = [
        structure
        for bundle in backbone_output.get("designed_structures", [])
        for structure in bundle.get("structures", [])
    ]
    candidates: list[CandidateRecord] = []
    sequence_config = manifest.executions[1].config
    af2_config = manifest.executions[2].config
    rmsd_config = manifest.executions[3].config

    for backbone_index, designed_row in enumerate(designed):
        backbone_id = f"bb-{backbone_index:03d}"
        structure_payload = designed_row["structure"]
        designed_path = output_dir / f"{backbone_id}.pdb"
        designed_path.write_text(_structure_text(structure_payload))
        mpnn_structure_payload = (
            {**structure_payload, "structure": str(designed_path)}
            if isinstance(structure_payload, Mapping)
            else str(designed_path)
        )
        mpnn_output = runner.run(
            "proteinmpnn-sample",
            {
                "inputs": [
                    {
                        "structure": mpnn_structure_payload,
                        "chains_to_redesign": spec["generation"]["sequence_provider"][
                            "chains_to_redesign"
                        ],
                    }
                ]
            },
            sequence_config,
        )
        mpnn_path = output_dir / f"{backbone_id}.proteinmpnn.json"
        mpnn_path.write_text(json.dumps(mpnn_output, indent=2) + "\n")
        complexes = [
            row
            for design_set in mpnn_output.get("design_sets", [])
            for row in design_set.get("complexes", [])
        ]
        for sequence_index, complex_row in enumerate(complexes):
            candidate_id = f"{backbone_id}-seq-{sequence_index:02d}"
            designed_chains = [
                chain
                for chain, redesigned in zip(
                    complex_row["chains"],
                    complex_row.get("designed", []),
                    strict=False,
                )
                if redesigned and chain.get("entity_type", "protein") == "protein"
            ]
            if not designed_chains:
                continue
            candidate_sequence = "".join(chain["sequence"] for chain in designed_chains)
            monomer_output = runner.run(
                "alphafold2-prediction",
                {"complexes": [{"chains": designed_chains}]},
                af2_config,
            )
            monomer_path = output_dir / f"{candidate_id}.alphafold2.monomer.json"
            monomer_path.write_text(json.dumps(monomer_output, indent=2) + "\n")
            complex_output = runner.run(
                "alphafold2-prediction",
                {"complexes": [{"chains": complex_row["chains"]}]},
                af2_config,
            )
            complex_path = output_dir / f"{candidate_id}.alphafold2.complex.json"
            complex_path.write_text(json.dumps(complex_output, indent=2) + "\n")
            monomer_structures = monomer_output.get("structures", [])
            complex_structures = complex_output.get("structures", [])
            refold_path: Path | None = None
            refolded_complex_path: Path | None = None
            metrics: dict[str, Any] = {}
            if monomer_structures:
                refold_path = output_dir / f"{candidate_id}.pdb"
                refold_path.write_text(_structure_text(monomer_structures[0]))
                monomer_metrics = monomer_structures[0].get("metrics", {})
                if isinstance(monomer_metrics, dict):
                    monomer_plddt = _numeric_metric(
                        monomer_metrics,
                        ("monomer_plddt", "mean_plddt", "avg_plddt", "plddt"),
                    )
                    if monomer_plddt is not None:
                        metrics["monomer_plddt"] = monomer_plddt
                    metrics["monomer_metrics"] = monomer_metrics
            if complex_structures:
                refolded_complex_path = output_dir / f"{candidate_id}.complex.pdb"
                refolded_complex_path.write_text(_structure_text(complex_structures[0]))
                complex_metrics = complex_structures[0].get("metrics", {})
                if not isinstance(complex_metrics, dict):
                    complex_metrics = {}
                metrics["complex_metrics"] = complex_metrics
                interface_pae = _interface_pae_from_matrix(
                    complex_metrics,
                    complex_row["chains"],
                    complex_row.get("designed", []),
                )
                if interface_pae is not None:
                    metrics["interface_pae"] = interface_pae
                clash_count = _interchain_clash_count(refolded_complex_path)
                metrics["clash_count"] = clash_count
                metrics["clash_check_passed"] = clash_count == 0
            if refold_path:
                try:
                    rmsd_output = runner.run(
                        "pymol-rmsd-alignment",
                        {
                            "target_structure": str(designed_path),
                            "mobile_structure": str(refold_path),
                        },
                        rmsd_config,
                    )
                    rmsd = rmsd_output.get("metrics", {}).get("rmsd")
                    if isinstance(rmsd, (int, float)):
                        metrics["binder_rmsd"] = float(rmsd)
                except Exception as exc:
                    metrics["binder_rmsd_error"] = f"{type(exc).__name__}: {exc}"
            candidates.append(
                CandidateRecord(
                    candidate_id=candidate_id,
                    backbone_id=backbone_id,
                    sequence=candidate_sequence,
                    proteinmpnn_perplexity=complex_row.get("metrics", {}).get("perplexity"),
                    designed_structure_path=_relative(designed_path),
                    refolded_structure_path=_relative(refold_path) if refold_path else None,
                    complex_structure_path=(
                        _relative(refolded_complex_path) if refolded_complex_path else None
                    ),
                    validation_metrics=metrics,
                )
            )

    fasta_path = output_dir / "screening_candidates.fasta"
    fasta_path.write_text(
        "".join(f">{row.candidate_id}\n{row.sequence}\n" for row in candidates)
    )
    manifest.candidates = candidates
    evaluate_candidate_validation(candidates, spec["validation"])
    eligible = [row for row in candidates if row.screening_status == "eligible"]
    fasta_path.write_text(
        "".join(f">{row.candidate_id}\n{row.sequence}\n" for row in eligible)
    )
    manifest.candidate_counts = {
        "generated": len(candidates),
        "validation_passed": len(eligible),
        "validation_failed": len(candidates) - len(eligible),
    }
    manifest.screening_fasta = _relative(fasta_path)
    if not candidates:
        manifest.status = "no_candidates"
    elif eligible:
        manifest.status = "ready_for_screening"
    else:
        manifest.status = "validation_blocked"
    manifest.executions[1].status = "completed"
    manifest.executions[2].status = "completed"
    manifest.executions[3].status = "completed"
    manifest_path.write_text(manifest.model_dump_json(indent=2) + "\n")
    return manifest
