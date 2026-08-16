"""Inspectable scientific tools exposed through LangGraph."""

from __future__ import annotations

import asyncio
import csv
import hashlib
import json
import re
import shlex
import subprocess
import urllib.request
from pathlib import Path
from typing import TYPE_CHECKING, Any

from langchain_core.tools import tool
from langsmith import traceable

from re_agent.agent.mcp_runtime import (
    RuntimeMCPError,
    invoke_runtime_mcp,
    runtime_mcp_status,
)
from re_agent.agent.run_store import (
    ROOT,
    latest_run_artifact,
    list_runs,
    new_run_id,
    write_run_artifact,
)
from re_agent.design.campaign import ProtoLocalRunner, run_campaign
from re_agent.immuno.adapters import ArtifactResponseAdapter, UnavailableResponseAdapter
from re_agent.immuno.contracts import CandidateAssessment
from re_agent.immuno.iedb_mhc import IEDBNetMHCIIpanProvider
from re_agent.immuno.pipeline import ImmunogenicityScreeningAgent
from re_agent.immuno.registry import MHCProviderRegistry, ResponseModelRegistry
from re_agent.immuno.structure import structure_reference_from_pdb

if TYPE_CHECKING:
    from re_agent.immuno.e2e_pls_pickle import TeamE2EPLSAdapter

SKILL_ROOTS = (
    ROOT / "skills",
    ROOT / "vendor" / "scientific-agent-skills" / "skills",
)

DEFAULT_SPEC_PATH = ROOT / "docs" / "design_spec.json"
REFERENCE_FASTA = ROOT / "data" / "raw" / "IL7Ra_binders_sequences.fasta"
MAX_PDB_BYTES = 25_000_000


def _resolve_repo_file(path: str, *, allowed_roots: tuple[Path, ...]) -> Path:
    target = (ROOT / path).resolve()
    if not any(root.resolve() in target.parents for root in allowed_roots):
        raise ValueError(f"path must be below: {', '.join(str(root) for root in allowed_roots)}")
    if not target.is_file():
        raise FileNotFoundError(target)
    return target


def _resolve_repo_directory(path: str, *, allowed_roots: tuple[Path, ...]) -> Path:
    target = (ROOT / path).resolve()
    if not any(root.resolve() in target.parents for root in allowed_roots):
        raise ValueError(f"path must be below: {', '.join(str(root) for root in allowed_roots)}")
    if not target.is_dir():
        raise FileNotFoundError(target)
    return target


def _response_adapter(response_artifact: str | None):
    if response_artifact:
        response_path = _resolve_repo_file(response_artifact, allowed_roots=(ROOT,))
        return ArtifactResponseAdapter(response_path)
    return UnavailableResponseAdapter(
        "teammate-model",
        "No response artifact supplied; combined rank is intentionally withheld.",
    )


def _screening_agent(response_artifact: str | None) -> ImmunogenicityScreeningAgent:
    response_adapter = _response_adapter(response_artifact)
    mhc_provider = IEDBNetMHCIIpanProvider(
        ROOT / "data" / "processed" / "mhc_cache" / "iedb"
    )
    return ImmunogenicityScreeningAgent(
        response_registry=ResponseModelRegistry([response_adapter]),
        mhc_registry=MHCProviderRegistry([mhc_provider]),
        default_response_adapter=response_adapter.adapter_id,
        hla_panel_path=ROOT / "docs" / "hla_class_ii_panel.v1.json",
        fusion_rule_path=ROOT / "docs" / "fusion_rule.v1.json",
        self_proteome_path=ROOT / "data" / "processed" / "self_proteome.parquet",
    )


def _normalized_sequence(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper()
    if not normalized or any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in normalized):
        raise ValueError("sequence must contain canonical amino-acid letters only")
    if len(normalized) < 15:
        raise ValueError("sequence must contain at least one 15-mer")
    return normalized


def _cached_rcsb_pdb(pdb_id: str) -> Path:
    normalized = pdb_id.strip().upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", normalized):
        raise ValueError("RCSB PDB ID must contain four alphanumeric characters")
    output = ROOT / "data" / "structures" / f"{normalized.lower()}.pdb"
    if output.is_file():
        return output

    request = urllib.request.Request(
        f"https://files.rcsb.org/download/{normalized}.pdb",
        headers={"User-Agent": "re-agent-workbench/1.0"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        payload = response.read(MAX_PDB_BYTES + 1)
    if len(payload) > MAX_PDB_BYTES:
        raise ValueError(f"RCSB structure {normalized} exceeds the size limit")
    if b"ATOM  " not in payload:
        raise ValueError(f"RCSB structure {normalized} contains no atomic coordinates")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".pdb.tmp")
    temporary.write_bytes(payload)
    temporary.replace(output)
    return output


def _assess_candidate(
    agent: ImmunogenicityScreeningAgent,
    *,
    candidate_id: str,
    sequence: str,
    structure_path: str | None = None,
    structure_chain_id: str = "A",
) -> CandidateAssessment:
    normalized = _normalized_sequence(sequence)
    assessment = agent.assess(candidate_id, normalized)
    if structure_path:
        resolved_structure = _resolve_repo_file(
            structure_path,
            allowed_roots=(ROOT / "data", ROOT / "results"),
        )
        assessment.structure = structure_reference_from_pdb(
            resolved_structure,
            sequence=normalized,
            chain_id=structure_chain_id,
            repository_root=ROOT,
        )
    return assessment


def _mhci_surrogates(
    checkpoint: str | None,
    checkpoints: list[str] | None = None,
    netmhcpan_checkpoint: str | None = None,
) -> list[TeamE2EPLSAdapter]:
    requested = [path for path in [checkpoint, *(checkpoints or [])] if path]
    if not requested:
        if netmhcpan_checkpoint:
            raise ValueError(
                "mhci_netmhcpan_checkpoint also requires a legacy MHC-I checkpoint "
                "for the cleavage and TAP factors"
            )
        return []
    try:
        from re_agent.immuno.e2e_pls_pickle import TeamE2EPLSAdapter
    except ImportError as exc:
        raise RuntimeError(
            "team MHC-I surrogate requires the ML dependencies: uv sync --extra ml"
        ) from exc
    unique_paths = list(dict.fromkeys(requested))
    netmhcpan_checkpoint_dir = (
        _resolve_repo_directory(netmhcpan_checkpoint, allowed_roots=(ROOT / "models",))
        if netmhcpan_checkpoint
        else None
    )
    return [
        TeamE2EPLSAdapter(
            _resolve_repo_file(path, allowed_roots=(ROOT / "models",)),
            netmhcpan_checkpoint_dir=netmhcpan_checkpoint_dir,
        )
        for path in unique_paths
    ]


def _attach_mhci_surrogates(
    assessment: CandidateAssessment,
    adapters: list[TeamE2EPLSAdapter],
) -> None:
    for adapter in adapters:
        result = adapter.predict(assessment.sequence)
        assessment.mhc_i_surrogate_results.append(result)
        for track_name, values in result.spatial_tracks.items():
            output_name = (
                f"{result.adapter_id}_{track_name}" if len(adapters) > 1 else track_name
            )
            assessment.spatial_tracks[output_name] = values


def _path_status(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(ROOT)),
        "available": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() and path.is_file() else None,
    }


def _json_safe(value: Any) -> Any:
    return json.loads(json.dumps(value, default=str))


def _unwrap_mcp_payload(value: Any) -> Any:
    if isinstance(value, list) and len(value) == 1:
        return _unwrap_mcp_payload(value[0])
    if isinstance(value, dict) and value.get("type") == "text" and "text" in value:
        text = value["text"]
        try:
            return _unwrap_mcp_payload(json.loads(text))
        except (TypeError, json.JSONDecodeError):
            return text
    if isinstance(value, dict) and "result" in value and isinstance(value["result"], dict):
        return value["result"]
    return value


class ProtoMCPRunner:
    """Dispatch campaign tools to the authenticated Proto Modal backend."""

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(
        self,
        tool_key: str,
        inputs: dict[str, Any],
        config: dict[str, Any],
    ) -> dict[str, Any]:
        result = asyncio.run(
            invoke_runtime_mcp(
                "proto",
                ("run_tool",),
                {
                    "tool_key": tool_key,
                    "inputs": inputs,
                    "config": config,
                    "output_dir": str(self.output_dir),
                    "run_on": "modal",
                },
            )
        )
        payload = _unwrap_mcp_payload(_json_safe(result))
        if not isinstance(payload, dict):
            raise RuntimeMCPError(
                f"Proto {tool_key} returned an unexpected payload: {type(payload).__name__}"
            )
        return payload


def _paperclip_cli(command: str) -> dict[str, Any]:
    parts = shlex.split(command)
    allowed = {
        "skill",
        "routines",
        "search",
        "map",
        "reduce",
        "results",
        "filter",
        "lookup",
        "grep",
        "scan",
        "sql",
        "cat",
        "head",
    }
    if not parts or parts[0] not in allowed:
        raise ValueError(f"unsupported Paperclip command: {parts[0] if parts else '<empty>'}")
    try:
        completed = subprocess.run(
            ["paperclip", *parts],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
    except OSError as exc:
        raise RuntimeMCPError(
            "Paperclip is not authenticated for the LangGraph worker. Configure "
            "PAPERCLIP_MCP_BEARER_TOKEN or install and authenticate the Paperclip CLI."
        ) from exc
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeMCPError(f"Paperclip CLI failed ({completed.returncode}): {detail}")
    return {"transport": "cli", "command": command, "output": completed.stdout}


async def _paperclip_command(command: str) -> dict[str, Any]:
    if runtime_mcp_status()["paperclip"]["authentication"] == "bearer":
        try:
            result = await invoke_runtime_mcp("paperclip", ("paperclip",), {"command": command})
            return {"transport": "mcp", "command": command, "output": _json_safe(result)}
        except RuntimeMCPError as exc:
            mcp_warning = str(exc)
    else:
        mcp_warning = "PAPERCLIP_MCP_BEARER_TOKEN is not configured."

    cli_result = await asyncio.to_thread(_paperclip_cli, command)
    cli_result["mcp_warning"] = mcp_warning
    return cli_result


def _paperclip_result_id(payload: Any, prefix: str) -> str | None:
    """Extract a saved Paperclip result id from MCP or CLI output."""

    matches = re.findall(rf"\b{re.escape(prefix)}_[A-Za-z0-9]+\b", json.dumps(payload))
    return matches[-1] if matches else None


def _paperclip_map_citations(payload: Any) -> list[dict[str, str]]:
    """Convert Paperclip CLI map line markers into durable citation URLs."""

    text = str(payload.get("output", "")) if isinstance(payload, dict) else str(payload)
    citations: list[dict[str, str]] = []
    for section in re.split(r"\n\s*✓\s+", text):
        paper_match = re.search(r"\b(PMC\d+)\b", section)
        if paper_match is None:
            continue
        line_markers = list(
            dict.fromkeys(
                re.findall(r"\bL\d+(?:-L?\d+)?\b", section)
            )
        )
        if not line_markers:
            continue
        paper_id = paper_match.group(1)
        citations.append(
            {
                "claim": f"Paperclip mapped line evidence for {paper_id}.",
                "url": (
                    f"https://paperclip.gxl.ai/citations/papers/{paper_id}"
                    f"#{','.join(line_markers)}"
                ),
            }
        )
    return citations


def _citation_rows(payload: Any) -> list[dict[str, str]]:
    strings: list[str] = []

    def collect(value: Any) -> None:
        if isinstance(value, str):
            strings.append(value)
        elif isinstance(value, dict):
            for item in value.values():
                collect(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                collect(item)

    collect(payload)
    text = "\n".join(strings)
    urls = sorted(
        set(
            re.findall(
                r"https?://(?:paperclip\.gxl\.ai|pmc\.ncbi\.nlm\.nih\.gov|"
                r"www\.ncbi\.nlm\.nih\.gov|www\.rcsb\.org|doi\.org)/[^\s\"'<>]+",
                text,
            )
        )
    )
    return [{"claim": "Paperclip research result", "url": url.rstrip(".,);")} for url in urls]


@tool
async def research_design_objective(objective: str) -> dict[str, Any]:
    """Research a protein-design objective through Paperclip before drafting a design spec."""

    objective = objective.strip()
    if len(objective) < 10:
        raise ValueError("objective must describe a target and intended design")

    await _paperclip_command("skill")
    routed: dict[str, Any] | None = None
    try:
        routed = await _paperclip_command(
            f"routines route {shlex.quote('protein binder target interface design')}"
        )
    except RuntimeMCPError:
        routed = None

    literature = await _paperclip_command(
        f"search -s pmc -n 8 {shlex.quote(objective + ' structure interface hotspot')}"
    )
    literature_id = _paperclip_result_id(literature, "s")
    if literature_id is None:
        raise RuntimeMCPError(
            "Paperclip literature search did not return a saved result set; "
            "the application cannot build line-pinned design evidence."
        )
    mapped_literature = await _paperclip_command(
        "map --from "
        f"{literature_id} "
        + shlex.quote(
            "Extract target identity, experimentally resolved structures, interface "
            "residues or hotspots, binder constraints, and limitations. Cite every "
            "claim with Paperclip line-pinned citation URLs."
        )
    )
    await _paperclip_command("skill proteins")
    proteins = await _paperclip_command(
        f"search -s proteins -n 8 {shlex.quote(objective)}"
    )
    searches = [literature, mapped_literature, proteins]
    citations = _citation_rows(searches)
    citations.extend(_paperclip_map_citations(mapped_literature))
    citations = list({citation["url"]: citation for citation in citations}.values())
    line_pinned = [
        citation
        for citation in citations
        if "paperclip.gxl.ai/citations/" in citation["url"] and "#L" in citation["url"]
    ]
    if not line_pinned:
        raise RuntimeMCPError(
            "Paperclip returned no line-pinned evidence URLs; design-spec drafting "
            "is blocked rather than substituting an uncited claim."
        )
    return write_run_artifact(
        run_id=new_run_id("target-research"),
        kind="target_research",
        payload={
            "title": "Paperclip target and interface research",
            "objective": objective,
            "status": "complete",
            "research": {
                "routine": routed,
                "searches": searches,
                "literature_result_id": literature_id,
            },
            "claim_boundary": (
                "Search results support design-spec drafting; target residues still require "
                "verification against the pinned structure."
            ),
            "citations": citations,
        },
    )


@tool
async def paperclip_evidence_command(command: str) -> dict[str, Any]:
    """Run a read-only Paperclip map, result, lookup, or line-reading command."""

    parts = shlex.split(command)
    if not parts:
        raise ValueError("command cannot be empty")
    if parts[0] != "skill":
        await _paperclip_command("skill")
    result = await _paperclip_command(command)
    citations = _citation_rows(result)
    return write_run_artifact(
        run_id=new_run_id("paperclip-evidence"),
        kind="paperclip_evidence",
        payload={
            "title": f"Paperclip evidence: {parts[0]}",
            "command": command,
            "result": result,
            "claim_boundary": (
                "This artifact retains retrieved evidence. Only line-pinned passages "
                "may support final design-spec claims."
            ),
            "citations": citations,
        },
    )


@tool
async def inspect_proto_design_tools() -> dict[str, Any]:
    """Verify Proto workspace and exact RFdiffusion3, ProteinMPNN, and AlphaFold2 schemas."""

    status = runtime_mcp_status()
    if not status["proto"]["configured"]:
        raise RuntimeMCPError("Proto runtime MCP is not configured")
    workspace = await invoke_runtime_mcp("proto", ("workspace_info",), {})
    searches = {}
    for query in (
        "RFdiffusion3 binder design",
        "ProteinMPNN inverse folding",
        "AlphaFold2",
        "PyMOL RMSD alignment",
    ):
        searches[query] = _json_safe(
            await invoke_runtime_mcp(
                "proto",
                ("search_tools",),
                {"query": query, "deployed_only": False, "limit": 5},
            )
        )
    schemas = {}
    deployments = {}
    for tool_key in (
        "rfdiffusion3-design",
        "proteinmpnn-sample",
        "alphafold2-prediction",
        "pymol-rmsd-alignment",
    ):
        schemas[tool_key] = _json_safe(
            await invoke_runtime_mcp(
                "proto",
                ("get_tool_schema",),
                {"tool_key": tool_key},
            )
        )
        deployments[tool_key] = _json_safe(
            await invoke_runtime_mcp(
                "proto",
                ("get_tool_info",),
                {"tool_key": tool_key},
            )
        )
    return write_run_artifact(
        run_id=new_run_id("proto-catalog"),
        kind="proto_tool_validation",
        payload={
            "title": "Proto design tool validation",
            "status": "complete",
            "runtime": status["proto"],
            "workspace": _json_safe(workspace),
            "searches": searches,
            "schemas": schemas,
            "deployments": deployments,
            "claim_boundary": "Tool discovery validates contracts but does not start compute.",
            "citations": [
                {
                    "claim": "Proto tool contracts were discovered through runtime MCP.",
                    "url": "https://mcp.evodesign.org/mcp",
                }
            ],
        },
    )


@tool
def immuno_architecture_status() -> dict[str, Any]:
    """Inspect required immunogenicity artifacts without inventing provider output."""

    checks = {
        "hla_panel": _path_status(ROOT / "docs" / "hla_class_ii_panel.v1.json"),
        "fusion_rule": _path_status(ROOT / "docs" / "fusion_rule.v1.json"),
        "self_proteome": _path_status(
            ROOT / "data" / "processed" / "self_proteome.parquet"
        ),
        "design_spec": _path_status(ROOT / "docs" / "design_spec.json"),
        "response_contract": _path_status(ROOT / "docs" / "RESPONSE_ADAPTER_CONTRACT.md"),
        "team_mhci_chao1_checkpoint": _path_status(
            ROOT / "models" / "chao1" / "cv5_heads.pkl 2"
        ),
    }
    response_ready = False
    response_reason = "No immutable teammate response artifact was supplied."
    status = {
        "status": "ready" if all(row["available"] for row in checks.values()) else "partial",
        "checks": checks,
        "providers": {
            "runtime_mcp": runtime_mcp_status(),
            "netmhciipan": {
                "version": "4.3",
                "channels": ["EL", "BA"],
                "mode": "IEDB Tools API with immutable local cache",
            },
            "response_model": {
                "available": response_ready,
                "fallback": "UnavailableResponseAdapter",
                "reason": response_reason,
            },
            "mhc_i_processing_surrogate": {
                "available": checks["team_mhci_chao1_checkpoint"]["available"],
                "checkpoint": checks["team_mhci_chao1_checkpoint"]["path"],
                "fusion_policy": "separate evidence lane; excluded from MHC-II fusion",
            },
        },
        "claim_boundary": (
            "Computational triage only; no measured immune-response probability is claimed."
        ),
    }
    return write_run_artifact(
        run_id=new_run_id("architecture"),
        kind="architecture_status",
        payload={
            "title": "Immunogenicity architecture readiness",
            "claim_boundary": status["claim_boundary"],
            "architecture": status,
            "citations": [
                {
                    "claim": "NetMHCIIpan 4.3 EL and BA are exposed through IEDB Tools.",
                    "url": "https://tools.iedb.org/main/tools-api/",
                }
            ],
        },
    )


@tool
@traceable(name="screen_candidate_immunogenicity", run_type="tool")
def screen_candidate(
    sequence: str,
    candidate_id: str = "workbench-candidate",
    response_artifact: str | None = None,
    mhci_surrogate_checkpoint: str | None = None,
    mhci_surrogate_checkpoints: list[str] | None = None,
    mhci_netmhcpan_checkpoint: str | None = None,
    source_metadata: dict[str, Any] | None = None,
    structure_path: str | None = None,
    structure_pdb_id: str | None = None,
    structure_chain_id: str = "A",
) -> dict[str, Any]:
    """Screen one sequence and optionally map residue evidence to a PDB chain."""

    agent = _screening_agent(response_artifact)
    if structure_pdb_id:
        if structure_path:
            raise ValueError("provide structure_path or structure_pdb_id, not both")
        structure_path = str(_cached_rcsb_pdb(structure_pdb_id).relative_to(ROOT))
    assessment = _assess_candidate(
        agent,
        candidate_id=candidate_id,
        sequence=sequence,
        structure_path=structure_path,
        structure_chain_id=structure_chain_id,
    )
    _attach_mhci_surrogates(
        assessment,
        _mhci_surrogates(
            mhci_surrogate_checkpoint,
            mhci_surrogate_checkpoints,
            mhci_netmhcpan_checkpoint,
        ),
    )
    payload = {
        "title": f"Immunogenicity screen: {candidate_id}",
        "schema_version": "1.0.0",
        "claim_boundary": (
            "Computational candidate triage from separate response, MHC-II EL, "
            "MHC-II BA, processing, tolerance, and optional MHC-I surrogate evidence."
        ),
        "source_metadata": source_metadata or {},
        "citations": [
            {
                "claim": "NetMHCIIpan integrates binding-affinity and eluted-ligand data.",
                "url": "https://paperclip.gxl.ai/citations/papers/PMC7319546#L11,L16,L36-L42",
            },
            {
                "claim": "EL and BA are called through the IEDB Tools API.",
                "url": "https://tools.iedb.org/main/tools-api/",
            },
        ],
        "assessment": assessment.model_dump(mode="json"),
    }
    return write_run_artifact(
        run_id=new_run_id("immuno"),
        kind="candidate_assessment",
        payload=payload,
    )


@tool
def discover_scientific_skills(query: str = "") -> dict[str, Any]:
    """List installed scientific skills and their brief descriptions."""

    needle = query.strip().lower()
    rows: list[dict[str, str]] = []
    for root in SKILL_ROOTS:
        if not root.exists():
            continue
        for path in sorted(root.glob("*/SKILL.md")):
            text = path.read_text(errors="replace")
            first = next(
                (
                    line.strip("# ").strip()
                    for line in text.splitlines()
                    if line.strip() and not line.startswith("---")
                ),
                path.parent.name,
            )
            haystack = f"{path.parent.name} {first}".lower()
            if needle and needle not in haystack:
                continue
            rows.append(
                {
                    "name": path.parent.name,
                    "description": first,
                    "path": str(path.relative_to(ROOT)),
                    "source": "project" if root == SKILL_ROOTS[0] else "extended",
                }
            )
    return {"query": query, "count": len(rows), "skills": rows[:100]}


@tool
def read_scientific_skill(path: str) -> dict[str, Any]:
    """Read an installed SKILL.md after constraining access to approved skill roots."""

    target = (ROOT / path).resolve()
    allowed_root = any(root.resolve() in target.parents for root in SKILL_ROOTS)
    if target.name != "SKILL.md" or not allowed_root:
        raise ValueError("path must reference an installed SKILL.md")
    content = target.read_text()
    return {
        "path": str(target.relative_to(ROOT)),
        "content": content,
    }


@tool
def list_scientific_runs(limit: int = 20) -> dict[str, Any]:
    """List local run manifests retained as the fallback to LangSmith."""

    bounded = max(1, min(limit, 100))
    return {"runs": list_runs(bounded)}


@tool
def replay_latest_design_campaign() -> dict[str, Any]:
    """Replay the latest hash-valid completed campaign without external calls."""

    source = latest_run_artifact("design_to_screen_campaign")
    payload = json.loads(json.dumps(source["payload"]))
    manifest = payload.get("manifest", {})
    if manifest.get("status") != "completed":
        raise ValueError("latest design campaign is not complete and cannot be replayed")
    payload.update(
        {
            "title": f"Replay: {payload.get('title', source['title'])}",
            "execution_mode": "replay",
            "replay_of": source["id"],
            "replay_source_sha256": source["sha256"],
        }
    )
    return write_run_artifact(
        run_id=new_run_id("design-replay"),
        kind="design_to_screen_campaign_replay",
        payload=payload,
        parent_run_id=source["id"],
    )


@tool
def read_design_spec(spec_path: str = "docs/design_spec.json") -> dict[str, Any]:
    """Read the versioned protein design objective and generation constraints."""

    path = _resolve_repo_file(
        spec_path,
        allowed_roots=(ROOT / "docs", ROOT / "results" / "design_specs"),
    )
    payload = json.loads(path.read_text())
    return {
        "path": str(path.relative_to(ROOT)),
        "spec": payload,
    }


def _pdb_residues(path: Path) -> dict[tuple[str, int], str]:
    residues: dict[tuple[str, int], str] = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            continue
        chain = line[21].strip()
        residue_number = line[22:26].strip()
        residue_name = line[17:20].strip()
        if not chain or not residue_number.lstrip("-").isdigit():
            continue
        residues[(chain, int(residue_number))] = residue_name
    return residues


@tool
def draft_cited_design_spec(
    objective: str,
    target_name: str,
    gene: str,
    uniprot_accession: str,
    pdb_id: str,
    target_chain: str,
    structure_path: str,
    hotspots: list[str],
    source_urls: list[str],
    binder_length: int = 80,
    n_backbones: int = 8,
    sequences_per_backbone: int = 4,
) -> dict[str, Any]:
    """Create a versioned design spec from researched target facts and a pinned PDB."""

    if not source_urls or any(not url.startswith("https://") for url in source_urls):
        raise ValueError("source_urls must contain at least one HTTPS evidence citation")
    structure = _resolve_repo_file(
        structure_path,
        allowed_roots=(ROOT / "data" / "raw" / "targets",),
    )
    if not 30 <= binder_length <= 200:
        raise ValueError("binder_length must be between 30 and 200 residues")
    if not 1 <= n_backbones <= 64 or not 1 <= sequences_per_backbone <= 16:
        raise ValueError("generation budget exceeds the bounded campaign limits")

    available_residues = _pdb_residues(structure)
    hotspot_rows: list[dict[str, Any]] = []
    for hotspot in hotspots:
        match = re.fullmatch(r"([A-Za-z0-9])(-?\d+)", hotspot.strip())
        if not match:
            raise ValueError(f"invalid hotspot {hotspot!r}; expected a chain and residue, e.g. B62")
        chain, residue_text = match.groups()
        residue_number = int(residue_text)
        residue_name = available_residues.get((chain, residue_number))
        if residue_name is None:
            raise ValueError(f"hotspot {hotspot} is absent from the pinned structure")
        hotspot_rows.append(
            {
                "chain": chain,
                "residue_number": residue_number,
                "residue_name": residue_name,
            }
        )
    if not hotspot_rows or any(row["chain"] != target_chain for row in hotspot_rows):
        raise ValueError("all hotspots must be present on the requested target chain")

    run_id = new_run_id("design-spec")
    spec_id = f"{gene.lower()}-{pdb_id.lower()}-{run_id[-8:]}"
    spec = {
        "schema_version": "1.0.0",
        "spec_id": spec_id,
        "status": "awaiting_compute_approval",
        "objective": {
            "type": "protein_binder_design",
            "description": objective.strip(),
            "claim_boundary": (
                "Computational prioritization only; generated structures and scores are not "
                "experimental evidence of folding, binding, or human immune response."
            ),
        },
        "target": {
            "name": target_name,
            "gene": gene,
            "uniprot_accession": uniprot_accession,
            "pdb_id": pdb_id.upper(),
            "chain_id": target_chain,
            "structure_path": str(structure.relative_to(ROOT)),
            "structure_sha256": hashlib.sha256(structure.read_bytes()).hexdigest(),
        },
        "interface": {
            "selection_basis": "Paperclip evidence, verified against the pinned PDB",
            "hotspots": hotspot_rows,
            "hotspot_validation": "Every residue is present on the selected target chain.",
        },
        "generation": {
            "backbone_provider": {
                "tool_key": "rfdiffusion3-design",
                "input_structure": str(structure.relative_to(ROOT)),
                "target_chain": target_chain,
                "binder_length": str(binder_length),
                "select_hotspots": ",".join(hotspots),
                "n_backbones": n_backbones,
                "seed": 20260816,
            },
            "sequence_provider": {
                "tool_key": "proteinmpnn-sample",
                "chains_to_redesign": ["A"],
                "fixed_context_chains": [target_chain],
                "sequences_per_backbone": sequences_per_backbone,
                "temperature": 0.1,
                "model_choice": "proteinmpnn",
                "seed": 20260816,
            },
            "maximum_candidates": n_backbones * sequences_per_backbone,
        },
        "validation": {
            "required_before_screening": [
                "unique_sequence",
                "proteinmpnn_perplexity_recorded",
                "monomer_structure_confidence_recorded",
                "designed_vs_refolded_binder_rmsd_recorded",
                "interface_confidence_recorded",
                "clash_check_recorded",
            ],
            "literature_prior_thresholds": {
                "monomer_plddt_min": 80.0,
                "interface_pae_max": 10.0,
                "binder_rmsd_angstrom_max": 1.0,
            },
            "threshold_policy": (
                "Thresholds are computational triage priors and must not be described as "
                "proof of binding."
            ),
        },
        "screening": {
            "window_length": 15,
            "response_model_registry": True,
            "mhc_primary": ["netmhciipan_el-4.3", "netmhciipan_ba-4.3"],
            "mhc_second_opinion": "MixMHC2pred-2.0",
            "mhc_optional_challenger": "Graph-pMHC",
            "preserve_component_scores": True,
        },
        "approval": {"required": True, "approved": False},
        "sources": [
            {"title": f"Design evidence {index + 1}", "url": url}
            for index, url in enumerate(source_urls)
        ],
    }
    spec_dir = ROOT / "results" / "design_specs"
    spec_dir.mkdir(parents=True, exist_ok=True)
    spec_path = spec_dir / f"{run_id}.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n")
    return write_run_artifact(
        run_id=run_id,
        kind="design_spec",
        payload={
            "title": f"{target_name} cited design specification",
            "spec_path": str(spec_path.relative_to(ROOT)),
            "spec": spec,
            "claim_boundary": spec["objective"]["claim_boundary"],
            "citations": [
                {"claim": row["title"], "url": row["url"]} for row in spec["sources"]
            ],
        },
    )


@tool
def plan_design_campaign(
    spec_path: str = "docs/design_spec.json",
    device: str = "cuda",
) -> dict[str, Any]:
    """Build the exact Proto campaign payloads without starting GPU compute."""

    resolved_spec = _resolve_repo_file(
        spec_path,
        allowed_roots=(ROOT / "docs", ROOT / "results" / "design_specs"),
    )
    spec = json.loads(resolved_spec.read_text())
    run_id = new_run_id("design-plan")
    output_dir = ROOT / "results" / "design_campaigns" / run_id
    manifest = run_campaign(
        resolved_spec,
        output_dir,
        ProtoLocalRunner(),
        approved=False,
        device=device,
    )
    return write_run_artifact(
        run_id=run_id,
        kind="design_campaign_plan",
        payload={
            "title": f"{spec['target']['name']} design campaign plan",
            "claim_boundary": manifest.claim_boundary,
            "spec_path": str(resolved_spec.relative_to(ROOT)),
            "manifest": manifest.model_dump(mode="json"),
            "citations": [
                {"claim": row.get("title", "Design source"), "url": row["url"]}
                for row in spec.get("sources", [])
                if row.get("url")
            ],
        },
    )


@tool
@traceable(name="execute_proto_design_campaign", run_type="tool")
def execute_design_campaign(
    spec_path: str = "docs/design_spec.json",
    device: str = "cuda",
    campaign_run_id: str | None = None,
    response_artifact: str | None = None,
    mhci_surrogate_checkpoint: str | None = None,
    mhci_surrogate_checkpoints: list[str] | None = None,
    mhci_netmhcpan_checkpoint: str | None = None,
    max_screen_candidates: int = 32,
) -> dict[str, Any]:
    """Run approved Proto design, structural gates, and immunogenicity screening."""

    resolved_spec = _resolve_repo_file(
        spec_path,
        allowed_roots=(ROOT / "docs", ROOT / "results" / "design_specs"),
    )
    spec = json.loads(resolved_spec.read_text())
    run_id = campaign_run_id or new_run_id("design")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", run_id):
        raise ValueError(
            "campaign_run_id may contain only letters, numbers, dot, dash, and underscore"
        )
    output_dir = ROOT / "results" / "design_campaigns" / run_id
    manifest = run_campaign(
        resolved_spec,
        output_dir,
        ProtoMCPRunner(output_dir / "proto_outputs"),
        approved=True,
        device=device,
    )
    bounded = max(1, min(max_screen_candidates, 32))
    eligible = [
        candidate
        for candidate in manifest.candidates
        if candidate.screening_status == "eligible"
    ][:bounded]
    assessments: list[dict[str, Any]] = []
    manifest.phase_status["immunogenicity"] = "running" if eligible else "blocked"
    manifest_path = output_dir / "campaign_manifest.json"
    manifest_path.write_text(manifest.model_dump_json(indent=2) + "\n")
    if eligible:
        agent = _screening_agent(response_artifact)
        surrogates = _mhci_surrogates(
            mhci_surrogate_checkpoint,
            mhci_surrogate_checkpoints,
            mhci_netmhcpan_checkpoint,
        )
        for candidate in eligible:
            assessment = agent.assess(
                candidate.candidate_id,
                _normalized_sequence(candidate.sequence),
            )
            if candidate.refolded_structure_path:
                try:
                    structure = _resolve_repo_file(
                        candidate.refolded_structure_path,
                        allowed_roots=(ROOT / "results",),
                    )
                    assessment.structure = structure_reference_from_pdb(
                        structure,
                        sequence=candidate.sequence,
                        chain_id="A",
                        repository_root=ROOT,
                    )
                except (FileNotFoundError, ValueError) as exc:
                    assessment.warnings.append(
                        f"3D heatmap withheld because structure mapping failed: {exc}"
                    )
            candidate.screening_status = "screened"
            _attach_mhci_surrogates(assessment, surrogates)
            assessments.append(assessment.model_dump(mode="json"))

    manifest.candidate_counts["screened"] = len(assessments)
    manifest.phase_status["immunogenicity"] = "completed" if eligible else "blocked"
    if eligible:
        manifest.status = "completed"
    manifest_path.write_text(manifest.model_dump_json(indent=2) + "\n")
    citations = [
        {"claim": row.get("title", "Design source"), "url": row["url"]}
        for row in spec.get("sources", [])
        if row.get("url")
    ]
    citations.extend(
        [
            {
                "claim": "NetMHCIIpan integrates binding-affinity and eluted-ligand data.",
                "url": "https://paperclip.gxl.ai/citations/papers/PMC7319546#L11,L16,L36-L42",
            },
            {
                "claim": "NetMHCIIpan EL and BA are called through IEDB Tools.",
                "url": "https://tools.iedb.org/main/tools-api/",
            },
        ]
    )
    return write_run_artifact(
        run_id=run_id,
        kind="design_to_screen_campaign",
        payload={
            "title": f"{spec['target']['name']} design-to-screen campaign",
            "claim_boundary": manifest.claim_boundary,
            "spec_path": str(resolved_spec.relative_to(ROOT)),
            "manifest": manifest.model_dump(mode="json"),
            "assessment": assessments[0] if assessments else None,
            "assessments": assessments,
            "citations": citations,
        },
    )


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                rows.append((name, "".join(sequence)))
            name = line[1:].strip().split()[0]
            sequence = []
        else:
            sequence.append(line.strip())
    if name is not None:
        rows.append((name, "".join(sequence)))
    return rows


def _reference_campaign_metadata() -> dict[str, dict[str, Any]]:
    screening_path = ROOT / "data" / "raw" / "il7ra" / "IL7Ra_binders_screening_results.csv"
    bli_path = ROOT / "data" / "raw" / "il7ra" / "bli_binding_data.csv"
    by_prediction: dict[str, dict[str, Any]] = {}
    design_to_prediction: dict[str, str] = {}
    with screening_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            prediction = row.get("AF2_prediction", "")
            design = row.get("RFdiffusion_output", "")
            if prediction:
                by_prediction[prediction] = {
                    "reference_design_id": design,
                    "reference_monomeric": row.get("Monomeric") == "True",
                    "reference_binder_label": row.get("Binder") == "True",
                    "noise_level": row.get("Noise_level"),
                    "fold_conditioned": row.get("Fold_conditioned") == "True",
                }
            if design and prediction:
                design_to_prediction[design] = prediction

    kd_values: dict[str, list[float]] = {}
    with bli_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            prediction = design_to_prediction.get(row.get("sample_id", ""))
            if not prediction:
                continue
            try:
                kd = float(row["KD (M)"])
            except (KeyError, TypeError, ValueError):
                continue
            kd_values.setdefault(prediction, []).append(kd)
    for prediction, values in kd_values.items():
        ordered = sorted(values)
        by_prediction.setdefault(prediction, {})["reference_bli_kd_median_m"] = ordered[
            len(ordered) // 2
        ]
    return by_prediction


@tool
@traceable(name="reference_campaign_preflight", run_type="tool")
def run_reference_campaign_preflight(
    limit: int = 3,
    response_artifact: str | None = None,
    mhci_surrogate_checkpoint: str | None = None,
    mhci_surrogate_checkpoints: list[str] | None = None,
    mhci_netmhcpan_checkpoint: str | None = None,
) -> dict[str, Any]:
    """Screen a bounded subset of the 95 historical sequences without GPU generation."""

    bounded = max(1, min(limit, 10))
    sequences = _read_fasta(REFERENCE_FASTA)
    metadata = _reference_campaign_metadata()
    agent = _screening_agent(response_artifact)
    surrogates = _mhci_surrogates(
        mhci_surrogate_checkpoint,
        mhci_surrogate_checkpoints,
        mhci_netmhcpan_checkpoint,
    )
    assessments = []
    for candidate_id, sequence in sequences[:bounded]:
        assessment = agent.assess(candidate_id, _normalized_sequence(sequence))
        _attach_mhci_surrogates(assessment, surrogates)
        assessments.append(
            {
                "reference": metadata.get(candidate_id, {}),
                **assessment.model_dump(mode="json"),
            }
        )
    spec = json.loads(DEFAULT_SPEC_PATH.read_text())
    return write_run_artifact(
        run_id=new_run_id("reference-preflight"),
        kind="reference_campaign_preflight",
        payload={
            "title": "IL-7Rα no-GPU screening preflight",
            "claim_boundary": (
                "Historical sequences validate the screening and reporting handoff only; "
                "they are not inputs to RFdiffusion3 or evidence of new design generation."
            ),
            "reference_dataset": {
                "fasta": str(REFERENCE_FASTA.relative_to(ROOT)),
                "total_sequences": len(sequences),
                "screened_sequences": len(assessments),
            },
            "assessment": assessments[0] if assessments else None,
            "assessments": assessments,
            "citations": [
                {"claim": row.get("title", "Reference source"), "url": row["url"]}
                for row in spec.get("sources", [])
                if row.get("url")
            ],
        },
    )


SCIENTIFIC_TOOLS = [
    research_design_objective,
    paperclip_evidence_command,
    inspect_proto_design_tools,
    immuno_architecture_status,
    screen_candidate,
    discover_scientific_skills,
    read_scientific_skill,
    list_scientific_runs,
    replay_latest_design_campaign,
    read_design_spec,
    draft_cited_design_spec,
    plan_design_campaign,
    execute_design_campaign,
    run_reference_campaign_preflight,
]
