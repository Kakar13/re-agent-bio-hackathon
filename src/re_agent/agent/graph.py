"""LangGraph scientific-agent runtime with deterministic review gates."""

from __future__ import annotations

import json
import os
import re
from typing import Any, Literal
from uuid import uuid4

from langchain_anthropic import ChatAnthropic
from langchain_core.messages import AIMessage, SystemMessage, ToolMessage
from langgraph.graph import END, START, StateGraph
from langgraph.types import interrupt

from re_agent.agent.run_store import ROOT, new_run_id, write_run_artifact
from re_agent.agent.state import ScientificAgentState
from re_agent.agent.tools import SCIENTIFIC_TOOLS

os.environ.setdefault("LANGSMITH_PROJECT", "re-agent-scientific-workbench")

SYSTEM_PROMPT = """You are re:AGENT, an inspectable protein-design and immunogenicity
screening scientist. Use tools before making empirical claims. Keep response propensity,
NetMHCIIpan EL, NetMHCIIpan BA, processing, tolerance, and optional MHC-I surrogate
evidence separate. Never invent missing predictor output. The primary chao1 checkpoint
under models/chao1 is an HLA-A*02:01 MHC-I processing surrogate, not a response model;
pass it as mhci_surrogate_checkpoints when selected and never fuse it into MHC-II.
A missing response artifact must remain unavailable and the combined rank must be
withheld. Summarize tool rationale, measured output, uncertainty, citations, and the
deterministic reviewer result. Do not expose hidden chain-of-thought; provide concise
evidence-backed rationale summaries instead. Do not make clinical or FDA claims. Ask
for explicit approval before expensive RFdiffusion, ProteinMPNN, AlphaFold2, Boltz, or
other remote compute.

For a natural-language protein-design objective, follow this artifact chain:
1. Call research_design_objective to gather Paperclip target, structure, interface, and
   constraint evidence.
2. Use paperclip_evidence_command to map the returned literature set and read the exact
   line-pinned passages that support target, PDB, interface, hotspot, and constraint claims.
3. Call inspect_proto_design_tools to verify exact Proto contracts.
4. Call draft_cited_design_spec with only line-supported values and a pinned local PDB.
5. Call plan_design_campaign with the returned spec_path.
6. Explain the plan and call execute_design_campaign only when the user asks to proceed;
   the runtime will pause for explicit compute approval.
7. The execution tool applies structural gates before automatically screening eligible
   candidates and returns separate risk channels plus exact-mapped 3D tracks.

The 95 IL-7Ralpha sequences are historical preflight/evaluation data. Use
run_reference_campaign_preflight to validate the no-GPU screening handoff, never as
RFdiffusion3 generation input.
"""

TOOL_MAP = {scientific_tool.name: scientific_tool for scientific_tool in SCIENTIFIC_TOOLS}
CUSTOM_MODEL_CHECKPOINTS = {
    "chao1": "models/chao1/cv5_heads.pkl 2",
}
PROFILE_CHECKPOINTS = {
    "mhc_ii_plus_chao1": [CUSTOM_MODEL_CHECKPOINTS["chao1"]],
}
# The NetMHCpan student replaces chao1's MHCflurry-derived MHC lane. It is passed
# alongside the chao1 checkpoint rather than instead of it, because cleavage and
# TAP still come from chao1; only the binding and presentation lanes change.
NETMHCPAN_STUDENT_CHECKPOINT = "models/a0201-netmhcpan-pda-cv5-v4/checkpoint"


def _netmhcpan_student_checkpoint() -> str | None:
    """Return the student ensemble path when it is present on disk.

    Absent checkpoints degrade to the chao1-only lane instead of failing the
    run, so a fresh clone without the downloaded weights still screens.
    """

    path = ROOT / NETMHCPAN_STUDENT_CHECKPOINT
    if (path / "deployment_manifest.json").exists():
        return NETMHCPAN_STUDENT_CHECKPOINT
    return None


def _selected_custom_checkpoints(state: ScientificAgentState) -> list[str]:
    return PROFILE_CHECKPOINTS.get(state.get("screening_profile", ""), [])


def _prompt_custom_checkpoints(content: str) -> list[str]:
    if "chao1" in content or any(
        term in content for term in ("team model", "main model", "mhc-i", "mhci surrogate")
    ):
        return [CUSTOM_MODEL_CHECKPOINTS["chao1"]]
    return []


def _profiled_tool_args(
    state: ScientificAgentState,
    tool_name: str,
    args: dict[str, Any],
) -> dict[str, Any]:
    profiled = dict(args)
    checkpoints = _selected_custom_checkpoints(state)
    if checkpoints and tool_name in {
        "screen_candidate",
        "run_reference_campaign_preflight",
        "execute_design_campaign",
    }:
        profiled.setdefault("mhci_surrogate_checkpoints", checkpoints)
        student = _netmhcpan_student_checkpoint()
        if student:
            profiled.setdefault("mhci_netmhcpan_checkpoint", student)
    return profiled


def _agent_node(state: ScientificAgentState) -> dict[str, Any]:
    pipeline_request = state.get("pipeline_request")
    if pipeline_request:
        objective = str(pipeline_request.get("objective", "")).strip()
        if len(objective) < 20:
            return {
                "pipeline_request": None,
                "messages": [
                    AIMessage(
                        content=(
                            "Describe the binder target and at least one design constraint "
                            "(for example length, interface, or excluded residues)."
                        )
                    )
                ],
            }
        return {
            "pipeline_request": None,
            "messages": [
                AIMessage(
                    content="",
                    tool_calls=[
                        {
                            "name": "research_design_objective",
                            "args": {"objective": objective},
                            "id": f"pipeline-research-{uuid4().hex[:10]}",
                            "type": "tool_call",
                        }
                    ],
                )
            ],
        }

    direct_request = state.get("direct_screen_request")
    if direct_request:
        sequence = re.sub(r"\s+", "", str(direct_request.get("sequence", ""))).upper()
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{15,}", sequence):
            return {
                "direct_screen_request": None,
                "messages": [
                    AIMessage(
                        content=(
                            "Chao1 screening requires at least 15 standard amino-acid "
                            "letters."
                        )
                    )
                ],
            }
        student_checkpoint = _netmhcpan_student_checkpoint()
        return {
            "direct_screen_request": None,
            "messages": [
                AIMessage(
                    content="",
                    tool_calls=[
                        {
                            "name": "screen_candidate",
                            "args": {
                                "sequence": sequence,
                                "candidate_id": direct_request.get(
                                    "candidate_id",
                                    "workbench-chao1-candidate",
                                ),
                                "mhci_surrogate_checkpoints": [
                                    CUSTOM_MODEL_CHECKPOINTS["chao1"]
                                ],
                                **(
                                    {"mhci_netmhcpan_checkpoint": student_checkpoint}
                                    if student_checkpoint
                                    else {}
                                ),
                                "source_metadata": {
                                    "ui_action": "direct_chao1_screen",
                                    "selected_model": "chao1",
                                    "mhci_binding_lane": (
                                        "netmhcpan_student"
                                        if student_checkpoint
                                        else "chao1_mhcflurry"
                                    ),
                                },
                                **(
                                    {
                                        "structure_pdb_id": "9S14",
                                        "structure_chain_id": "A",
                                    }
                                    if direct_request.get("candidate_id") == "pda:9s14:0"
                                    else {}
                                ),
                            },
                            "id": f"direct-chao1-{uuid4().hex[:10]}",
                            "type": "tool_call",
                        }
                    ],
                )
            ],
        }

    force_keyless = os.environ.get("RE_AGENT_FORCE_KEYLESS", "").lower() in {
        "1",
        "true",
        "yes",
    }
    if not os.getenv("ANTHROPIC_API_KEY") or force_keyless:
        last = state["messages"][-1]
        if isinstance(last, ToolMessage):
            reviewer = state.get("reviews", [])
            suffix = (
                f" Deterministic reviewer: {reviewer[-1]['status']}."
                if reviewer
                else ""
            )
            return {
                "messages": [
                    AIMessage(
                        content=(
                            f"The local scientific tool `{last.name}` completed and its "
                            f"versioned artifact is open in the canvas.{suffix}"
                        )
                    )
                ]
            }
        raw_content = str(getattr(last, "content", ""))
        content = raw_content.lower()
        fallback_tool: tuple[str, dict[str, Any]] | None = None
        sequence_match = re.search(r"\b[ACDEFGHIKLMNPQRSTVWY]{15,}\b", raw_content.upper())
        if "screen" in content and sequence_match:
            structure_match = re.search(
                r"(?:structure|pdb)\s+([A-Za-z0-9_./-]+\.pdb)\b",
                raw_content,
                flags=re.IGNORECASE,
            )
            chain_match = re.search(r"\bchain\s+([A-Za-z0-9]+)\b", raw_content, re.IGNORECASE)
            args: dict[str, Any] = {
                "sequence": sequence_match.group(0),
                "candidate_id": "keyless-workbench-candidate",
            }
            checkpoints = _selected_custom_checkpoints(state) or _prompt_custom_checkpoints(
                content
            )
            if checkpoints:
                args["mhci_surrogate_checkpoints"] = checkpoints
            if structure_match:
                args["structure_path"] = structure_match.group(1)
                args["structure_chain_id"] = chain_match.group(1) if chain_match else "A"
            fallback_tool = (
                "screen_candidate",
                args,
            )
        elif "architecture" in content or "readiness" in content:
            fallback_tool = ("immuno_architecture_status", {})
        elif "reference" in content and ("preflight" in content or "95" in content):
            args = {"limit": 1}
            checkpoints = _selected_custom_checkpoints(state) or _prompt_custom_checkpoints(
                content
            )
            if checkpoints:
                args["mhci_surrogate_checkpoints"] = checkpoints
            fallback_tool = ("run_reference_campaign_preflight", args)
        elif "replay" in content and ("design" in content or "campaign" in content):
            fallback_tool = ("replay_latest_design_campaign", {})
        elif "proto" in content and ("inspect" in content or "tool" in content):
            fallback_tool = ("inspect_proto_design_tools", {})
        elif "research" in content and ("target" in content or "design" in content):
            fallback_tool = ("research_design_objective", {"objective": raw_content})
        elif "campaign plan" in content:
            fallback_tool = ("plan_design_campaign", {"device": "cuda"})
        elif "design spec" in content or "design specification" in content:
            fallback_tool = ("read_design_spec", {})
        elif "skill" in content:
            fallback_tool = ("discover_scientific_skills", {"query": "protein structure"})
        if fallback_tool:
            name, args = fallback_tool
            return {
                "messages": [
                    AIMessage(
                        content="",
                        tool_calls=[
                            {
                                "name": name,
                                "args": args,
                                "id": f"fallback-{uuid4().hex[:10]}",
                                "type": "tool_call",
                            }
                        ],
                    )
                ]
            }
        return {
            "messages": [
                AIMessage(
                    content=(
                        "ANTHROPIC_API_KEY is not configured. Keyless mode supports the "
                        "readiness, reference-preflight, target-research, Proto inspection, "
                        "design-specification, skills, and campaign-plan actions; "
                        "open-ended agent planning requires the key in the repo-root .env."
                    )
                )
            ]
        }
    model = ChatAnthropic(
        model=os.getenv("RE_AGENT_MODEL", "claude-sonnet-4-6"),
        temperature=0,
        max_tokens=4096,
    ).bind_tools(SCIENTIFIC_TOOLS)
    context: list[Any] = [SystemMessage(content=SYSTEM_PROMPT)]
    selected_checkpoints = _selected_custom_checkpoints(state)
    if selected_checkpoints:
        context.append(
            SystemMessage(
                content=(
                    "The user selected a custom-model screening profile. Pass "
                    f"mhci_surrogate_checkpoints={selected_checkpoints!r} to candidate, "
                    "reference-preflight, and campaign-execution screening tools. Keep "
                    "every custom MHC-I lane outside MHC-II fusion."
                )
            )
        )
    if state.get("reviews"):
        context.append(
            SystemMessage(
                content="Latest deterministic reviewer result:\n"
                + json.dumps(state["reviews"][-1], default=str)
            )
        )
    response = model.invoke(context + state["messages"])
    return {"messages": [response]}


def _route_agent(state: ScientificAgentState) -> Literal["tools", "__end__"]:
    last = state["messages"][-1]
    if isinstance(last, AIMessage) and last.tool_calls:
        return "tools"
    return END


def _compute_approval_summary(args: dict[str, Any]) -> str:
    spec_path = (ROOT / str(args.get("spec_path", "docs/design_spec.json"))).resolve()
    try:
        if ROOT not in spec_path.parents:
            raise ValueError("spec outside repository")
        spec = json.loads(spec_path.read_text())
        backbone_count = int(spec["generation"]["backbone_provider"]["n_backbones"])
        sequences_per_backbone = int(
            spec["generation"]["sequence_provider"]["sequences_per_backbone"]
        )
        candidates = backbone_count * sequences_per_backbone
        return (
            f"Run {backbone_count} RFdiffusion3 backbone(s), up to {candidates} "
            f"ProteinMPNN candidate(s), and up to {candidates * 2} AlphaFold2 "
            "monomer/complex predictions on Modal. Structural gates run before "
            "immunogenicity screening. This may incur cost."
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return (
            "Run RFdiffusion3, ProteinMPNN, AlphaFold2 monomer/complex prediction, "
            "and structural validation on Modal before screening. This may incur cost."
        )


async def _tools_node(state: ScientificAgentState) -> dict[str, Any]:
    last = state["messages"][-1]
    if not isinstance(last, AIMessage):
        return {"messages": [], "artifacts": []}
    messages: list[ToolMessage] = []
    artifacts: list[dict[str, Any]] = []
    execute_call = next(
        (call for call in last.tool_calls if call["name"] == "execute_design_campaign"),
        None,
    )
    compute_approved = True
    if execute_call is not None:
        execute_args = _profiled_tool_args(
            state,
            execute_call["name"],
            execute_call["args"],
        )
        approval = interrupt(
            {
                "type": "compute_approval",
                "tool": execute_call["name"],
                "args": execute_args,
                "summary": _compute_approval_summary(execute_args),
            }
        )
        compute_approved = isinstance(approval, dict) and bool(approval.get("approved"))

    for call in last.tool_calls:
        name = call["name"]
        tool_call_id = call["id"]
        selected = TOOL_MAP.get(name)
        args = _profiled_tool_args(state, name, call["args"])
        if selected is None:
            result: Any = {"error": f"unknown tool: {name}"}
        else:
            try:
                if name == "execute_design_campaign" and not compute_approved:
                    result = {
                        "status": "cancelled",
                        "reason": "GPU compute was not approved.",
                    }
                else:
                    result = await selected.ainvoke(args)
            except Exception as exc:
                result = {"error": f"{type(exc).__name__}: {exc}"}
        if isinstance(result, dict) and {"id", "kind", "payload"} <= result.keys():
            artifacts.append(result)
        messages.append(
            ToolMessage(
                content=json.dumps(result, default=str),
                name=name,
                tool_call_id=tool_call_id,
            )
        )
    return {"messages": messages, "artifacts": artifacts}


def _review_artifact(artifact: dict[str, Any]) -> dict[str, Any]:
    payload = artifact.get("payload", {})
    checks: list[dict[str, Any]] = []

    def add(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": passed, "detail": detail})

    add("citations_present", bool(payload.get("citations")), "Artifact must cite methods.")
    text = json.dumps(payload).lower()
    add("claim_boundary", "fda" not in text, "Artifact must avoid prohibited FDA framing.")

    assessment = payload.get("assessment")
    if assessment:
        mhc_results = assessment.get("mhc_results", [])
        primary = next(
            (row for row in mhc_results if row.get("provider_id") == "netmhciipan"),
            None,
        )
        add(
            "netmhciipan_ok",
            bool(primary and primary.get("status") == "ok"),
            "Primary NetMHCIIpan evidence should be available.",
        )
        add(
            "hla_panel_complete",
            bool(primary and len(primary.get("supported_alleles", [])) == 18),
            "The frozen panel contains 18 class-II alleles.",
        )
        response_results = assessment.get("response_results", [])
        response_unavailable = any(
            row.get("status") == "unavailable" for row in response_results
        )
        add(
            "placeholder_response_gating",
            not response_unavailable or assessment.get("combined_rank_score") is None,
            "Unavailable response evidence must force a null combined rank.",
        )
        surrogate_results = assessment.get("mhc_i_surrogate_results", [])
        if surrogate_results:
            sequence_length = len(assessment.get("sequence", ""))
            add(
                "mhci_surrogate_separate",
                all(
                    row.get("provenance", {}).get("capability")
                    == "mhc_i_processing_surrogate"
                    for row in surrogate_results
                )
                and not any(
                    key.lower().startswith("mhci_") or "surrogate" in key.lower()
                    for key in assessment.get("component_summary", {})
                    .get("proxy_values", {})
                ),
                "The MHC-I surrogate must remain outside MHC-II late fusion.",
            )
            add(
                "mhci_spatial_tracks_aligned",
                all(
                    len(values) == sequence_length
                    for row in surrogate_results
                    for values in row.get("spatial_tracks", {}).values()
                ),
                "MHC-I surrogate tracks must match the parent sequence length.",
            )

    manifest = payload.get("manifest")
    if manifest:
        candidates = manifest.get("candidates", [])
        screened = [row for row in candidates if row.get("screening_status") == "screened"]
        improperly_screened = [
            row
            for row in screened
            if row.get("validation_status") != "pass"
        ]
        add(
            "structural_gate_enforced",
            not improperly_screened,
            "Only candidates that pass every recorded structural gate may be screened.",
        )
        add(
            "candidate_counts_reconcile",
            manifest.get("candidate_counts", {}).get("screened", len(screened))
            == len(payload.get("assessments", [])),
            "Campaign and immunogenicity artifact counts must agree.",
        )
        if candidates:
            add(
                "validation_checks_recorded",
                all(bool(row.get("validation_checks")) for row in candidates),
                "Every generated candidate must retain its structural validation checks.",
            )

    failed = [row for row in checks if not row["passed"]]
    return {
        "artifact_id": artifact.get("id"),
        "status": "pass" if not failed else "fail",
        "checks": checks,
        "summary": (
            "All deterministic scientific gates passed."
            if not failed
            else f"{len(failed)} deterministic gate(s) failed."
        ),
    }


def _reviewer_node(state: ScientificAgentState) -> dict[str, Any]:
    artifacts = state.get("artifacts", [])
    if not artifacts:
        return {"reviews": []}
    source = artifacts[-1]
    if any(review.get("artifact_id") == source.get("id") for review in state.get("reviews", [])):
        return {"reviews": []}
    review = _review_artifact(source)
    reviewer_artifact = write_run_artifact(
        run_id=new_run_id("review"),
        kind="review",
        payload={
            "title": f"Reviewer: {source.get('title', source.get('id'))}",
            **review,
        },
        parent_run_id=source.get("id"),
    )
    return {"reviews": [review], "artifacts": [reviewer_artifact]}


builder = StateGraph(ScientificAgentState)
builder.add_node("agent", _agent_node)
builder.add_node("tools", _tools_node)
builder.add_node("reviewer", _reviewer_node)
builder.add_edge(START, "agent")
builder.add_conditional_edges("agent", _route_agent, {"tools": "tools", END: END})
builder.add_edge("tools", "reviewer")
builder.add_edge("reviewer", "agent")

graph = builder.compile()
