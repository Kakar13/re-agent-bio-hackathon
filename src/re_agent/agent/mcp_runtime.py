"""Runtime MCP clients used by the LangGraph process.

Cursor and Pi own separate MCP sessions.  This module deliberately creates
LangGraph-local sessions so deployed graph workers do not depend on an IDE.
"""

from __future__ import annotations

import os
import shutil
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from langchain_mcp_adapters.client import MultiServerMCPClient

ROOT = Path(__file__).resolve().parents[3]


class RuntimeMCPError(RuntimeError):
    """Raised when a configured runtime MCP server cannot satisfy a tool call."""


def runtime_mcp_connections() -> dict[str, dict[str, Any]]:
    """Build secret-safe MCP connection settings from the graph environment."""

    connections: dict[str, dict[str, Any]] = {}

    paperclip_headers: dict[str, str] = {}
    paperclip_token = os.getenv("PAPERCLIP_MCP_BEARER_TOKEN")
    if paperclip_token:
        paperclip_headers["Authorization"] = f"Bearer {paperclip_token}"
    paperclip: dict[str, Any] = {
        "transport": "http",
        "url": os.getenv("PAPERCLIP_MCP_URL", "https://paperclip.gxl.ai/mcp"),
    }
    if paperclip_headers:
        paperclip["headers"] = paperclip_headers
    connections["paperclip"] = paperclip

    proto_repo_value = os.getenv("PROTO_TOOLS_REPO")
    proto_repo = Path(proto_repo_value).expanduser().resolve() if proto_repo_value else None
    proto_key = os.getenv("PROTO_API_KEY")
    if proto_repo is not None:
        pyproject = proto_repo / "pyproject.toml"
        if not pyproject.is_file():
            raise RuntimeMCPError(
                f"PROTO_TOOLS_REPO does not contain pyproject.toml: {proto_repo}"
            )
        uv = shutil.which("uv")
        if uv is None:
            raise RuntimeMCPError("PROTO_TOOLS_REPO requires the `uv` executable")
        connections["proto"] = {
            "transport": "stdio",
            "command": uv,
            "args": [
                "run",
                "--directory",
                str(proto_repo),
                "--extra",
                "mcp",
                "proto-tools-mcp",
                "--device",
                os.getenv("RE_AGENT_PROTO_CATALOG_DEVICE", "modal"),
            ],
        }
    elif proto_key:
        connections["proto"] = {
            "transport": "http",
            "url": os.getenv("PROTO_MCP_URL", "https://mcp.evodesign.org/mcp"),
            "headers": {"Authorization": f"Bearer {proto_key}"},
        }
    else:
        executable = ROOT / ".venv" / "bin" / "proto-tools-mcp"
        if executable.exists():
            connections["proto"] = {
                "transport": "stdio",
                "command": str(executable),
                "args": ["--device", os.getenv("RE_AGENT_PROTO_CATALOG_DEVICE", "modal")],
                "env": {"PYTHONPATH": str(ROOT / "src")},
            }

    return connections


def runtime_mcp_status() -> dict[str, Any]:
    """Return connection readiness without exposing headers or credentials."""

    connections = runtime_mcp_connections()
    return {
        "paperclip": {
            "configured": bool(os.getenv("PAPERCLIP_MCP_BEARER_TOKEN"))
            or shutil.which("paperclip") is not None,
            "mcp_endpoint_configured": "paperclip" in connections,
            "cli_available": shutil.which("paperclip") is not None,
            "transport": connections.get("paperclip", {}).get("transport"),
            "authentication": (
                "bearer"
                if os.getenv("PAPERCLIP_MCP_BEARER_TOKEN")
                else "external_or_cli_fallback"
            ),
        },
        "proto": {
            "configured": "proto" in connections,
            "transport": connections.get("proto", {}).get("transport"),
            "authentication": (
                "local_modal"
                if os.getenv("PROTO_TOOLS_REPO")
                else "bearer"
                if os.getenv("PROTO_API_KEY")
                else "local_stdio"
            ),
            "source": (
                "monorepo"
                if os.getenv("PROTO_TOOLS_REPO")
                else "hosted"
                if os.getenv("PROTO_API_KEY")
                else "installed"
            ),
        },
    }


async def invoke_runtime_mcp(
    server_name: str,
    preferred_tool_names: Iterable[str],
    arguments: dict[str, Any],
) -> Any:
    """Discover and invoke one tool through a LangGraph-owned MCP session."""

    connections = runtime_mcp_connections()
    connection = connections.get(server_name)
    if connection is None:
        raise RuntimeMCPError(f"runtime MCP server is not configured: {server_name}")

    client = MultiServerMCPClient({server_name: connection}, handle_tool_errors=False)
    try:
        tools = await client.get_tools(server_name=server_name)
    except BaseException as exc:
        raise RuntimeMCPError(
            f"{server_name} MCP discovery failed: {type(exc).__name__}: {exc}"
        ) from exc

    preferred = tuple(preferred_tool_names)
    selected = next(
        (
            candidate
            for candidate in tools
            if candidate.name in preferred
            or any(candidate.name.endswith(f"_{name}") for name in preferred)
        ),
        None,
    )
    if selected is None:
        available = ", ".join(sorted(tool.name for tool in tools)) or "<none>"
        raise RuntimeMCPError(
            f"{server_name} MCP does not expose any of {preferred}; available: {available}"
        )

    try:
        return await selected.ainvoke(arguments)
    except BaseException as exc:
        raise RuntimeMCPError(
            f"{server_name}.{selected.name} failed: {type(exc).__name__}: {exc}"
        ) from exc
