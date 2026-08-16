from __future__ import annotations

import json

from re_agent.agent.mcp_runtime import runtime_mcp_connections, runtime_mcp_status


def test_runtime_mcp_status_never_exposes_credentials(monkeypatch) -> None:
    monkeypatch.delenv("PROTO_TOOLS_REPO", raising=False)
    monkeypatch.setenv("PROTO_API_KEY", "proto-secret")
    monkeypatch.setenv("PAPERCLIP_MCP_BEARER_TOKEN", "paperclip-secret")

    connections = runtime_mcp_connections()
    status = runtime_mcp_status()

    assert connections["proto"]["transport"] == "http"
    assert connections["paperclip"]["transport"] == "http"
    rendered = json.dumps(status)
    assert "proto-secret" not in rendered
    assert "paperclip-secret" not in rendered
    assert status["paperclip"]["authentication"] == "bearer"
    assert status["proto"]["authentication"] == "bearer"


def test_runtime_mcp_prefers_local_proto_monorepo(monkeypatch, tmp_path) -> None:
    (tmp_path / "pyproject.toml").write_text("[project]\nname='proto-tools'\nversion='0'\n")
    monkeypatch.setenv("PROTO_TOOLS_REPO", str(tmp_path))
    monkeypatch.setenv("PROTO_API_KEY", "hosted-fallback")

    connections = runtime_mcp_connections()
    status = runtime_mcp_status()

    assert connections["proto"]["transport"] == "stdio"
    assert connections["proto"]["args"] == [
        "run",
        "--directory",
        str(tmp_path),
        "--extra",
        "mcp",
        "proto-tools-mcp",
        "--device",
        "modal",
    ]
    assert status["proto"]["authentication"] == "local_modal"
    assert status["proto"]["source"] == "monorepo"
