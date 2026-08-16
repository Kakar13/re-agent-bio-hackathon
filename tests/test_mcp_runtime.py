from __future__ import annotations

import json

from re_agent.agent.mcp_runtime import runtime_mcp_connections, runtime_mcp_status


def test_runtime_mcp_status_never_exposes_credentials(monkeypatch) -> None:
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
