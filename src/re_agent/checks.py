"""Local environment checks for the weekend build."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass

from re_agent.config import settings


@dataclass(frozen=True)
class Check:
    name: str
    ok: bool
    detail: str


def _which(cmd: str) -> Check:
    path = shutil.which(cmd)
    if path:
        return Check(cmd, True, path)
    return Check(cmd, False, "not on PATH")


def _paperclip() -> Check:
    if not shutil.which("paperclip"):
        return Check("paperclip", False, "run: curl -fsSL https://paperclip.gxl.ai/install.sh | bash")
    try:
        result = subprocess.run(
            ["paperclip", "config"],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
        output = (result.stdout or result.stderr).strip()
        ok = result.returncode == 0 and "Auth" in output
        return Check("paperclip", ok, output.splitlines()[0] if output else "installed")
    except (OSError, subprocess.TimeoutExpired) as exc:
        return Check("paperclip", False, str(exc))


def _import_pkg(name: str, module: str, hint: str) -> Check:
    try:
        __import__(module)
        return Check(name, True, "import ok")
    except ImportError:
        return Check(name, False, hint)


def _modal() -> Check:
    path = shutil.which("modal")
    if path:
        return Check("modal", True, path)
    # uv run puts tools on PATH for the venv; try import as fallback
    try:
        import modal  # noqa: F401

        return Check("modal", True, "import ok — run: uv run modal setup")
    except ImportError:
        return Check("modal", False, "run: uv sync --extra proto && uv run modal setup")


def run_checks() -> list[Check]:
    hint = "run: uv sync --extra proto"
    return [
        _which("python3"),
        _which("uv"),
        Check(
            "pi",
            bool(shutil.which("pi")),
            shutil.which("pi") or "run: curl -fsSL https://pi.dev/install.sh | sh",
        ),
        _which("claude"),
        _paperclip(),
        _import_pkg("proto-client", "proto_client", hint),
        _import_pkg("proto-tools", "proto_tools", hint),
        _import_pkg("proto-language", "proto_language", hint),
        _modal(),
        Check(
            "ANTHROPIC_API_KEY",
            bool(settings.anthropic_api_key),
            "set in .env" if settings.anthropic_api_key else "missing — grab credits at lightning talks",
        ),
        Check(
            "PROTO_API_KEY",
            bool(settings.proto_api_key),
            "set in .env" if settings.proto_api_key else "missing — Proto workspace key (hosted SDK / Cursor MCP)",
        ),
    ]
