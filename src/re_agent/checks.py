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


def _proto_sdk() -> Check:
    try:
        import proto_client  # noqa: F401

        return Check("proto-client", True, "import ok")
    except ImportError:
        return Check("proto-client", False, "run: uv pip install proto-client")


def run_checks() -> list[Check]:
    return [
        _which("python3"),
        _which("uv"),
        _which("claude"),
        _paperclip(),
        _proto_sdk(),
        Check(
            "ANTHROPIC_API_KEY",
            bool(settings.anthropic_api_key),
            "set in .env" if settings.anthropic_api_key else "missing — grab credits at lightning talks",
        ),
        Check(
            "PROTO_API_KEY",
            bool(settings.proto_api_key),
            "set in .env" if settings.proto_api_key else "missing — Proto workspace key",
        ),
    ]
