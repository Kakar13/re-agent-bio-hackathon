"""Versioned local run manifests retained alongside LangSmith traces."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[3]
RUNS_DIR = ROOT / "results" / "workbench" / "runs"
SCHEMA_VERSION = "1.0"


def _git_revision() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def new_run_id(prefix: str) -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return f"{prefix}-{stamp}-{uuid4().hex[:8]}"


def write_run_artifact(
    *,
    run_id: str,
    kind: str,
    payload: dict[str, Any],
    filename: str = "artifact.json",
    parent_run_id: str | None = None,
) -> dict[str, Any]:
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = run_dir / filename
    artifact_path.write_text(json.dumps(payload, indent=2, default=str) + "\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "parent_run_id": parent_run_id,
        "kind": kind,
        "created_at": datetime.now(UTC).isoformat(),
        "git_revision": _git_revision(),
        "artifact": {
            "path": str(artifact_path.relative_to(ROOT)),
            "sha256": _sha256(artifact_path),
        },
        "observability": {
            "provider": "langsmith",
            "project": os.getenv("LANGSMITH_PROJECT", "re-agent-scientific-workbench"),
            "tracing_enabled": os.getenv("LANGSMITH_TRACING", "").lower() == "true",
        },
    }
    manifest_path = run_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return {
        "id": run_id,
        "kind": kind,
        "title": payload.get("title", run_id),
        "path": manifest["artifact"]["path"],
        "manifest_path": str(manifest_path.relative_to(ROOT)),
        "sha256": manifest["artifact"]["sha256"],
        "payload": payload,
    }


def list_runs(limit: int = 50) -> list[dict[str, Any]]:
    if not RUNS_DIR.exists():
        return []
    manifests: list[dict[str, Any]] = []
    for path in sorted(RUNS_DIR.glob("*/manifest.json"), reverse=True):
        try:
            manifest = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        manifest["manifest_path"] = str(path.relative_to(ROOT))
        manifests.append(manifest)
        if len(manifests) >= limit:
            break
    return manifests
