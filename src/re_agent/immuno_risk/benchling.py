"""Benchling handoff: pull AA sequences, publish idempotent run summaries."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)


def _client():
    url = os.environ.get("BENCHLING_TENANT_URL") or os.environ.get("BENCHLING_URL")
    api_key = os.environ.get("BENCHLING_API_KEY")
    if not url or not api_key:
        raise RuntimeError(
            "Benchling not configured. Set BENCHLING_TENANT_URL and BENCHLING_API_KEY "
            "(from Benchling booth / Discord). Never commit credentials."
        )
    try:
        from benchling_sdk.auth.api_key_auth import ApiKeyAuth
        from benchling_sdk.benchling import Benchling
    except ImportError as exc:
        raise RuntimeError("Install benchling-sdk: uv sync --extra immuno") from exc
    return Benchling(url=url.rstrip("/"), auth_method=ApiKeyAuth(api_key))


def configured() -> bool:
    return bool(
        (os.environ.get("BENCHLING_TENANT_URL") or os.environ.get("BENCHLING_URL"))
        and os.environ.get("BENCHLING_API_KEY")
    )


def pull_candidates(
    *,
    ids: list[str] | None = None,
    name_includes: str | None = None,
    dry_run: bool = False,
    limit: int = 20,
) -> list[dict[str, Any]]:
    """Pull AA Sequence entities as pipeline candidates."""
    if dry_run or not configured():
        return [
            {
                "id": "aa_dry_run_example",
                "name": "DRY_RUN_CANDIDATE",
                "amino_acids": "GILGFVFTLAAAKKKLLLGGGGSSSTTTT",
                "dry_run": True,
                "note": "Configure BENCHLING_* env vars for live pull",
            }
        ]
    benchling = _client()
    folder_id = os.environ.get("BENCHLING_FOLDER_ID")
    out: list[dict[str, Any]] = []
    if ids:
        for aa_id in ids:
            seq = benchling.aa_sequences.get_by_id(aa_id)
            out.append(
                {
                    "id": seq.id,
                    "name": seq.name,
                    "amino_acids": getattr(seq, "amino_acids", None) or getattr(seq, "bases", ""),
                    "web_url": getattr(seq, "web_url", None),
                }
            )
        return out

    kwargs: dict[str, Any] = {}
    if folder_id:
        kwargs["folder_id"] = folder_id
    if name_includes:
        kwargs["name_includes"] = name_includes
    pages = benchling.aa_sequences.list(**kwargs)
    for page in pages:
        for seq in page:
            out.append(
                {
                    "id": seq.id,
                    "name": seq.name,
                    "amino_acids": getattr(seq, "amino_acids", None) or getattr(seq, "bases", ""),
                    "web_url": getattr(seq, "web_url", None),
                }
            )
            if len(out) >= limit:
                return out
    return out


def _run_payload(run_dir: Path) -> dict[str, Any]:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    summary = json.loads((run_dir / "summary.json").read_text())
    return {
        "run_id": manifest["run_id"],
        "sequence_id": manifest["sequence_id"],
        "delivery_mode": manifest.get("delivery_mode"),
        "overall_risk": manifest.get("overall_risk"),
        "score0to100": manifest.get("score0to100"),
        "confidence": manifest.get("confidence"),
        "aggregation_overall": manifest.get("aggregation_overall"),
        "predictor_versions": manifest.get("predictor_versions"),
        "checksums_sha256": manifest.get("checksums_sha256"),
        "caveats": manifest.get("caveats"),
        "top_flagged": summary.get("risk", {}).get("peptides_flagged", [])[:15],
        "local_artifact_dir": str(run_dir),
        "trace_url": os.environ.get("LANGSMITH_TRACE_URL") or "",
    }


def publish_run(run_dir: Path, *, dry_run: bool = False) -> dict[str, Any]:
    """Publish idempotent run record linked to source sequence."""
    run_dir = Path(run_dir)
    payload = _run_payload(run_dir)
    run_id = payload["run_id"]

    if dry_run or not configured():
        return {
            "status": "dry_run",
            "run_id": run_id,
            "payload_preview": payload,
            "note": "Set BENCHLING_* and omit --dry-run to publish",
        }

    benchling = _client()
    folder_id = os.environ.get("BENCHLING_FOLDER_ID")
    if not folder_id:
        raise RuntimeError("BENCHLING_FOLDER_ID required to publish run records")

    # Idempotency: search custom entities / entries by name containing run_id
    name = f"immuno-risk-{run_id}"
    schema_id = os.environ.get("BENCHLING_RUN_SCHEMA_ID")

    # Prefer custom entity when schema present; else create a notebook entry
    try:
        if schema_id:
            from benchling_sdk.helpers.serialization_helpers import fields
            from benchling_sdk.models import CustomEntityCreate

            # Check existing by name
            existing = list(benchling.custom_entities.list(name_includes=run_id, folder_id=folder_id))
            for page in existing:
                for ent in page:
                    if run_id in (ent.name or ""):
                        return {
                            "status": "exists",
                            "run_id": run_id,
                            "id": ent.id,
                            "web_url": getattr(ent, "web_url", None),
                        }

            entity = CustomEntityCreate(
                name=name,
                folder_id=folder_id,
                schema_id=schema_id,
                fields=fields(
                    {
                        k: {"value": json.dumps(v) if isinstance(v, (dict, list)) else v}
                        for k, v in {
                            "Run ID": run_id,
                            "Overall Risk": payload["overall_risk"],
                            "Score": payload["score0to100"],
                            "Confidence": payload["confidence"],
                            "Aggregation": payload["aggregation_overall"],
                            "Artifact Dir": payload["local_artifact_dir"],
                        }.items()
                        if v is not None
                    }
                ),
            )
            created = benchling.custom_entities.create(entity)
            return {
                "status": "created",
                "run_id": run_id,
                "id": created.id,
                "web_url": getattr(created, "web_url", None),
            }

        from benchling_sdk.models import EntryCreate

        entry = EntryCreate(name=name, folder_id=folder_id)
        created = benchling.entries.create(entry)
        return {
            "status": "created_entry_reduced_mode",
            "run_id": run_id,
            "id": created.id,
            "web_url": getattr(created, "web_url", None),
            "note": "No BENCHLING_RUN_SCHEMA_ID — reduced notebook entry mode",
            "payload": payload,
        }
    except Exception as exc:  # noqa: BLE001
        log.exception("Benchling publish failed")
        return {"status": "error", "run_id": run_id, "error": str(exc)}
