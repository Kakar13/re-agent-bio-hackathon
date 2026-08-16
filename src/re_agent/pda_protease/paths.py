"""Shared filesystem layout.

Everything cacheable lands under ``data/`` so a re-run costs nothing, and every
run-scoped artifact lands under ``results/pda_protease/<run-id>/``. Both trees
are gitignored except for their ``.gitkeep``.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]

DATA_RAW = REPO_ROOT / "data" / "raw"
DATA_PROCESSED = REPO_ROOT / "data" / "processed"

PDA_CACHE = DATA_RAW / "pda"
STRUCTURE_CACHE = DATA_RAW / "structures"
MEROPS_CACHE = DATA_RAW / "merops"
IEDB_CACHE = DATA_RAW / "iedb"
UNIPROT_CACHE = DATA_RAW / "uniprot"
BOLTZ_CACHE = DATA_RAW / "boltz"

RESULTS_ROOT = REPO_ROOT / "results" / "pda_protease"


def run_dir(run_id: str) -> Path:
    d = RESULTS_ROOT / run_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def ensure_caches() -> None:
    for d in (
        PDA_CACHE,
        STRUCTURE_CACHE,
        MEROPS_CACHE,
        IEDB_CACHE,
        UNIPROT_CACHE,
        BOLTZ_CACHE,
        DATA_PROCESSED,
        RESULTS_ROOT,
    ):
        d.mkdir(parents=True, exist_ok=True)


def env_flag(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}
