#!/usr/bin/env python3
"""Build the NetMHCpan corpus as a detachable Modal CPU job.

Example:

    uv run --extra proto modal run --detach scripts/build_netmhcpan_corpus_modal.py \
      --run-name a0201-netmhcpan-pda-full-v1 \
      --pda-only-full \
      --parent-batch-size 20
"""

from __future__ import annotations

import json
from pathlib import Path

import modal

APP_NAME = "mhci-netmhcpan-corpus"
ARTIFACT_VOLUME_NAME = "mhci-netmhcpan-artifacts"
CACHE_VOLUME_NAME = "mhci-netmhcpan-cache"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "httpx>=0.28",
        "numpy>=1.26",
        "pandas>=2.2",
        "pyarrow>=17",
    )
    .add_local_python_source("re_agent")
)

app = modal.App(APP_NAME)
artifact_volume = modal.Volume.from_name(ARTIFACT_VOLUME_NAME, create_if_missing=True)
cache_volume = modal.Volume.from_name(CACHE_VOLUME_NAME, create_if_missing=True)


@app.function(
    image=image,
    cpu=4,
    memory=16384,
    timeout=6 * 60 * 60,
    volumes={"/artifacts": artifact_volume, "/cache": cache_volume},
)
def build_corpus(
    *,
    run_name: str,
    pda_designs_bytes: bytes | None,
    target_rows: int,
    pda_challenge_rows: int,
    n_human: int,
    n_viral: int,
    n_bacterial: int,
    n_de_novo: int,
    api_batch_size: int,
    pda_only_full: bool,
    parent_batch_size: int,
    seed: int,
) -> dict:
    from re_agent.e2e_pls.netmhcpan_corpus import (
        build_corpus_artifacts,
        build_pda_corpus_artifacts,
    )

    pda_path = Path("/tmp/pda_designs.parquet")
    if pda_designs_bytes is not None:
        pda_path.write_bytes(pda_designs_bytes)
    elif pda_challenge_rows or pda_only_full:
        raise ValueError("PDA rows requested but no PDA designs parquet was uploaded")

    output_dir = Path("/artifacts/corpora") / run_name
    if pda_only_full:
        manifest = build_pda_corpus_artifacts(
            output_dir=output_dir,
            cache_dir=Path("/cache"),
            pda_designs_path=pda_path,
            api_batch_size=api_batch_size,
            parent_batch_size=parent_batch_size,
            seed=seed,
        )
    else:
        manifest = build_corpus_artifacts(
            output_dir=output_dir,
            cache_dir=Path("/cache"),
            target_rows=target_rows,
            pda_challenge_rows=pda_challenge_rows,
            n_human_proteins=n_human,
            n_viral_proteins=n_viral,
            n_bacterial_proteins=n_bacterial,
            n_de_novo_proteins=n_de_novo,
            pda_designs_path=pda_path,
            api_batch_size=api_batch_size,
            seed=seed,
        )
    cache_volume.commit()
    artifact_volume.commit()
    manifest["modal"] = {
        "artifact_volume": ARTIFACT_VOLUME_NAME,
        "artifact_path": f"/corpora/{run_name}",
        "cache_volume": CACHE_VOLUME_NAME,
    }
    return manifest


@app.local_entrypoint()
def main(
    run_name: str = "a0201-netmhcpan-10k",
    pda_designs: str = "data/processed/pda_designs.parquet",
    target_rows: int = 10_000,
    pda_challenge_rows: int = 1_000,
    n_human: int = 120,
    n_viral: int = 60,
    n_bacterial: int = 80,
    n_de_novo: int = 40,
    api_batch_size: int = 500,
    pda_only_full: bool = False,
    parent_batch_size: int = 20,
    seed: int = 0,
) -> None:
    pda_path = Path(pda_designs)
    pda_bytes = None
    if pda_challenge_rows or pda_only_full:
        if not pda_path.exists():
            raise SystemExit(f"PDA designs not found: {pda_path}")
        pda_bytes = pda_path.read_bytes()
    result = build_corpus.remote(
        run_name=run_name,
        pda_designs_bytes=pda_bytes,
        target_rows=target_rows,
        pda_challenge_rows=pda_challenge_rows,
        n_human=n_human,
        n_viral=n_viral,
        n_bacterial=n_bacterial,
        n_de_novo=n_de_novo,
        api_batch_size=api_batch_size,
        pda_only_full=pda_only_full,
        parent_batch_size=parent_batch_size,
        seed=seed,
    )
    print(json.dumps(result, indent=2))
    print(
        "\nInspect cloud artifacts with:\n"
        f"  uv run --extra proto modal volume ls {ARTIFACT_VOLUME_NAME} "
        f"corpora/{run_name}"
    )
