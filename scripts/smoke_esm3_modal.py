#!/usr/bin/env python3
"""Run ESM3-open on Modal and print live progress + embedding stats."""

from __future__ import annotations

import sys
import time

from proto_tools.tools.masked_models.esm3.esm3_embeddings import (
    ESM3EmbeddingsConfig,
    run_esm3_embeddings,
)
from proto_tools.tools.masked_models.shared_data_models import MaskedModelInput

DEFAULT_SEQ = "MKTL"


def log(msg: str) -> None:
    print(msg, flush=True)


def main() -> int:
    sequence = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SEQ
    log("ESM3 Modal smoke")
    log("-" * 48)
    log(f"sequence : {sequence}")
    log(f"length   : {len(sequence)}")
    log("device   : modal")
    log("model    : esm3_sm_open_v1")
    log("-" * 48)
    log("Connecting to Modal and loading ESM3 (first call is slow)...")

    started = time.perf_counter()
    out = run_esm3_embeddings(
        MaskedModelInput(sequences=[sequence]),
        ESM3EmbeddingsConfig(device="modal", verbose=1),
    )
    elapsed = time.perf_counter() - started

    log("-" * 48)
    log(f"success         : {out.success}")
    log(f"errors          : {out.errors or '[]'}")
    log(f"tool_id         : {out.tool_id}")
    log(f"wall_time_s     : {elapsed:.1f}")
    log(f"tool_time_s     : {getattr(out, 'execution_time', None)}")

    if not out.success or not out.results:
        log("No embedding returned.")
        return 1

    result = out.results[0]
    embedding = result.mean_embedding
    log(f"embedding_dim   : {len(embedding)}")
    log(f"first_8         : {[round(x, 4) for x in embedding[:8]]}")
    log(f"mean            : {sum(embedding) / len(embedding):.6f}")
    log(f"l2_norm         : {sum(x * x for x in embedding) ** 0.5:.4f}")
    log("done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
