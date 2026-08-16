"""Frozen ESM-2 per-residue embeddings with an on-disk cache.

The encoder is never trained: embeddings are computed once and reused by every
experiment, which is what makes repeated Mean Teacher runs cheap on a laptop.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.immuno.config import (
    ESM_DIM,
    ESM_LAYER,
    ESM_MODEL,
    PATHS,
    WINDOW,
    ensure_dirs,
    torch_device,
)

_MODEL_CACHE: dict[str, tuple] = {}


def load_esm(device: str | None = None):
    """Load ESM-2 once per process and keep it in eval mode."""
    import esm

    device = device or torch_device()
    key = f"{ESM_MODEL}:{device}"
    if key in _MODEL_CACHE:
        return _MODEL_CACHE[key]

    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM_MODEL)
    model = model.to(device).eval()
    for p in model.parameters():
        p.requires_grad_(False)
    converter = alphabet.get_batch_converter()
    _MODEL_CACHE[key] = (model, alphabet, converter, device)
    return _MODEL_CACHE[key]


def embed_sequences(
    seqs: list[str], batch_size: int = 512, device: str | None = None, progress: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """Embed isolated windows.

    Returns `(n, WINDOW, ESM_DIM)` float16 embeddings zero-padded on the right and
    an `(n, WINDOW)` float32 mask marking real residues.
    """
    import torch

    model, alphabet, converter, device = load_esm(device)
    n = len(seqs)
    out = np.zeros((n, WINDOW, ESM_DIM), dtype=np.float16)
    mask = np.zeros((n, WINDOW), dtype=np.float32)

    start = time.time()
    for i in range(0, n, batch_size):
        chunk = seqs[i : i + batch_size]
        _, _, tokens = converter([(str(j), s) for j, s in enumerate(chunk)])
        tokens = tokens.to(device)
        with torch.no_grad():
            rep = model(tokens, repr_layers=[ESM_LAYER])["representations"][ESM_LAYER]
        rep = rep.float().cpu().numpy()
        for j, s in enumerate(chunk):
            length = min(len(s), WINDOW)
            # Strip the BOS token; ESM appends EOS after the sequence.
            out[i + j, :length] = rep[j, 1 : length + 1].astype(np.float16)
            mask[i + j, :length] = 1.0
        if progress and (i // batch_size) % 20 == 0:
            done = min(i + batch_size, n)
            rate = done / max(time.time() - start, 1e-6)
            print(f"  embedded {done}/{n} ({rate:.0f}/s)", flush=True)
    return out, mask


def cache_paths(name: str) -> tuple[Path, Path, Path]:
    base = PATHS.embed_cache / f"{name}__{ESM_MODEL}"
    return (
        base.with_suffix(".emb.npy"),
        base.with_suffix(".mask.npy"),
        base.with_suffix(".index.parquet"),
    )


def embed_table(df: pd.DataFrame, name: str, force: bool = False, batch_size: int = 512) -> Path:
    """Embed a window table, writing arrays aligned to `df` row order."""
    ensure_dirs()
    emb_path, mask_path, idx_path = cache_paths(name)
    if emb_path.exists() and not force:
        cached = pd.read_parquet(idx_path)
        if len(cached) == len(df) and (cached["seq"].to_numpy() == df["seq"].to_numpy()).all():
            print(f"embeddings cached: {name} ({len(df)})")
            return emb_path

    seqs = df["seq"].tolist()
    print(f"embedding {name}: {len(seqs)} windows on {torch_device()}")
    emb, mask = embed_sequences(seqs, batch_size=batch_size, progress=True)
    np.save(emb_path, emb)
    np.save(mask_path, mask)
    df[["seq"]].to_parquet(idx_path, index=False)
    print(f"  wrote {emb_path.name} {emb.shape}")
    return emb_path


def load_embeddings(name: str) -> tuple[np.ndarray, np.ndarray]:
    emb_path, mask_path, _ = cache_paths(name)
    if not emb_path.exists():
        raise FileNotFoundError(f"missing embedding cache for '{name}' — run embed.main() first")
    return np.load(emb_path, mmap_mode="r"), np.load(mask_path, mmap_mode="r")


def main() -> None:
    from re_agent.immuno import data

    for name, loader in (
        ("labeled", data.build_labeled),
        ("unlabeled", data.build_unlabeled),
        ("reference", data.build_reference),
    ):
        embed_table(loader(), name)


if __name__ == "__main__":
    main()
