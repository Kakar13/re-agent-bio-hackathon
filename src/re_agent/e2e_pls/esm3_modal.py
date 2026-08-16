"""ESM3-open runtime: Modal GPU deployment, pooling recipes, and a
deterministic local embedding cache.

Two client modes, same interface:
  - `mode="modal"`: calls the deployed Modal app (`e2e-pls-esm3`), which
    runs `esm3-sm-open-v1` in half precision on an A10G. Requires the
    user's own `HF_TOKEN` (gated weights) and Modal auth -- see
    START_HERE.md / docs/SETUP.md. Not runnable in this session; nothing
    here invents or reads a Modal/HF credential.
  - `mode="mock"`: deterministic pseudo-embeddings/generations seeded from
    the input sequence. No GPU, no network. This is what fixtures, tests,
    and local dashboard dev run against; every other Track 2 module is
    built and verified against this path.

The Modal-side call surface (`LogitsConfig(return_embeddings=True)`,
masked `generate`) follows the public ESM3 SDK as of this plan's writing;
confirm against the installed `esm` package's API before the first real
`modal deploy`, since this path has not been execution-tested here.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Literal

import numpy as np

from re_agent.e2e_pls.schema import CANONICAL_RESIDUES, EMBEDDING_DIM, ESM3_MODEL_ID

# Human-facing model id (schema.py, dataset rows) -> real pretrained registry name
# used by the `esm` SDK's `ESM3.from_pretrained(...)`.
_PRETRAINED_NAME = {"esm3-sm-open-v1": "esm3_sm_open_v1"}

_RESIDUES = sorted(CANONICAL_RESIDUES)
MASK_TOKEN = "_"

ClientMode = Literal["modal", "mock"]


# --------------------------------------------------------------------------
# Deterministic local cache
# --------------------------------------------------------------------------


class EmbeddingCache:
    """Memory-mapped float32 cache, keyed by an arbitrary string key.

    Backs both Track 2's call-level cache (avoid re-hitting Modal for a
    sequence already embedded) and the row-id-keyed embedding store the
    plan calls for alongside the Track 1 parquet -- same access pattern,
    `get_or_compute(key, fn)` either way.
    """

    def __init__(self, path: str | Path, dim: int = EMBEDDING_DIM, initial_capacity: int = 1024):
        self.path = Path(path)
        self.dim = dim
        self.index_path = self.path.with_suffix(".index.json")
        self.path.parent.mkdir(parents=True, exist_ok=True)

        if self.index_path.exists():
            index = json.loads(self.index_path.read_text())
            self._key_to_row: dict[str, int] = index["key_to_row"]
            self._capacity = index["capacity"]
        else:
            self._key_to_row = {}
            self._capacity = initial_capacity
            self._write_index()

        self._mmap = self._open_mmap(self._capacity)

    def _open_mmap(self, capacity: int) -> np.memmap:
        mode = "r+" if self.path.exists() and self.path.stat().st_size > 0 else "w+"
        return np.memmap(self.path, dtype="float32", mode=mode, shape=(capacity, self.dim))

    def _write_index(self) -> None:
        self.index_path.write_text(
            json.dumps(
                {"dim": self.dim, "capacity": self._capacity, "key_to_row": self._key_to_row}
            )
        )

    def _grow(self, new_capacity: int) -> None:
        old_rows = len(self._key_to_row)
        old_data = np.array(self._mmap[:old_rows])
        del self._mmap
        new_mmap = np.memmap(self.path, dtype="float32", mode="w+", shape=(new_capacity, self.dim))
        new_mmap[:old_rows] = old_data
        new_mmap.flush()
        self._mmap = new_mmap
        self._capacity = new_capacity

    def __contains__(self, key: str) -> bool:
        return key in self._key_to_row

    def get(self, key: str) -> np.ndarray | None:
        row = self._key_to_row.get(key)
        if row is None:
            return None
        return np.array(self._mmap[row])

    def put(self, key: str, vector: np.ndarray) -> None:
        if vector.shape != (self.dim,):
            raise ValueError(f"expected shape ({self.dim},), got {vector.shape}")
        if key not in self._key_to_row:
            row = len(self._key_to_row)
            if row >= self._capacity:
                self._grow(self._capacity * 2)
            self._key_to_row[key] = row
            self._write_index()
        self._mmap[self._key_to_row[key]] = vector.astype("float32")
        self._mmap.flush()

    def get_or_compute(self, key: str, compute_fn) -> np.ndarray:
        cached = self.get(key)
        if cached is not None:
            return cached
        vector = compute_fn()
        self.put(key, vector)
        return vector


# --------------------------------------------------------------------------
# Pooling recipes (operate on a raw (L, EMBEDDING_DIM) residue embedding
# for the full n_flank + peptide + c_flank window)
# --------------------------------------------------------------------------


def pool_cleavage_termini(
    residue_embeddings: np.ndarray, n_flank_len: int, peptide_len: int, window: int = 3
) -> tuple[np.ndarray, np.ndarray]:
    """Mean-pool a small window straddling each cleavage site.

    N-terminus site sits between n_flank and the peptide's first residue;
    C-terminus site sits between the peptide's last residue and c_flank.
    """
    length = residue_embeddings.shape[0]
    n_lo, n_hi = max(0, n_flank_len - window), min(length, n_flank_len + window)
    c_site = n_flank_len + peptide_len
    c_lo, c_hi = max(0, c_site - window), min(length, c_site + window)
    return residue_embeddings[n_lo:n_hi].mean(axis=0), residue_embeddings[c_lo:c_hi].mean(axis=0)


def pool_9mer(residue_embeddings: np.ndarray, n_flank_len: int, peptide_len: int) -> np.ndarray:
    """Mean-pool just the peptide residues (excludes flanks)."""
    return residue_embeddings[n_flank_len : n_flank_len + peptide_len].mean(axis=0)


# --------------------------------------------------------------------------
# Client
# --------------------------------------------------------------------------


class ESM3Client:
    def __init__(
        self,
        mode: ClientMode = "mock",
        model_id: str = ESM3_MODEL_ID,
        cache: EmbeddingCache | None = None,
        modal_app_name: str = "e2e-pls-esm3",
    ):
        self.mode = mode
        self.model_id = model_id
        self.cache = cache
        self._modal_fns = None
        if mode == "modal":
            self._modal_fns = _load_modal_functions(modal_app_name)

    def embed_residues(self, sequence: str) -> np.ndarray:
        """Per-residue embeddings, shape (len(sequence), EMBEDDING_DIM)."""
        if not sequence:
            return np.zeros((0, EMBEDDING_DIM), dtype="float32")
        # Per-residue arrays vary in length across sequences, so they don't
        # fit the fixed-row EmbeddingCache; only pooled (fixed-dim) vectors
        # go through `self.cache` (see `embed_pooled`).
        if self.mode == "mock":
            return self._embed_residues_mock(sequence)
        return self._embed_residues_modal(sequence)

    def embed_pooled(self, peptide: str, n_flank: str, c_flank: str, recipe: str) -> np.ndarray:
        """Pooled, fixed-`EMBEDDING_DIM` vector for one of the head recipes.

        `recipe` in {"cleave_n", "cleave_c", "mean_9mer"}.
        """
        full_seq = f"{n_flank}{peptide}{c_flank}"
        cache_key = f"pooled|{recipe}|{self.model_id}|{full_seq}"

        def _compute() -> np.ndarray:
            residues = self.embed_residues(full_seq)
            if recipe == "mean_9mer":
                return pool_9mer(residues, len(n_flank), len(peptide))
            n_vec, c_vec = pool_cleavage_termini(residues, len(n_flank), len(peptide))
            return n_vec if recipe == "cleave_n" else c_vec

        if self.cache is None:
            return _compute()
        return self.cache.get_or_compute(cache_key, _compute)

    def generate_masked(
        self, sequence: str, mask_positions: list[int], num_candidates: int = 8
    ) -> list[str]:
        """Fill `mask_positions` via ESM3 sequence-track masked generation.

        Returns up to `num_candidates` full-length candidate sequences.
        Never mutates positions outside `mask_positions`.
        """
        if self.mode == "mock":
            return self._generate_masked_mock(sequence, mask_positions, num_candidates)
        return self._generate_masked_modal(sequence, mask_positions, num_candidates)

    def sequence_log_likelihood(self, sequence: str) -> float:
        """Scalar ESM3 sequence-track likelihood, used as a steering guardrail."""
        if self.mode == "mock":
            return self._sequence_log_likelihood_mock(sequence)
        return self._modal_fns["log_likelihood"].remote(sequence, self.model_id)

    # -- mock implementations --------------------------------------------

    def _seeded_rng(self, *parts: str) -> np.random.Generator:
        digest = hashlib.sha256("|".join(parts).encode()).digest()
        seed = int.from_bytes(digest[:8], "little")
        return np.random.default_rng(seed)

    def _embed_residues_mock(self, sequence: str) -> np.ndarray:
        rng = self._seeded_rng("residues", self.model_id, sequence)
        vecs = rng.standard_normal((len(sequence), EMBEDDING_DIM)).astype("float32")
        # per-residue identity bias so different residues get distinguishable,
        # reproducible directions instead of pure noise
        for i, aa in enumerate(sequence):
            aa_rng = self._seeded_rng("residue-identity", self.model_id, aa)
            vecs[i] += 0.5 * aa_rng.standard_normal(EMBEDDING_DIM).astype("float32")
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        return vecs / np.clip(norms, 1e-6, None)

    def _generate_masked_mock(
        self, sequence: str, mask_positions: list[int], num_candidates: int
    ) -> list[str]:
        candidates = []
        for c in range(num_candidates):
            rng = self._seeded_rng("mask-gen", self.model_id, sequence, str(mask_positions), str(c))
            chars = list(sequence)
            for pos in mask_positions:
                choices = [r for r in _RESIDUES if r != sequence[pos]]
                chars[pos] = choices[rng.integers(0, len(choices))]
            candidates.append("".join(chars))
        return candidates

    def _sequence_log_likelihood_mock(self, sequence: str) -> float:
        rng = self._seeded_rng("likelihood", self.model_id, sequence)
        # centered around a plausible per-residue log-prob; deterministic per sequence
        return float(rng.normal(-2.5, 0.3) * len(sequence))

    # -- real Modal implementations ---------------------------------------

    def _embed_residues_modal(self, sequence: str) -> np.ndarray:
        result = self._modal_fns["embed_residues"].remote([sequence], self.model_id)
        return np.array(result[0], dtype="float32")

    def _generate_masked_modal(
        self, sequence: str, mask_positions: list[int], num_candidates: int
    ) -> list[str]:
        masked = list(sequence)
        for pos in mask_positions:
            masked[pos] = MASK_TOKEN
        masked_seq = "".join(masked)
        return self._modal_fns["generate_masked"].remote(masked_seq, num_candidates, self.model_id)


def _load_modal_functions(app_name: str) -> dict:
    """Lazily resolve the deployed Modal functions. Requires `modal` to be
    installed and `modal token set` / `modal setup` to have been run by the
    user -- see docs/SETUP.md. Raises with a clear message otherwise.
    """
    try:
        import modal
    except ImportError as exc:
        raise RuntimeError(
            "mode='modal' requires the `modal` package: uv sync --extra e2e_pls (or --extra proto)"
        ) from exc
    return {
        "embed_residues": modal.Function.from_name(app_name, "embed_residues"),
        "generate_masked": modal.Function.from_name(app_name, "generate_masked"),
        "log_likelihood": modal.Function.from_name(app_name, "log_likelihood"),
    }


# --------------------------------------------------------------------------
# Modal app definition (deploy with: uv run modal deploy src/re_agent/e2e_pls/esm3_modal.py)
# --------------------------------------------------------------------------


def _build_modal_app():
    import modal

    image = (
        modal.Image.debian_slim(python_version="3.11")
        .pip_install("esm>=3.1.0", "torch>=2.2", "huggingface_hub>=0.24")
        .env({"HF_HUB_ENABLE_HF_TRANSFER": "0"})
    )
    app = modal.App("e2e-pls-esm3")
    weights_volume = modal.Volume.from_name("e2e-pls-esm3-weights", create_if_missing=True)
    volumes = {"/root/.cache/huggingface": weights_volume}
    # The user's own HF_TOKEN (gated ESM3 weights) must exist as a Modal
    # secret named "huggingface-secret" -- created via `modal secret create`,
    # never inlined here. See docs/SETUP.md.
    secrets = [modal.Secret.from_name("huggingface-secret")]

    def _load_client(model_id: str):
        from esm.models.esm3 import ESM3

        pretrained_name = _PRETRAINED_NAME.get(model_id, model_id)
        client = ESM3.from_pretrained(pretrained_name).to("cuda").half()
        return client

    @app.function(image=image, gpu="A10G", volumes=volumes, secrets=secrets, timeout=900)
    def embed_residues(
        sequences: list[str], model_id: str = ESM3_MODEL_ID
    ) -> list[list[list[float]]]:
        from esm.sdk.api import ESMProtein, LogitsConfig

        client = _load_client(model_id)
        # length-bucketed batching: group same-length sequences to avoid
        # padding waste, run each bucket together
        by_length: dict[int, list[tuple[int, str]]] = {}
        for i, seq in enumerate(sequences):
            by_length.setdefault(len(seq), []).append((i, seq))

        out: list[list[list[float]] | None] = [None] * len(sequences)
        for _length, items in by_length.items():
            for i, seq in items:
                protein = ESMProtein(sequence=seq)
                protein_tensor = client.encode(protein)
                logits_output = client.logits(protein_tensor, LogitsConfig(return_embeddings=True))
                embeddings = logits_output.embeddings.squeeze(0).float().cpu().numpy()
                out[i] = embeddings.tolist()
        return out

    @app.function(image=image, gpu="A10G", volumes=volumes, secrets=secrets, timeout=900)
    def generate_masked(
        masked_sequence: str, num_candidates: int, model_id: str = ESM3_MODEL_ID
    ) -> list[str]:
        from esm.sdk.api import ESMProtein, GenerationConfig

        client = _load_client(model_id)
        candidates = []
        for _ in range(num_candidates):
            protein = ESMProtein(sequence=masked_sequence)
            generated = client.generate(
                protein, GenerationConfig(track="sequence", num_steps=8, temperature=0.7)
            )
            candidates.append(generated.sequence)
        return candidates

    @app.function(image=image, gpu="A10G", volumes=volumes, secrets=secrets, timeout=300)
    def log_likelihood(sequence: str, model_id: str = ESM3_MODEL_ID) -> float:
        from esm.sdk.api import ESMProtein, LogitsConfig

        client = _load_client(model_id)
        protein = ESMProtein(sequence=sequence)
        protein_tensor = client.encode(protein)
        logits_output = client.logits(protein_tensor, LogitsConfig(sequence=True))
        log_probs = logits_output.logits.sequence.log_softmax(-1)
        token_ids = protein_tensor.sequence
        return float(log_probs.gather(-1, token_ids.unsqueeze(-1)).sum().item())

    return app


if __name__ == "__main__":
    app = _build_modal_app()
