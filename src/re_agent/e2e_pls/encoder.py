"""Protein sequence encoder: ESM-2 locally, plus a mock and optional Modal
ESM3 backend, all behind a single `ProteinEncoder` interface.

Three modes, same interface:
  - `mode="esm2"` (default, real): loads `esm2_t33_650M_UR50D` locally via
    the `fair-esm` package. First use downloads ~2.5 GB of weights from
    Meta's public URL (no HF token, no license acceptance). Runs on CUDA,
    Apple MPS, or CPU -- picked automatically. This is what training,
    inference, and the dashboard use for real embeddings.
  - `mode="mock"`: deterministic pseudo-embeddings/generations seeded from
    the input sequence. No download, no GPU. This is what tests and the
    dev fixture use to stay fast and offline.
  - `mode="modal"`: calls the deployed Modal app for ESM3-open weights.
    Requires the user's own `HF_TOKEN` (gated) and Modal auth. Never
    executed in this session; kept as an optional upgrade path.

`EMBEDDING_DIM` in `schema.py` (1280) is tied to ESM-2 t33's hidden size.
If you swap encoders, update it there and rebuild the dataset so
`encoder_model_id` in each row matches what actually produced the vectors.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Literal

import numpy as np

from re_agent.e2e_pls.schema import CANONICAL_RESIDUES, EMBEDDING_DIM, ENCODER_MODEL_ID

_RESIDUES = sorted(CANONICAL_RESIDUES)
MASK_TOKEN = "_"

ClientMode = Literal["esm2", "modal", "mock"]


# --------------------------------------------------------------------------
# Deterministic local cache
# --------------------------------------------------------------------------


class EmbeddingCache:
    """Memory-mapped float32 cache, keyed by an arbitrary string key.

    Backs both the encoder's call-level cache (avoid re-embedding a sequence
    already seen) and the row-id-keyed embedding store alongside the Track 1
    parquet -- same access pattern, `get_or_compute(key, fn)` either way.
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
# ESM-2 backend
# --------------------------------------------------------------------------

# Human-facing model id (schema.ENCODER_MODEL_ID) -> `fair-esm` registry name.
_ESM2_PRETRAINED = {
    "esm2_t33_650M_UR50D": "esm2_t33_650M_UR50D",
    "esm2_t30_150M_UR50D": "esm2_t30_150M_UR50D",
    "esm2_t12_35M_UR50D": "esm2_t12_35M_UR50D",
}
_ESM2_REPR_LAYER = {
    "esm2_t33_650M_UR50D": 33,
    "esm2_t30_150M_UR50D": 30,
    "esm2_t12_35M_UR50D": 12,
}


def _pick_torch_device():
    import torch

    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


class _Esm2Backend:
    """Lazy singleton wrapping fair-esm's ESM-2 model + alphabet + batch converter."""

    _cache: dict[str, _Esm2Backend] = {}

    def __init__(self, model_id: str):
        import esm as fair_esm
        import torch

        if model_id not in _ESM2_PRETRAINED:
            raise ValueError(f"unknown ESM-2 model_id: {model_id}. Known: {list(_ESM2_PRETRAINED)}")
        pretrained_name = _ESM2_PRETRAINED[model_id]
        loader = getattr(fair_esm.pretrained, pretrained_name)
        model, alphabet = loader()
        model.eval()
        self.device = _pick_torch_device()
        model = model.to(self.device)
        self.model = model
        self.alphabet = alphabet
        self.batch_converter = alphabet.get_batch_converter()
        self.repr_layer = _ESM2_REPR_LAYER[model_id]
        self._torch = torch

    @classmethod
    def get(cls, model_id: str) -> _Esm2Backend:
        if model_id not in cls._cache:
            cls._cache[model_id] = cls(model_id)
        return cls._cache[model_id]

    def embed_residues(self, sequence: str) -> np.ndarray:
        """Per-residue embeddings from the final transformer layer, shape (L, D)."""
        torch = self._torch
        _, _, tokens = self.batch_converter([("s", sequence)])
        tokens = tokens.to(self.device)
        with torch.no_grad():
            out = self.model(tokens, repr_layers=[self.repr_layer])
        # fair-esm prepends a BOS token and appends an EOS token, so we strip both.
        reps = out["representations"][self.repr_layer][0, 1 : 1 + len(sequence)]
        return reps.float().cpu().numpy()

    def fill_masks(
        self, sequence: str, mask_positions: list[int], num_candidates: int
    ) -> list[str]:
        """Sample from the ESM-2 masked-language-model head at each masked position.

        For each of `num_candidates` draws, replace `mask_positions` with the
        `<mask>` token, run the model, and sample from the softmax distribution
        at each position independently. Never mutates unmasked positions.
        """
        torch = self._torch
        rng = np.random.default_rng()
        masked = list(sequence)
        for pos in mask_positions:
            masked[pos] = "<mask>"
        # fair-esm's batch converter expects the sequence as a single string with
        # `<mask>` tokens spelled out; we tokenize ourselves for control.
        raw = "".join(c if c != "<mask>" else "<mask>" for c in masked)
        _, _, tokens = self.batch_converter([("s", raw)])
        tokens = tokens.to(self.device)
        with torch.no_grad():
            out = self.model(tokens)
        logits = out["logits"][0]  # (L, vocab)
        # positions in the token tensor are shifted by 1 for the prepended BOS
        token_positions = [p + 1 for p in mask_positions]

        # only sample among the 20 canonical amino acid tokens, not special tokens
        aa_indices = [self.alphabet.get_idx(aa) for aa in _RESIDUES]
        original_indices = [self.alphabet.get_idx(sequence[p]) for p in mask_positions]

        candidates = []
        for _ in range(num_candidates):
            chars = list(sequence)
            for pos, tok_pos, orig_idx in zip(
                mask_positions, token_positions, original_indices, strict=True
            ):
                pos_logits = logits[tok_pos, aa_indices].float().cpu().numpy()
                pos_logits[aa_indices.index(orig_idx)] = -np.inf  # force a mutation
                probs = np.exp(pos_logits - pos_logits.max())
                probs = probs / probs.sum()
                pick = rng.choice(len(aa_indices), p=probs)
                chars[pos] = _RESIDUES[pick]
            candidates.append("".join(chars))
        return candidates

    def sequence_log_likelihood(self, sequence: str) -> float:
        """Sum of log-probs of the observed residues under the LM head."""
        torch = self._torch
        _, _, tokens = self.batch_converter([("s", sequence)])
        tokens = tokens.to(self.device)
        with torch.no_grad():
            out = self.model(tokens)
        logits = out["logits"][0, 1 : 1 + len(sequence)]  # strip BOS/EOS
        log_probs = torch.log_softmax(logits.float(), dim=-1)
        aa_indices = torch.tensor(
            [self.alphabet.get_idx(aa) for aa in sequence], device=logits.device
        )
        return float(log_probs.gather(-1, aa_indices.unsqueeze(-1)).sum().cpu().item())


# --------------------------------------------------------------------------
# Client
# --------------------------------------------------------------------------


class ProteinEncoder:
    def __init__(
        self,
        mode: ClientMode = "esm2",
        model_id: str = ENCODER_MODEL_ID,
        cache: EmbeddingCache | None = None,
        modal_app_name: str = "e2e-pls-esm3",
    ):
        self.mode = mode
        self.model_id = model_id
        self.cache = cache
        self._modal_fns = None
        self._esm2 = None
        if mode == "modal":
            self._modal_fns = _load_modal_functions(modal_app_name)
        elif mode == "esm2":
            self._esm2 = _Esm2Backend.get(model_id)

    def embed_residues(self, sequence: str) -> np.ndarray:
        """Per-residue embeddings, shape (len(sequence), EMBEDDING_DIM)."""
        if not sequence:
            return np.zeros((0, EMBEDDING_DIM), dtype="float32")
        # Per-residue arrays vary in length across sequences, so they don't
        # fit the fixed-row EmbeddingCache; only pooled (fixed-dim) vectors
        # go through `self.cache` (see `embed_pooled`).
        if self.mode == "mock":
            return self._embed_residues_mock(sequence)
        if self.mode == "esm2":
            return self._esm2.embed_residues(sequence).astype("float32")
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
        """Fill `mask_positions` via the encoder's masked-language-model head.

        Returns up to `num_candidates` full-length candidate sequences.
        Never mutates positions outside `mask_positions`.
        """
        if self.mode == "mock":
            return self._generate_masked_mock(sequence, mask_positions, num_candidates)
        if self.mode == "esm2":
            return self._esm2.fill_masks(sequence, mask_positions, num_candidates)
        return self._generate_masked_modal(sequence, mask_positions, num_candidates)

    def sequence_log_likelihood(self, sequence: str) -> float:
        """Scalar language-model likelihood, used as a steering guardrail."""
        if self.mode == "mock":
            return self._sequence_log_likelihood_mock(sequence)
        if self.mode == "esm2":
            return self._esm2.sequence_log_likelihood(sequence)
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
            rng = self._seeded_rng(
                "mask-gen", self.model_id, sequence, str(mask_positions), str(c)
            )
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

    # -- real Modal (ESM3) implementations, unused unless mode="modal" ---

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
        return self._modal_fns["generate_masked"].remote(
            masked_seq, num_candidates, self.model_id
        )


def _load_modal_functions(app_name: str) -> dict:
    """Lazily resolve the deployed Modal ESM3 functions. Requires `modal` to
    be installed and `modal setup` to have been run by the user; the ESM3
    weights are gated so an HF_TOKEN secret must exist too. Not runnable in
    this repo out of the box -- see project README for the opt-in path.
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
# Modal app definition for the optional ESM3 backend.
# Deploy with:  uv run modal deploy src/re_agent/e2e_pls/encoder.py
# --------------------------------------------------------------------------

# Map schema-level model id -> real pretrained registry name used by the ESM SDK.
_ESM3_PRETRAINED = {"esm3-sm-open-v1": "esm3_sm_open_v1"}


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
    secrets = [modal.Secret.from_name("huggingface-secret")]

    def _load_client(model_id: str):
        from esm.models.esm3 import ESM3

        pretrained_name = _ESM3_PRETRAINED.get(model_id, model_id)
        client = ESM3.from_pretrained(pretrained_name).to("cuda").half()
        return client

    @app.function(image=image, gpu="A10G", volumes=volumes, secrets=secrets, timeout=900)
    def embed_residues(
        sequences: list[str], model_id: str = "esm3-sm-open-v1"
    ) -> list[list[list[float]]]:
        from esm.sdk.api import ESMProtein, LogitsConfig

        client = _load_client(model_id)
        by_length: dict[int, list[tuple[int, str]]] = {}
        for i, seq in enumerate(sequences):
            by_length.setdefault(len(seq), []).append((i, seq))

        out: list[list[list[float]] | None] = [None] * len(sequences)
        for _length, items in by_length.items():
            for i, seq in items:
                protein = ESMProtein(sequence=seq)
                protein_tensor = client.encode(protein)
                logits_output = client.logits(
                    protein_tensor, LogitsConfig(return_embeddings=True)
                )
                embeddings = logits_output.embeddings.squeeze(0).float().cpu().numpy()
                out[i] = embeddings.tolist()
        return out

    @app.function(image=image, gpu="A10G", volumes=volumes, secrets=secrets, timeout=900)
    def generate_masked(
        masked_sequence: str, num_candidates: int, model_id: str = "esm3-sm-open-v1"
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
    def log_likelihood(sequence: str, model_id: str = "esm3-sm-open-v1") -> float:
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
