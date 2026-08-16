"""Shared paths, constants, and hyperparameters for the immunogenicity model."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

# MHC class II presents a ~15-residue register; every prediction unit is a 15-mer.
WINDOW = 15

# ESM-2 150M: 30 layers, 640-dim residue embeddings. Small enough to embed the
# full corpus on a laptop, big enough to carry real structural signal.
ESM_MODEL = "esm2_t30_150M_UR50D"
ESM_LAYER = 30
ESM_DIM = 640

# The 20 canonical amino acids; anything else in a source sequence is dropped.
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = frozenset(AA)


@dataclass(frozen=True)
class TrainConfig:
    """Hyperparameters shared by the labeled-only baseline and Mean Teacher."""

    hidden: int = 256
    heads: int = 4
    dropout: float = 0.3
    lr: float = 3e-4
    weight_decay: float = 1e-4
    epochs: int = 30
    batch_labeled: int = 128
    batch_unlabeled: int = 128
    seed: int = 0

    # Mean Teacher. Consistency is ramped in so the student learns real labels
    # before it is asked to agree with a teacher that is still noise.
    consistency_weight: float = 10.0
    rampup_epochs: int = 10
    ema_decay_rampup: float = 0.99
    ema_decay_final: float = 0.999

    # Stochastic perturbations applied to the frozen embeddings.
    input_noise: float = 0.15
    feature_dropout: float = 0.1
    residue_mask_prob: float = 0.1


@dataclass(frozen=True)
class Paths:
    """Every artifact the pipeline reads or writes, in one place."""

    labeled: Path = field(default=DATA_PROCESSED / "iedb_class2_windows.parquet")
    unlabeled: Path = field(default=DATA_PROCESSED / "denovo_windows.parquet")
    reference: Path = field(default=DATA_PROCESSED / "reference_windows.parquet")
    embed_cache: Path = field(default=DATA_PROCESSED / "esm_cache")
    models: Path = field(default=RESULTS / "models")
    figures: Path = field(default=RESULTS / "figures")
    reports: Path = field(default=RESULTS / "reports")


PATHS = Paths()


def ensure_dirs() -> None:
    for path in (
        DATA_RAW,
        DATA_PROCESSED,
        PATHS.embed_cache,
        PATHS.models,
        PATHS.figures,
        PATHS.reports,
    ):
        path.mkdir(parents=True, exist_ok=True)


def torch_device() -> str:
    import torch

    if torch.backends.mps.is_available():
        return "mps"
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"
