"""The Track 1 / Track 2 integration boundary.

This module is the single source of truth for the aligned 9-mer row
contract: field names, units, allowed missingness, split names, and the
dataset-version hash. Track 1 (dataset) and Track 2 (model) both import
this module instead of re-declaring field names, so a change here is a
breaking-change to both sides at once.

Row unit: one peptide/HLA pair, normally a 9-mer, with its parent-protein
context (flanks) and provenance. See `e2e-pls-track-c` plan, "Exact row
contract".
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

import pandas as pd
import pyarrow as pa

SCHEMA_VERSION = "1"

ENCODER_MODEL_ID = "esm2_t33_650M_UR50D"
EMBEDDING_DIM = 1280

CANONICAL_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")

SOURCE_DOMAINS = ("natural_human", "natural_viral", "natural_bacterial", "de_novo", "demo")
LABEL_ORIGINS = ("measured", "teacher")
SPLIT_NAMES = ("train", "val", "test")

# name -> (pyarrow type, nullable)
REQUIRED_FIELDS: dict[str, tuple[pa.DataType, bool]] = {
    # Identity / context
    "row_id": (pa.string(), False),
    "peptide": (pa.string(), False),
    "length": (pa.int32(), False),
    "parent_sequence_id": (pa.string(), False),
    "parent_sequence_hash": (pa.string(), False),
    "start": (pa.int32(), False),
    "end": (pa.int32(), False),
    "n_flank": (pa.string(), True),
    "c_flank": (pa.string(), True),
    "source_domain": (pa.string(), False),
    # Targets (nullable: not every row carries every label at first pass)
    "cleave_n_prob": (pa.float32(), True),
    "cleave_c_prob": (pa.float32(), True),
    "tap_log_ic50_relative": (pa.float32(), True),
    "mhc_affinity_nm": (pa.float32(), True),
    "mhc_percentile": (pa.float32(), True),
    "hla_allele": (pa.string(), False),
    # Provenance
    "label_origin": (pa.string(), False),
    "label_model_version": (pa.string(), True),
    "source_uri": (pa.string(), True),
    "license": (pa.string(), False),
    "generated_at": (pa.string(), False),
    # Leakage controls
    "study_id": (pa.string(), False),
    "protein_cluster_id": (pa.string(), False),
    "peptide_cluster_id": (pa.string(), False),
    "split": (pa.string(), False),
    # Encoder artifacts. Embeddings themselves live in a separate memmap keyed
    # by embedding_cache_key -- see encoder.EmbeddingCache -- not duplicated here.
    "encoder_model_id": (pa.string(), False),
    "pooling_recipe": (pa.string(), False),
    "embedding_cache_key": (pa.string(), False),
    "has_structure_track": (pa.bool_(), False),
}

PYARROW_SCHEMA = pa.schema(
    [
        pa.field(name, dtype, nullable=nullable)
        for name, (dtype, nullable) in REQUIRED_FIELDS.items()
    ]
)

REQUIRED_COLUMNS = tuple(REQUIRED_FIELDS)


@dataclass
class ValidationResult:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def validate_dataframe(df: pd.DataFrame) -> ValidationResult:
    """Enforce the row contract and the build's quality gates.

    Returns a `ValidationResult`; `result.ok` is False if any hard gate
    (duplicate row IDs, coordinate/flank inconsistency, missing
    provenance, split-group overlap, noncanonical residues, invalid
    enums) fails. Missing target labels are a warning, not an error --
    not every source populates all three targets for every row.
    """
    result = ValidationResult()

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        result.errors.append(f"missing required columns: {missing_cols}")
        return result  # nothing else is checkable without the columns

    _check_required_non_null(df, result)
    _check_duplicate_row_ids(df, result)
    _check_coordinates(df, result)
    _check_residues(df, result)
    _check_enums(df, result)
    _check_split_leakage(df, result)
    _check_target_missingness(df, result)

    return result


def _check_required_non_null(df: pd.DataFrame, result: ValidationResult) -> None:
    for name, (_, nullable) in REQUIRED_FIELDS.items():
        if nullable:
            continue
        n_missing = df[name].isna().sum()
        if n_missing:
            result.errors.append(f"{name}: {n_missing} required (non-nullable) values are missing")


def _check_duplicate_row_ids(df: pd.DataFrame, result: ValidationResult) -> None:
    dupes = df["row_id"][df["row_id"].duplicated(keep=False)]
    if len(dupes):
        sample = sorted(dupes.unique())[:5]
        result.errors.append(f"duplicate row_id values ({dupes.nunique()} distinct, e.g. {sample})")


def _check_coordinates(df: pd.DataFrame, result: ValidationResult) -> None:
    bad_length = df[df["length"] != df["peptide"].str.len()]
    if len(bad_length):
        result.errors.append(f"{len(bad_length)} rows: length != len(peptide)")

    bad_span = df[(df["end"] - df["start"]) != df["length"]]
    if len(bad_span):
        result.errors.append(f"{len(bad_span)} rows: end - start != length")

    bad_order = df[df["start"] < 0]
    if len(bad_order):
        result.errors.append(f"{len(bad_order)} rows: start < 0")


def _check_residues(df: pd.DataFrame, result: ValidationResult) -> None:
    def _noncanonical(seq: str) -> bool:
        return bool(seq) and not set(seq.upper()) <= CANONICAL_RESIDUES

    bad_peptide = df[df["peptide"].map(_noncanonical)]
    if len(bad_peptide):
        result.errors.append(f"{len(bad_peptide)} rows: peptide has noncanonical residues")

    for flank_col in ("n_flank", "c_flank"):
        bad_flank = df[df[flank_col].fillna("").map(_noncanonical)]
        if len(bad_flank):
            result.errors.append(f"{len(bad_flank)} rows: {flank_col} has noncanonical residues")


def _check_enums(df: pd.DataFrame, result: ValidationResult) -> None:
    checks = {
        "source_domain": SOURCE_DOMAINS,
        "label_origin": LABEL_ORIGINS,
        "split": SPLIT_NAMES,
    }
    for col, allowed in checks.items():
        bad = df.loc[~df[col].isin(allowed), col]
        if len(bad):
            bad_values = sorted(bad.unique())[:5]
            result.errors.append(f"{col}: {bad.nunique()} rows outside {allowed}: {bad_values}")


def _check_split_leakage(df: pd.DataFrame, result: ValidationResult) -> None:
    for cluster_col in ("protein_cluster_id", "peptide_cluster_id"):
        splits_per_cluster = df.groupby(cluster_col)["split"].nunique()
        leaking = splits_per_cluster[splits_per_cluster > 1]
        if len(leaking):
            sample = list(leaking.index[:5])
            result.errors.append(
                f"split-group overlap on {cluster_col}: {len(leaking)} clusters span multiple "
                f"splits (e.g. {sample})"
            )


def _check_target_missingness(df: pd.DataFrame, result: ValidationResult) -> None:
    targets = (
        "cleave_n_prob",
        "cleave_c_prob",
        "tap_log_ic50_relative",
        "mhc_affinity_nm",
        "mhc_percentile",
    )
    for col in targets:
        n_missing = df[col].isna().sum()
        if n_missing:
            result.warnings.append(f"{col}: {n_missing}/{len(df)} rows missing this target")


def dataset_version_hash(df: pd.DataFrame) -> str:
    """Content hash for the dataset build: schema version + a stable row fingerprint.

    Deterministic given the same rows regardless of input ordering, so it
    can be recomputed and compared across a rebuild. Uses pandas' vectorized
    hasher so a multi-million-row frame does not walk Python row objects.
    """
    hasher = hashlib.sha256()
    hasher.update(f"schema={SCHEMA_VERSION}".encode())
    subset = df.loc[:, list(REQUIRED_COLUMNS)].sort_values("row_id", kind="mergesort")
    hashed = pd.util.hash_pandas_object(subset, index=False)
    hasher.update(hashed.to_numpy(dtype="uint64", copy=False).tobytes())
    return hasher.hexdigest()[:16]


def embedding_cache_key(
    peptide: str, n_flank: str, c_flank: str, model_id: str = ENCODER_MODEL_ID
) -> str:
    """Deterministic key for the memmap embedding cache in esm3_modal.EmbeddingCache."""
    raw = f"{model_id}|{n_flank or ''}|{peptide}|{c_flank or ''}"
    return hashlib.sha256(raw.encode()).hexdigest()[:24]
