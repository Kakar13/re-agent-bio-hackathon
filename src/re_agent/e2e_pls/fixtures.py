"""Small, checked-in dev fixture conforming to `schema.py`.

Lets Track 2 (model/product) build and test against the final row
contract before Track 1's real parquet exists. All labels here are
synthetic (`label_origin="teacher"`, `label_model_version="fixture-synthetic-v1"`)
-- do not train a checkpoint intended for the dashboard on this fixture
alone, it exists for API/shape development and tests only.

Demo seed sequences (for the dashboard's "select a sequence" input) are
kept separate from the fixture table: they are scoring inputs, not
training rows, per the plan ("these are inputs, not training labels").
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.e2e_pls import schema

FIXTURE_DIR = Path(__file__).parent / "fixtures"
FIXTURE_PATH = FIXTURE_DIR / "dev_fixture.json"

_FIXTURE_HLA = "HLA-A*02:01"
_FIXTURE_TIMESTAMP = (
    "1970-01-01T00:00:00Z"  # placeholder: fixture is synthetic, not really generated then
)
_FLANK_LEN = 4
_PEPTIDE_LEN = 9

# Deterministic toy "parent proteins" -- not real sequences, just enough
# residue diversity to exercise tiling/flank/pooling code paths.
_NATURAL_PARENTS = {
    "natural_0": "MSTAVLENPGLGRKLSDFGQETSYIED",
    "natural_1": "MKTAYIAKQRQISFVKSHFSRQLEERL",
    "natural_2": "MADEEKLPPGWEKRMSRSSGRVYYFNH",
}
_DE_NOVO_PARENTS = {
    "de_novo_0": "GSHMSEEELKEAVKLLKKAEELVKKGD",
}

# protein_cluster_id -> split, assigned whole-protein so no cluster spans splits
_SPLIT_BY_PARENT = {
    "natural_0": "train",
    "natural_1": "train",
    "natural_2": "val",
    "de_novo_0": "test",
}

# Literature-cited HLA-A*02:01 epitopes + one placeholder de novo binder,
# for the dashboard's demo-seed picker. Not part of the training fixture.
DEMO_SEQUENCES: dict[str, dict[str, str]] = {
    "cmv_pp65_nlv": {
        "sequence": "NLVPMVATV",
        "kind": "viral_epitope",
        "description": "CMV pp65 495-503, canonical HLA-A*02:01 epitope (literature reference, not new data).",  # noqa: E501
    },
    "flu_m1_gil": {
        "sequence": "GILGFVFTL",
        "kind": "viral_epitope",
        "description": "Influenza A M1 58-66, canonical HLA-A*02:01 epitope (literature reference).",  # noqa: E501
    },
    "synthetic_neoantigen_demo": {
        "sequence": "MASNTVSAQGGSNRKVFVGGLSPDTSEEQIREYFGQFGEVESIELPMDPKLNKRRGFCFITF",
        "kind": "neoantigen",
        "description": "Synthetic placeholder sequence for the demo, not a real patient-derived neoantigen.",  # noqa: E501
    },
    "synthetic_designed_binder_demo": {
        "sequence": _DE_NOVO_PARENTS["de_novo_0"],
        "kind": "designed_binder",
        "description": "Synthetic placeholder de novo-style sequence for the demo, not an actual RFdiffusion output.",  # noqa: E501
    },
}


def _tile_parent(
    parent_id: str, sequence: str, source_domain: str, rng: np.random.Generator
) -> list[dict]:
    rows = []
    n = len(sequence)
    for start in range(0, n - _PEPTIDE_LEN + 1):
        end = start + _PEPTIDE_LEN
        peptide = sequence[start:end]
        n_flank = sequence[max(0, start - _FLANK_LEN) : start]
        c_flank = sequence[end : end + _FLANK_LEN]

        cleave_n = float(np.clip(rng.beta(2, 2), 0, 1))
        cleave_c = float(np.clip(rng.beta(2, 2), 0, 1))
        tap_score = float(rng.normal(0, 1))
        affinity_nm = float(np.exp(rng.uniform(np.log(1), np.log(50000))))
        percentile = float(np.clip(100 * affinity_nm / 50000, 0.01, 100))

        rows.append(
            {
                "row_id": f"{parent_id}:{start}-{end}",
                "peptide": peptide,
                "length": len(peptide),
                "parent_sequence_id": parent_id,
                "parent_sequence_hash": f"fixture-{parent_id}",
                "start": start,
                "end": end,
                "n_flank": n_flank,
                "c_flank": c_flank,
                "source_domain": source_domain,
                "cleave_n_prob": cleave_n,
                "cleave_c_prob": cleave_c,
                "tap_log_ic50_relative": tap_score,
                "mhc_affinity_nm": affinity_nm,
                "mhc_percentile": percentile,
                "hla_allele": _FIXTURE_HLA,
                "label_origin": "teacher",
                "label_model_version": "fixture-synthetic-v1",
                "source_uri": None,
                "license": "internal-fixture",
                "generated_at": _FIXTURE_TIMESTAMP,
                "study_id": f"fixture_study_{source_domain}",
                "protein_cluster_id": parent_id,
                "peptide_cluster_id": f"{parent_id}_nbhd{start // 3}",
                "split": _SPLIT_BY_PARENT[parent_id],
                "encoder_model_id": schema.ENCODER_MODEL_ID,
                "pooling_recipe": "mean_9mer",
                "embedding_cache_key": schema.embedding_cache_key(peptide, n_flank, c_flank),
                "has_structure_track": False,
            }
        )
    return rows


def generate_fixture(seed: int = 7) -> pd.DataFrame:
    """Regenerate the dev fixture deterministically. Does not write to disk."""
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for parent_id, sequence in _NATURAL_PARENTS.items():
        rows.extend(_tile_parent(parent_id, sequence, "natural_human", rng))
    for parent_id, sequence in _DE_NOVO_PARENTS.items():
        rows.extend(_tile_parent(parent_id, sequence, "de_novo", rng))

    df = pd.DataFrame(rows)
    for col, (dtype, _) in schema.REQUIRED_FIELDS.items():
        if pa_is_float(dtype):
            df[col] = df[col].astype("float32")
    return df[list(schema.REQUIRED_COLUMNS)]


def pa_is_float(dtype) -> bool:
    return str(dtype) in ("float", "float32", "float64", "double")


def write_dev_fixture(df: pd.DataFrame | None = None) -> Path:
    df = df if df is not None else generate_fixture()
    result = schema.validate_dataframe(df)
    if not result.ok:
        raise ValueError(f"generated fixture fails schema validation: {result.errors}")
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    FIXTURE_PATH.write_text(json.dumps(df.to_dict(orient="records"), indent=2, default=str))
    return FIXTURE_PATH


def load_dev_fixture() -> pd.DataFrame:
    """Load the checked-in fixture, generating+writing it on first use."""
    if not FIXTURE_PATH.exists():
        write_dev_fixture()
    records = json.loads(FIXTURE_PATH.read_text())
    df = pd.DataFrame.from_records(records)
    for col, (dtype, _) in schema.REQUIRED_FIELDS.items():
        if pa_is_float(dtype):
            df[col] = df[col].astype("float32")
    return df[list(schema.REQUIRED_COLUMNS)]
