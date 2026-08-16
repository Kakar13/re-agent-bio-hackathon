"""Track 1 teacher labeling: Pepsickle (proteasomal cleavage), MHCflurry
(MHC-I binding affinity), and TAP imputation applied to the candidate pool.

Pepsickle pins `torch==1.13.1` / `scikit-learn==1.2.0`, which conflicts with
this project's main environment (`torch>=2.2`). It runs in an isolated venv
at `.tools/pepsickle` (create with `scripts/setup_pepsickle_env.sh`) and is
invoked as a subprocess CLI -- see that script for setup.

Pepsickle predicts a cleavage probability at every residue position of a
*whole protein* in one call, not per-9-mer, so it's run once per unique
parent sequence and then sliced per window -- far cheaper than one call per
peptide. MHCflurry runs in-process (its newer releases use PyTorch, no
env conflict).

TAP imputation: measured TAP data is scarce (613 DS613 peptides). To give
every row a TAP value for scoring/steering coverage, `impute_tap_labels()`
trains a ridge regression on one-hot 9-mer features of those 613 measured
points and applies it to the rest. Every imputed row is flagged in
`label_model_version` with `tap_ridge=ds613_v1`; downstream trainers of the
TAP head must still fit only on the measured rows (study_id == "ds613")
rather than on our own imputations, or they will learn to reproduce the
imputer instead of the underlying signal.
"""

from __future__ import annotations

import subprocess
import tempfile
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PEPSICKLE_PYTHON = REPO_ROOT / ".tools" / "pepsickle" / "bin" / "python"


def pepsickle_version(pepsickle_python: Path = DEFAULT_PEPSICKLE_PYTHON) -> str:
    result = subprocess.run(
        [
            str(pepsickle_python),
            "-c",
            "from importlib.metadata import version; print(version('pepsickle'))",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=True,
    )
    return result.stdout.strip()


def mhcflurry_version() -> str:
    import mhcflurry

    return mhcflurry.__version__


def run_pepsickle_on_proteins(
    sequences: dict[str, str], pepsickle_python: Path = DEFAULT_PEPSICKLE_PYTHON
) -> pd.DataFrame:
    """Returns a long table: columns [protein_id, position (1-indexed), cleav_prob]."""
    if not pepsickle_python.exists():
        raise FileNotFoundError(
            f"pepsickle isolated env not found at {pepsickle_python}. Set it up with:\n"
            f"  uv venv .tools/pepsickle --python 3.11 && "
            f"uv pip install --python .tools/pepsickle/bin/python pepsickle"
        )
    pepsickle_bin = pepsickle_python.parent / "pepsickle"

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = Path(tmp) / "proteins.fasta"
        out_path = Path(tmp) / "cleavage.tsv"
        with fasta_path.open("w") as f:
            for protein_id, seq in sequences.items():
                f.write(f">{protein_id}\n{seq}\n")
        result = subprocess.run(
            [str(pepsickle_bin), "-f", str(fasta_path), "-o", str(out_path)],
            capture_output=True,
            text=True,
            timeout=1800,
        )
        if result.returncode != 0:
            raise RuntimeError(f"pepsickle failed (exit {result.returncode}): {result.stderr}")
        table = pd.read_csv(out_path, sep="\t")
    return table.rename(columns={"protein_id": "parent_sequence_id"})


def attach_cleavage_labels(df: pd.DataFrame, cleavage_table: pd.DataFrame) -> pd.DataFrame:
    """Look up each window's N/C-terminal cleavage probability from the
    per-protein table. `start`/`end` are 0-indexed half-open bounds, which
    line up exactly with pepsickle's 1-indexed positions: the residue right
    before the window (last n_flank residue) is 1-indexed position `start`,
    and the window's own last residue is 1-indexed position `end`.
    """
    df = df.drop(columns=["cleave_n_prob", "cleave_c_prob"], errors="ignore").copy()
    n_lookup = cleavage_table.rename(columns={"position": "start", "cleav_prob": "cleave_n_prob"})[
        ["parent_sequence_id", "start", "cleave_n_prob"]
    ]
    c_lookup = cleavage_table.rename(columns={"position": "end", "cleav_prob": "cleave_c_prob"})[
        ["parent_sequence_id", "end", "cleave_c_prob"]
    ]
    df = df.merge(n_lookup, on=["parent_sequence_id", "start"], how="left")
    df = df.merge(c_lookup, on=["parent_sequence_id", "end"], how="left")
    return df


def attach_mhc_labels(df: pd.DataFrame, only_missing: bool = True) -> pd.DataFrame:
    """Batched MHCflurry binding-affinity prediction, grouped by allele.
    Binding-only output (not the presentation/processing predictor), per
    the plan's instruction to avoid double-counting antigen processing.
    """
    from mhcflurry import Class1AffinityPredictor

    df = df.copy()
    if "mhc_affinity_nm" not in df.columns:
        df["mhc_affinity_nm"] = np.nan
        df["mhc_percentile"] = np.nan

    target_mask = df["mhc_affinity_nm"].isna() if only_missing else pd.Series(True, index=df.index)
    if not target_mask.any():
        return df

    predictor = Class1AffinityPredictor.load()
    for allele, group in df[target_mask].groupby("hla_allele"):
        preds = predictor.predict_to_dataframe(peptides=group["peptide"].tolist(), allele=allele)
        df.loc[group.index, "mhc_affinity_nm"] = preds["prediction"].to_numpy()
        df.loc[group.index, "mhc_percentile"] = preds["prediction_percentile"].to_numpy()
    return df


def stamp_label_model_version(
    df: pd.DataFrame, pepsickle_python: Path = DEFAULT_PEPSICKLE_PYTHON
) -> pd.DataFrame:
    df = df.copy()
    try:
        pep_v = pepsickle_version(pepsickle_python)
    except (FileNotFoundError, subprocess.CalledProcessError):
        pep_v = "unavailable"
    version_string = f"pepsickle={pep_v};mhcflurry={mhcflurry_version()}"
    df["label_model_version"] = np.where(
        df["label_model_version"] == "", version_string, df["label_model_version"]
    )
    return df


TAP_RIDGE_TAG = "tap_ridge=ds613_v1"
_AMINO_ACID_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_AA_INDEX = {aa: i for i, aa in enumerate(_AMINO_ACID_ALPHABET)}


def _one_hot_9mer(peptide: str) -> np.ndarray:
    """9 positions x 20 amino acids -> length-180 vector. Classic featurization
    used across the TAP-prediction literature (SVMTAP, PREDTAP, TAPPred).

    Non-canonical residues are silently skipped (their position stays zero).
    """
    v = np.zeros(9 * 20, dtype=np.float32)
    for pos, aa in enumerate(peptide[:9]):
        idx = _AA_INDEX.get(aa)
        if idx is not None:
            v[pos * 20 + idx] = 1.0
    return v


def _featurize_peptides(peptides: pd.Series) -> np.ndarray:
    return np.stack([_one_hot_9mer(p) for p in peptides])


def impute_tap_labels(df: pd.DataFrame, ridge_alpha: float = 10.0) -> pd.DataFrame:
    """Fill missing `tap_log_ic50_relative` with predictions from a ridge model
    trained on the rows that already carry a measured TAP value.

    NOT a source of new data -- this only spreads the measured 613-peptide
    signal across the rest of the dataset via a linear model over one-hot
    9-mer features. Imputed rows are flagged in `label_model_version` with
    `TAP_RIDGE_TAG`. Downstream TapHead training in `train.py` must filter
    to measured rows only (`study_id == "ds613"`) to avoid learning from
    its own imputations.
    """
    from sklearn.linear_model import Ridge

    df = df.copy()
    measured_mask = df["tap_log_ic50_relative"].notna().values
    missing_mask = ~measured_mask

    if measured_mask.sum() == 0 or missing_mask.sum() == 0:
        return df

    x_train = _featurize_peptides(df.loc[measured_mask, "peptide"])
    y_train = df.loc[measured_mask, "tap_log_ic50_relative"].to_numpy(dtype=np.float64)
    ridge = Ridge(alpha=ridge_alpha).fit(x_train, y_train)

    x_missing = _featurize_peptides(df.loc[missing_mask, "peptide"])
    df.loc[missing_mask, "tap_log_ic50_relative"] = ridge.predict(x_missing).astype(np.float32)
    df.loc[missing_mask, "label_model_version"] = (
        df.loc[missing_mask, "label_model_version"].fillna("").astype(str) + ";" + TAP_RIDGE_TAG
    ).str.lstrip(";")
    return df


def has_measured_tap(df: pd.DataFrame) -> pd.Series:
    """Boolean mask: rows whose TAP value was NOT written by `impute_tap_labels`.

    In the real Track 1 build this is exactly the DS613 rows. In tests/fixtures
    that carry synthetic TAP values without going through the imputer, it also
    returns True -- the semantic is "TAP is not our own extrapolation of DS613,"
    which is what the downstream TapHead needs to avoid learning to recover
    the imputer instead of the underlying signal.
    """
    origin = df["label_model_version"].fillna("").astype(str)
    return ~origin.str.contains(TAP_RIDGE_TAG, regex=False)


def read_pepsickle_output(path: str | Path) -> pd.DataFrame:
    """Helper for tests/inspection: parse a saved pepsickle TSV without rerunning it."""
    return pd.read_csv(StringIO(Path(path).read_text()), sep="\t").rename(
        columns={"protein_id": "parent_sequence_id"}
    )
