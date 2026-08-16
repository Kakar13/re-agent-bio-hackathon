"""Per-residue solvent accessibility for joining structure to epitope risk.

Relative solvent accessibility (RSA) answers a design question rather than an
immunological one. MHC-I peptides come from proteasomal degradation of unfolded
protein, so burial does not decide whether an epitope is presented. What burial
decides is whether the designer can do anything about it: a flagged 9-mer on a
surface loop can be resampled freely, while one packed into the core cannot be
changed without risking the fold.
"""

from __future__ import annotations

import numpy as np

# Theoretical maximum accessible surface area per residue, Tien et al. (2013),
# PLoS ONE 8(11): e80635, table 1 "theoretical" column, in square angstroms.
MAX_ACCESSIBLE_SURFACE_AREA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "E": 223.0,
    "Q": 225.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}
# Rost and Sander's conventional cutoff for calling a residue buried.
BURIED_RSA_THRESHOLD = 0.25
MIN_OBSERVED_RESIDUES_PER_WINDOW = 5


def relative_accessibility(residue_codes: list[str], absolute_sasa: np.ndarray) -> np.ndarray:
    """Normalize absolute SASA by each residue's theoretical maximum."""

    absolute_sasa = np.asarray(absolute_sasa, dtype=np.float64)
    if len(residue_codes) != len(absolute_sasa):
        raise ValueError("residue codes and SASA values must align")
    maxima = np.array(
        [MAX_ACCESSIBLE_SURFACE_AREA.get(code, np.nan) for code in residue_codes],
        dtype=np.float64,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        rsa = absolute_sasa / maxima
    # Values slightly above 1 occur for terminal and highly extended residues.
    return np.clip(rsa, 0.0, 1.5)


def align_chain_to_design(design_sequence: str, chain_sequence: str) -> int | None:
    """Return the design offset at which the observed chain starts, or None.

    Deposited structures routinely omit disordered termini, so the observed
    sequence is usually a substring of the design rather than equal to it.
    """

    if not chain_sequence or not design_sequence:
        return None
    if chain_sequence == design_sequence:
        return 0
    index = design_sequence.find(chain_sequence)
    if index >= 0:
        return index
    # The reverse case: the deposited chain carries an expression tag the
    # design record does not. Anchor on the design instead.
    index = chain_sequence.find(design_sequence)
    if index >= 0:
        return -index
    return None


def project_onto_design(
    design_length: int,
    offset: int,
    chain_rsa: np.ndarray,
) -> np.ndarray:
    """Place observed per-residue RSA onto design coordinates, NaN elsewhere."""

    projected = np.full(design_length, np.nan, dtype=np.float64)
    for chain_index, value in enumerate(np.asarray(chain_rsa, dtype=np.float64)):
        design_index = offset + chain_index
        if 0 <= design_index < design_length:
            projected[design_index] = value
    return projected


def window_accessibility(
    projected_rsa: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    *,
    min_observed: int = MIN_OBSERVED_RESIDUES_PER_WINDOW,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Summarize RSA over each 9-mer window.

    Returns mean RSA, the fraction of observed residues that are buried, and the
    count of residues that had structural coverage. Windows with too little
    coverage return NaN rather than a mean over one or two residues.
    """

    projected_rsa = np.asarray(projected_rsa, dtype=np.float64)
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if starts.shape != ends.shape:
        raise ValueError("window starts and ends must align")

    mean_rsa = np.full(len(starts), np.nan, dtype=np.float64)
    buried_fraction = np.full(len(starts), np.nan, dtype=np.float64)
    observed = np.zeros(len(starts), dtype=np.int64)
    for index, (start, end) in enumerate(zip(starts, ends, strict=True)):
        values = projected_rsa[start:end]
        present = values[~np.isnan(values)]
        observed[index] = len(present)
        if len(present) >= min_observed:
            mean_rsa[index] = float(present.mean())
            buried_fraction[index] = float((present < BURIED_RSA_THRESHOLD).mean())
    return mean_rsa, buried_fraction, observed
