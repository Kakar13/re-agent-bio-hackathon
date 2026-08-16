"""Resolve cleavage sites before any presentation prediction happens.

Three independent lines of evidence are combined per candidate bond:

* **Sequence preference** - the MEROPS log-odds matrix for each protease.
* **Composition-controlled significance** - the same matrix scored against
  shuffles of that same sequence. De novo designs have skewed composition, so a
  raw log-odds score is not comparable across designs; an empirical p-value
  against a matched null is. Shuffling preserves composition exactly.
* **Physical reachability** - relative solvent accessibility and B-factor over
  the P4-P4' span, from the design's own crystal structure.

The proteasome is handled separately by NetChop, since class I processing
happens in the cytosol rather than the endolysosome.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass

import numpy as np

from .accessibility import ResidueStructure, SpanFeatures, span_features
from .merops import POSITIONS, STANDARD_AA, SpecificityMatrix

log = logging.getLogger(__name__)

AA_INDEX = {aa: i for i, aa in enumerate(STANDARD_AA)}
N_FLANK = 4
WINDOW = len(POSITIONS)

DEFAULT_N_SHUFFLES = 200
# A site must beat this fraction of composition-matched random windows.
DEFAULT_P_THRESHOLD = 0.01
DEFAULT_NETCHOP_THRESHOLD = 0.5


@dataclass
class CleavageSite:
    """One protease-bond hypothesis with its evidence."""

    design_id: str
    cut_index: int  # bond lies between residues cut_index-1 and cut_index
    window: str  # P4-P4'
    protease: str
    pwm_score: float
    p_value: float
    z_score: float
    netchop_score: float | None
    mean_rsa: float | None
    mean_bfactor: float | None
    coil_fraction: float
    observed_fraction: float
    accessibility: str

    @property
    def scissile_bond(self) -> str:
        """Human-readable P1|P1' description, 1-based."""
        p1 = self.window[N_FLANK - 1]
        p1p = self.window[N_FLANK]
        return f"{p1}{self.cut_index}|{p1p}{self.cut_index + 1}"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["scissile_bond"] = self.scissile_bond
        return d


def encode(sequence: str) -> np.ndarray:
    """Sequence to matrix-column indices; unknown residues become -1."""
    return np.array([AA_INDEX.get(aa, -1) for aa in sequence], dtype=np.int16)


def pwm_array(matrix: SpecificityMatrix) -> np.ndarray:
    """The 8 x 20 log-odds array, ordered P4..P4'."""
    arr = np.zeros((WINDOW, len(STANDARD_AA)), dtype=np.float64)
    for i, pos in enumerate(POSITIONS):
        row = matrix.pwm.get(pos, {})
        for aa, j in AA_INDEX.items():
            arr[i, j] = row.get(aa, 0.0)
    return arr


def score_all_windows(codes: np.ndarray, pwm: np.ndarray) -> np.ndarray:
    """Score every full P4-P4' window in one pass.

    Result index ``k`` is the window starting at residue ``k``, whose scissile
    bond sits after ``k + 3``.
    """
    n = len(codes) - WINDOW + 1
    if n <= 0:
        return np.zeros(0)
    out = np.zeros(n, dtype=np.float64)
    for pos in range(WINDOW):
        col = codes[pos : pos + n]
        valid = col >= 0
        contrib = np.zeros(n)
        contrib[valid] = pwm[pos, col[valid]]
        out += contrib
    return out


def null_distribution(
    codes: np.ndarray, pwm: np.ndarray, *, n_shuffles: int, rng: np.random.Generator
) -> np.ndarray:
    """Window scores from composition-matched shuffles of the same sequence."""
    pool = []
    scrambled = codes.copy()
    for _ in range(n_shuffles):
        rng.shuffle(scrambled)
        pool.append(score_all_windows(scrambled, pwm))
    return np.concatenate(pool) if pool else np.zeros(0)


def scan_design(
    design_id: str,
    sequence: str,
    matrices: dict[str, SpecificityMatrix],
    *,
    features: list[ResidueStructure] | None = None,
    netchop: list[float] | None = None,
    n_shuffles: int = DEFAULT_N_SHUFFLES,
    seed: int = 0,
) -> list[CleavageSite]:
    """Score every bond in one design under every protease."""
    codes = encode(sequence)
    if len(codes) < WINDOW:
        return []
    rng = np.random.default_rng(seed)

    sites: list[CleavageSite] = []
    span_cache: dict[int, SpanFeatures] = {}

    for name, matrix in matrices.items():
        pwm = pwm_array(matrix)
        observed = score_all_windows(codes, pwm)
        null = null_distribution(codes, pwm, n_shuffles=n_shuffles, rng=rng)
        mu, sigma = (float(null.mean()), float(null.std())) if null.size else (0.0, 1.0)
        sigma = sigma or 1.0
        null_sorted = np.sort(null)

        for k, raw in enumerate(observed):
            cut_index = k + N_FLANK  # bond precedes this residue index
            # Right-tail empirical p-value, with the +1 correction that keeps a
            # p-value from ever being exactly zero.
            n_ge = null.size - int(np.searchsorted(null_sorted, raw, side="left"))
            p_value = (n_ge + 1) / (null.size + 1) if null.size else 1.0

            if cut_index not in span_cache:
                span_cache[cut_index] = (
                    span_features(features, cut_index)
                    if features
                    else SpanFeatures(None, None, 0.0, 0.0, "no-structure")
                )
            sf = span_cache[cut_index]

            nc = None
            if netchop and 0 <= cut_index - 1 < len(netchop):
                # NetChop scores the residue whose C-terminal bond is cut, which
                # is P1 - the residue immediately before the scissile bond.
                nc = float(netchop[cut_index - 1])

            sites.append(
                CleavageSite(
                    design_id=design_id,
                    cut_index=cut_index,
                    window=sequence[k : k + WINDOW],
                    protease=name,
                    pwm_score=float(raw),
                    p_value=float(p_value),
                    z_score=float((raw - mu) / sigma),
                    netchop_score=nc,
                    mean_rsa=sf.mean_rsa,
                    mean_bfactor=sf.mean_bfactor,
                    coil_fraction=sf.coil_fraction,
                    observed_fraction=sf.observed_fraction,
                    accessibility=sf.classification,
                )
            )
    return sites


def best_site_per_bond(sites: list[CleavageSite]) -> dict[int, CleavageSite]:
    """The most likely protease for each bond, by p-value then score."""
    best: dict[int, CleavageSite] = {}
    for s in sites:
        cur = best.get(s.cut_index)
        if cur is None or (s.p_value, -s.pwm_score) < (cur.p_value, -cur.pwm_score):
            best[s.cut_index] = s
    return best


def attribution_margin(sites: list[CleavageSite], cut_index: int) -> float:
    """Score gap between the best and second-best protease at one bond.

    A large margin means the site is characteristic of one enzyme; a small one
    means the attribution is not really distinguishable between proteases.
    """
    at_bond = sorted(
        (s for s in sites if s.cut_index == cut_index), key=lambda s: -s.pwm_score
    )
    if len(at_bond) < 2:
        return 0.0
    return at_bond[0].pwm_score - at_bond[1].pwm_score


def high_confidence_cuts(
    sites: list[CleavageSite],
    *,
    p_threshold: float = DEFAULT_P_THRESHOLD,
    require_accessible: bool = False,
) -> set[int]:
    """Bonds a cathepsin is confidently predicted to cut.

    These define the fragment boundaries the class II digest may not cross.
    """
    out = set()
    for s in sites:
        if s.p_value > p_threshold:
            continue
        if require_accessible and s.accessibility == "unfolding-required":
            continue
        out.add(s.cut_index)
    return out


def netchop_cut_positions(
    netchop: list[float], *, threshold: float = DEFAULT_NETCHOP_THRESHOLD
) -> set[int]:
    """Residue indices (0-based) whose C-terminal bond the proteasome cuts."""
    return {i for i, v in enumerate(netchop) if v >= threshold}
