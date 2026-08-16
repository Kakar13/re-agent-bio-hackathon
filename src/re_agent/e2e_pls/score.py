"""Peptide-level scoring and protein-level risk aggregation.

Tiles a parent sequence into overlapping 9-mers (the biologically correct
unit per the plan's scope note), scores each with the three calibrated
heads, and rolls windows up into a protein-level risk summary.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from re_agent.e2e_pls.esm3_modal import ESM3Client
from re_agent.e2e_pls.model import ThreeHeadModel

DEFAULT_PEPTIDE_LEN = 9
DEFAULT_FLANK_LEN = 4


@dataclass
class Window:
    peptide: str
    start: int
    end: int
    n_flank: str
    c_flank: str


def tile_sequence(
    sequence: str, peptide_len: int = DEFAULT_PEPTIDE_LEN, flank_len: int = DEFAULT_FLANK_LEN
) -> list[Window]:
    n = len(sequence)
    windows = []
    for start in range(0, n - peptide_len + 1):
        end = start + peptide_len
        windows.append(
            Window(
                peptide=sequence[start:end],
                start=start,
                end=end,
                n_flank=sequence[max(0, start - flank_len) : start],
                c_flank=sequence[end : end + flank_len],
            )
        )
    return windows


@dataclass
class PeptideScore:
    peptide: str
    start: int
    end: int
    cleave_n_prob: float
    cleave_c_prob: float
    tap_log_ic50_relative: float
    tap_uncertainty: float
    mhc_cosine_similarity: float
    mhc_presentation_propensity: float
    composite_risk: float

    def to_dict(self) -> dict:
        return asdict(self)


def compute_peptide_risk(
    cleave_n_prob: float, cleave_c_prob: float, mhc_presentation_propensity: float
) -> float:
    """Geometric mean of the three calibrated-probability components.

    TAP is reported alongside but excluded here: the TAP head is a
    regularized ridge regressor with bootstrap uncertainty, not a
    calibrated probability (see plan: "no deep TAP head claim is
    supportable at this size"). This treats N-cleavage, C-cleavage, and
    presentation as independent gates -- an approximation, not a
    measured joint probability.
    """
    components = np.clip([cleave_n_prob, cleave_c_prob, mhc_presentation_propensity], 1e-6, 1.0)
    return float(np.exp(np.mean(np.log(components))))


def score_window(
    window: Window, hla_allele: str, esm_client: ESM3Client, heads: ThreeHeadModel
) -> PeptideScore:
    n_vec = esm_client.embed_pooled(window.peptide, window.n_flank, window.c_flank, "cleave_n")
    c_vec = esm_client.embed_pooled(window.peptide, window.n_flank, window.c_flank, "cleave_c")
    mer_vec = esm_client.embed_pooled(window.peptide, window.n_flank, window.c_flank, "mean_9mer")

    cleave_n, cleave_c = heads.cleavage.predict(n_vec[None, :], c_vec[None, :])
    tap_mean, tap_std = heads.tap.predict_with_uncertainty(mer_vec[None, :])
    mhc = heads.mhc.score(mer_vec, hla_allele)

    risk = compute_peptide_risk(
        float(cleave_n[0]), float(cleave_c[0]), mhc["presentation_propensity"]
    )
    return PeptideScore(
        peptide=window.peptide,
        start=window.start,
        end=window.end,
        cleave_n_prob=float(cleave_n[0]),
        cleave_c_prob=float(cleave_c[0]),
        tap_log_ic50_relative=float(tap_mean[0]),
        tap_uncertainty=float(tap_std[0]),
        mhc_cosine_similarity=mhc["cosine_similarity"],
        mhc_presentation_propensity=mhc["presentation_propensity"],
        composite_risk=risk,
    )


def score_sequence(
    sequence: str,
    hla_allele: str,
    esm_client: ESM3Client,
    heads: ThreeHeadModel,
    peptide_len: int = DEFAULT_PEPTIDE_LEN,
    flank_len: int = DEFAULT_FLANK_LEN,
) -> list[PeptideScore]:
    windows = tile_sequence(sequence, peptide_len, flank_len)
    return [score_window(w, hla_allele, esm_client, heads) for w in windows]


@dataclass
class ProteinRisk:
    n_windows: int
    top_k: int
    threshold: float
    top_k_mean_risk: float
    max_risk: float
    max_risk_window: PeptideScore | None
    count_above_threshold: int

    def to_dict(self) -> dict:
        d = asdict(self)
        d["max_risk_window"] = self.max_risk_window.to_dict() if self.max_risk_window else None
        return d


def aggregate_protein_risk(
    scores: list[PeptideScore], top_k: int = 5, threshold: float = 0.5
) -> ProteinRisk:
    if not scores:
        return ProteinRisk(0, top_k, threshold, 0.0, 0.0, None, 0)
    risks = np.array([s.composite_risk for s in scores])
    order = np.argsort(-risks)
    top_k_risks = risks[order[:top_k]]
    return ProteinRisk(
        n_windows=len(scores),
        top_k=top_k,
        threshold=threshold,
        top_k_mean_risk=float(top_k_risks.mean()),
        max_risk=float(risks[order[0]]),
        max_risk_window=scores[order[0]],
        count_above_threshold=int((risks >= threshold).sum()),
    )
