"""Peptide-level scoring and protein-level risk aggregation.

Tiles a parent sequence into overlapping 9-mers (the biologically correct
unit per the plan's scope note), scores each with the three calibrated
heads, and rolls windows up into a protein-level risk summary.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from re_agent.e2e_pls.encoder import ProteinEncoder
from re_agent.e2e_pls.model import ThreeHeadModel

DEFAULT_PEPTIDE_LEN = 9
DEFAULT_FLANK_LEN = 4
# TAP bootstrap std of 1.0 log-IC50 unit maps to 0.5 tap-certainty.
TAP_UNCERTAINTY_SCALE = 1.0
# Max std of three values in [0, 1] is ~0.471 (e.g. 0, 0, 1).
_HEAD_AGREEMENT_STD_MAX = float(np.std([0.0, 0.0, 1.0]))


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
    pathway_rank_score: float
    confidence: float
    confidence_tap: float
    confidence_mhc: float
    confidence_context: float
    confidence_agreement: float

    def to_dict(self) -> dict:
        return asdict(self)


def compute_peptide_risk(
    cleave_n_prob: float, cleave_c_prob: float, mhc_presentation_propensity: float
) -> float:
    """Equal-weight geometric mean of N-cleavage, C-cleavage, and MHC.

    Kept as an inspectable processing heuristic. Campaign ranking should use
    :func:`compute_pathway_rank_score` instead of averaging EL, BA, and
    affinity, which are not three independent biological events.
    """
    components = np.clip([cleave_n_prob, cleave_c_prob, mhc_presentation_propensity], 1e-6, 1.0)
    return float(np.exp(np.mean(np.log(components))))


# Presentation dominates because MHC-I binding/display is the selective step
# (NetMHCpan EL is the ligand-presentation endpoint). Proteasomal generation
# is necessary but more promiscuous, so it is a soft gate, not an equal vote.
# BA and affinity are omitted: affinity is a monotone transform of BA, and BA
# plus EL would count the same MHC event twice.
PATHWAY_EL_WEIGHT = 0.70
PATHWAY_GENERATION_WEIGHT = 0.30


def compute_generation_factor(cleave_n_prob: float, cleave_c_prob: float) -> float:
    """Both termini must be cut. Geometric mean of the two cleavage heads."""

    n_prob, c_prob = np.clip([cleave_n_prob, cleave_c_prob], 1e-6, 1.0)
    return float(np.sqrt(n_prob * c_prob))


def compute_pathway_rank_score(
    el_presentation_propensity: float,
    cleave_n_prob: float,
    cleave_c_prob: float,
    *,
    el_weight: float = PATHWAY_EL_WEIGHT,
    generation_weight: float = PATHWAY_GENERATION_WEIGHT,
) -> float:
    """Campaign rank for one 9-mer.

    .. math::

        R = p_{EL}^{0.70}\\,\\bigl(\\sqrt{p_N\\,p_C}\\bigr)^{0.30}

    EL is the presentation endpoint. Generation is a soft gate: a window that
    looks uncuttable is discounted even if EL is high. This is a monotonic
    triage heuristic, not a measured pathway probability.
    """

    if abs(el_weight + generation_weight - 1.0) > 1e-9:
        raise ValueError("pathway weights must sum to 1")
    el = float(np.clip(el_presentation_propensity, 1e-6, 1.0))
    generation = compute_generation_factor(cleave_n_prob, cleave_c_prob)
    return float((el**el_weight) * (generation**generation_weight))


def compute_confidence(
    tap_uncertainty: float,
    mhc_presentation_propensity: float,
    cleave_n_prob: float,
    cleave_c_prob: float,
    n_flank: str = "",
    c_flank: str = "",
    flank_len: int = DEFAULT_FLANK_LEN,
    tap_uncertainty_scale: float = TAP_UNCERTAINTY_SCALE,
) -> dict[str, float]:
    """Decomposed [0, 1] confidence for one window. Not a calibrated coverage
    probability -- four measurable sources, geometrically averaged:

    - tap: bootstrap std of the TAP ridge (lower std -> higher certainty)
    - mhc: how far presentation propensity sits from the 0.5 decision boundary
    - context: fraction of expected N/C flank residues present
    - agreement: how tightly the three probability heads cluster
    """
    tap = 1.0 / (1.0 + max(float(tap_uncertainty), 0.0) / tap_uncertainty_scale)
    mhc = float(np.clip(2.0 * abs(mhc_presentation_propensity - 0.5), 0.0, 1.0))
    context = 0.5 * (
        min(len(n_flank), flank_len) / flank_len + min(len(c_flank), flank_len) / flank_len
    )
    probs = np.clip([cleave_n_prob, cleave_c_prob, mhc_presentation_propensity], 0.0, 1.0)
    agreement = float(
        np.clip(1.0 - np.std(probs) / _HEAD_AGREEMENT_STD_MAX, 0.0, 1.0)
    )
    parts = np.clip([tap, mhc, context, agreement], 1e-6, 1.0)
    score = float(np.exp(np.mean(np.log(parts))))
    return {
        "confidence": score,
        "confidence_tap": float(tap),
        "confidence_mhc": mhc,
        "confidence_context": float(context),
        "confidence_agreement": agreement,
    }


def score_window(
    window: Window, hla_allele: str, esm_client: ProteinEncoder, heads: ThreeHeadModel
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
    pathway_rank = compute_pathway_rank_score(
        mhc["presentation_propensity"],
        float(cleave_n[0]),
        float(cleave_c[0]),
    )
    conf = compute_confidence(
        tap_uncertainty=float(tap_std[0]),
        mhc_presentation_propensity=mhc["presentation_propensity"],
        cleave_n_prob=float(cleave_n[0]),
        cleave_c_prob=float(cleave_c[0]),
        n_flank=window.n_flank,
        c_flank=window.c_flank,
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
        pathway_rank_score=pathway_rank,
        **conf,
    )


def score_sequence(
    sequence: str,
    hla_allele: str,
    esm_client: ProteinEncoder,
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
    mean_confidence: float
    max_risk_confidence: float

    def to_dict(self) -> dict:
        d = asdict(self)
        d["max_risk_window"] = self.max_risk_window.to_dict() if self.max_risk_window else None
        return d


def aggregate_protein_risk(
    scores: list[PeptideScore], top_k: int = 5, threshold: float = 0.5
) -> ProteinRisk:
    if not scores:
        return ProteinRisk(0, top_k, threshold, 0.0, 0.0, None, 0, 0.0, 0.0)
    risks = np.array([s.pathway_rank_score for s in scores])
    order = np.argsort(-risks)
    top_k_risks = risks[order[:top_k]]
    max_window = scores[order[0]]
    return ProteinRisk(
        n_windows=len(scores),
        top_k=top_k,
        threshold=threshold,
        top_k_mean_risk=float(top_k_risks.mean()),
        max_risk=float(risks[order[0]]),
        max_risk_window=max_window,
        count_above_threshold=int((risks >= threshold).sum()),
        mean_confidence=float(np.mean([s.confidence for s in scores])),
        max_risk_confidence=max_window.confidence,
    )


# --------------------------------------------------------------------------
# Residue heatmap: roll window scores onto each amino acid
# --------------------------------------------------------------------------


@dataclass
class ResidueAttribution:
    """Per-residue importance for the sequence visual.

    `risk_max` / `risk_mean` come from every 9-mer window that covers this
    position. That is the prediction the judges see, not transformer
    attention -- the heads have no attention maps. `mhc_importance` is the
    optional extra track: how aligned this residue's ESM vector is with the
    allele's binder centroid (same signal Steer-to-Safety uses to pick a
    mutation site).
    """

    position: int  # 0-indexed
    residue: str
    risk_max: float
    risk_mean: float
    confidence_at_max: float
    n_windows: int
    mhc_importance: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def compute_residue_heatmap(sequence: str, scores: list[PeptideScore]) -> list[ResidueAttribution]:
    """Paint each residue with the risk of the 9-mers that cover it."""
    n = len(sequence)
    risks: list[list[float]] = [[] for _ in range(n)]
    confs: list[list[float]] = [[] for _ in range(n)]
    for s in scores:
        for i in range(s.start, s.end):
            if 0 <= i < n:
                risks[i].append(s.composite_risk)
                confs[i].append(s.confidence)

    out: list[ResidueAttribution] = []
    for i, aa in enumerate(sequence):
        if not risks[i]:
            out.append(
                ResidueAttribution(
                    position=i,
                    residue=aa,
                    risk_max=0.0,
                    risk_mean=0.0,
                    confidence_at_max=0.0,
                    n_windows=0,
                )
            )
            continue
        r = np.asarray(risks[i], dtype=np.float64)
        c = np.asarray(confs[i], dtype=np.float64)
        imax = int(r.argmax())
        out.append(
            ResidueAttribution(
                position=i,
                residue=aa,
                risk_max=float(r.max()),
                risk_mean=float(r.mean()),
                confidence_at_max=float(c[imax]),
                n_windows=int(r.size),
            )
        )
    return out


def attach_mhc_residue_importance(
    sequence: str,
    attributions: list[ResidueAttribution],
    esm_client: ProteinEncoder,
    heads: ThreeHeadModel,
    hla_allele: str,
) -> list[ResidueAttribution]:
    """Fill `mhc_importance` from per-residue cosine to the binder centroid.

    Mapped from cosine [-1, 1] to [0, 1] so it shares a color scale with risk.
    """
    if hla_allele not in heads.mhc.centroids or not sequence:
        return attributions
    residue_vecs = esm_client.embed_residues(sequence)
    if residue_vecs.shape[0] != len(sequence):
        return attributions
    z = heads.mhc.project(residue_vecs)
    cosine = z @ heads.mhc.centroids[hla_allele]
    mapped = np.clip((cosine + 1.0) / 2.0, 0.0, 1.0)
    for i, attr in enumerate(attributions):
        attr.mhc_importance = float(mapped[i])
    return attributions


def explain_sequence(
    sequence: str,
    hla_allele: str,
    esm_client: ProteinEncoder,
    heads: ThreeHeadModel,
    peptide_len: int = DEFAULT_PEPTIDE_LEN,
    flank_len: int = DEFAULT_FLANK_LEN,
) -> tuple[list[PeptideScore], ProteinRisk, list[ResidueAttribution]]:
    """Score windows, roll up protein risk, and build the residue heatmap."""
    scores = score_sequence(sequence, hla_allele, esm_client, heads, peptide_len, flank_len)
    protein = aggregate_protein_risk(scores)
    heatmap = compute_residue_heatmap(sequence, scores)
    attach_mhc_residue_importance(sequence, heatmap, esm_client, heads, hla_allele)
    return scores, protein, heatmap
