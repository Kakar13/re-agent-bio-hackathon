"""Constrained latent-guided sequence steering ("Steer to Safety").

Latent-space distance-to-centroid guides *which* position to mutate next
and *ranks* candidates; final acceptance is always verified on the real,
re-embedded discrete sequence against the teacher-calibrated head scores
-- never on the latent target alone. At most `max_mutations` are ever
applied, and positions in `protected_positions` are never touched.

TAP regression direction assumes higher `tap_log_ic50_relative` is worse
transport (weaker TAP affinity), matching the convention used for
`mhc_affinity_nm`/IC50-style teacher scores; flag if DS613's actual sign
convention differs when Track 1's real labels land.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from re_agent.e2e_pls.esm3_modal import ESM3Client
from re_agent.e2e_pls.model import ThreeHeadModel
from re_agent.e2e_pls.score import DEFAULT_FLANK_LEN, PeptideScore, Window, score_window

CLAIM_DISCLAIMER = (
    "Experimental cleavage/TAP/MHC-I surrogate scores and predicted presentation propensity from "
    "calibrated teacher-trained heads. Not measured immunogenicity, T-cell response, safety, or "
    "retained function."
)


@dataclass
class SteeringConfig:
    hla_allele: str = "HLA-A*02:01"
    max_mutations: int = 3
    max_steps: int = 6
    mask_candidates_per_step: int = 8
    protected_positions: frozenset[int] = field(default_factory=frozenset)
    latent_movement_penalty: float = 0.25
    cleavage_regression_penalty: float = 1.0
    tap_regression_penalty: float = 0.5
    preservation_threshold: float = 0.5  # max allowed drop in per-residue sequence log-likelihood

    def to_dict(self) -> dict:
        d = asdict(self)
        d["protected_positions"] = sorted(self.protected_positions)
        return d


@dataclass
class CandidateRecord:
    position: int
    from_residue: str
    to_residue: str
    composite_risk: float
    cleave_n_prob: float
    cleave_c_prob: float
    tap_log_ic50_relative: float
    latent_distance: float
    sequence_log_likelihood: float
    objective: float
    accepted: bool
    rejection_reason: str | None


@dataclass
class StepRecord:
    step: int
    masked_position: int
    candidates: list[CandidateRecord]
    accepted_candidate: CandidateRecord | None
    duration_s: float


@dataclass
class SteeringTrace:
    input_sequence: str
    output_sequence: str
    target_window_start: int
    target_window_end: int
    hla_allele: str
    mutations: list[dict]
    initial_score: dict
    final_score: dict
    steps: list[StepRecord]
    model_version: str
    esm3_model_id: str
    dataset_version_hash: str
    constraints: dict
    total_duration_s: float
    claim_disclaimer: str = CLAIM_DISCLAIMER

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2))

    def optimization_path(self) -> list[dict]:
        """One row per step for the dashboard's dual-axis convergence chart."""
        rows = []
        for s in self.steps:
            c = s.accepted_candidate
            rows.append(
                {
                    "step": s.step,
                    "position": s.masked_position,
                    "accepted": c is not None,
                    "composite_risk": c.composite_risk if c else None,
                    "latent_distance": c.latent_distance if c else None,
                    "n_candidates": len(s.candidates),
                }
            )
        return rows


def _pick_next_position(
    sequence: str,
    target_start: int,
    target_end: int,
    candidate_positions: list[int],
    esm_client: ESM3Client,
    heads: ThreeHeadModel,
    hla_allele: str,
) -> int:
    """Greedy pick: the unprotected position whose residue embedding is most
    aligned with the binder centroid (i.e. most implicated in presentation risk).
    """
    window_start = max(0, target_start - DEFAULT_FLANK_LEN)
    window_end = min(len(sequence), target_end + DEFAULT_FLANK_LEN)
    residues = esm_client.embed_residues(sequence[window_start:window_end])
    centroid = heads.mhc.centroids[hla_allele]

    best_pos, best_score = candidate_positions[0], -np.inf
    for p in candidate_positions:
        local_idx = p - window_start
        z = heads.mhc.project(residues[local_idx][None, :])[0]
        s = float(z @ centroid)
        if s > best_score:
            best_score, best_pos = s, p
    return best_pos


def _window_at(sequence: str, target_start: int, target_end: int) -> Window:
    return Window(
        peptide=sequence[target_start:target_end],
        start=target_start,
        end=target_end,
        n_flank=sequence[max(0, target_start - DEFAULT_FLANK_LEN) : target_start],
        c_flank=sequence[target_end : target_end + DEFAULT_FLANK_LEN],
    )


def _evaluate_candidate(
    cand_seq: str,
    position: int,
    from_residue: str,
    target_start: int,
    target_end: int,
    hla_allele: str,
    esm_client: ESM3Client,
    heads: ThreeHeadModel,
    orig_score: PeptideScore,
    orig_z: np.ndarray,
    orig_log_likelihood: float,
    config: SteeringConfig,
) -> tuple[CandidateRecord, PeptideScore]:
    window = _window_at(cand_seq, target_start, target_end)
    window_score = score_window(window, hla_allele, esm_client, heads)

    mer_vec = esm_client.embed_pooled(window.peptide, window.n_flank, window.c_flank, "mean_9mer")
    z = heads.mhc.project(mer_vec[None, :])[0]
    latent_distance = float(np.linalg.norm(z - orig_z))

    cleavage_regression = max(0.0, orig_score.cleave_n_prob - window_score.cleave_n_prob) + max(
        0.0, orig_score.cleave_c_prob - window_score.cleave_c_prob
    )
    tap_regression = max(0.0, window_score.tap_log_ic50_relative - orig_score.tap_log_ic50_relative)

    objective = (
        window_score.composite_risk
        + config.latent_movement_penalty * latent_distance
        + config.cleavage_regression_penalty * cleavage_regression
        + config.tap_regression_penalty * tap_regression
    )

    log_likelihood = esm_client.sequence_log_likelihood(cand_seq)
    ll_delta_per_residue = (log_likelihood - orig_log_likelihood) / len(cand_seq)

    accepted, rejection_reason = True, None
    if window_score.composite_risk >= orig_score.composite_risk:
        accepted, rejection_reason = False, "did not improve composite risk vs. original"
    elif ll_delta_per_residue < -config.preservation_threshold:
        accepted, rejection_reason = (
            False,
            f"sequence log-likelihood dropped {ll_delta_per_residue:.3f}/residue, "
            "exceeds preservation threshold",
        )

    record = CandidateRecord(
        position=position,
        from_residue=from_residue,
        to_residue=cand_seq[position],
        composite_risk=window_score.composite_risk,
        cleave_n_prob=window_score.cleave_n_prob,
        cleave_c_prob=window_score.cleave_c_prob,
        tap_log_ic50_relative=window_score.tap_log_ic50_relative,
        latent_distance=latent_distance,
        sequence_log_likelihood=log_likelihood,
        objective=objective,
        accepted=accepted,
        rejection_reason=rejection_reason,
    )
    return record, window_score


def steer_to_safety(
    parent_sequence: str,
    target_start: int,
    target_end: int,
    esm_client: ESM3Client,
    heads: ThreeHeadModel,
    config: SteeringConfig | None = None,
) -> SteeringTrace:
    """Steer the peptide at `parent_sequence[target_start:target_end]` (usually
    the highest-risk window from `score.score_sequence`) toward lower presentation
    risk, subject to `config`'s mutation budget, protected positions, and
    preservation threshold.
    """
    config = config or SteeringConfig()
    t0 = time.monotonic()

    orig_window = _window_at(parent_sequence, target_start, target_end)
    orig_score = score_window(orig_window, config.hla_allele, esm_client, heads)
    orig_mer_vec = esm_client.embed_pooled(
        orig_window.peptide, orig_window.n_flank, orig_window.c_flank, "mean_9mer"
    )
    orig_z = heads.mhc.project(orig_mer_vec[None, :])[0]
    orig_log_likelihood = esm_client.sequence_log_likelihood(parent_sequence)

    current_sequence = parent_sequence
    current_score = orig_score
    mutations: list[dict] = []
    steps: list[StepRecord] = []
    remaining_positions = [
        p for p in range(target_start, target_end) if p not in config.protected_positions
    ]

    for step_i in range(config.max_steps):
        if len(mutations) >= config.max_mutations or not remaining_positions:
            break
        step_t0 = time.monotonic()

        position = _pick_next_position(
            current_sequence,
            target_start,
            target_end,
            remaining_positions,
            esm_client,
            heads,
            config.hla_allele,
        )
        remaining_positions.remove(position)
        from_residue = current_sequence[position]

        raw_candidates = esm_client.generate_masked(
            current_sequence, [position], num_candidates=config.mask_candidates_per_step
        )

        candidate_records: list[CandidateRecord] = []
        best: tuple[CandidateRecord, str, PeptideScore] | None = None
        for cand_seq in raw_candidates:
            if cand_seq[position] == from_residue:
                continue
            record, window_score = _evaluate_candidate(
                cand_seq,
                position,
                from_residue,
                target_start,
                target_end,
                config.hla_allele,
                esm_client,
                heads,
                orig_score,
                orig_z,
                orig_log_likelihood,
                config,
            )
            candidate_records.append(record)
            if record.accepted and (best is None or record.objective < best[0].objective):
                best = (record, cand_seq, window_score)

        steps.append(
            StepRecord(
                step=step_i,
                masked_position=position,
                candidates=candidate_records,
                accepted_candidate=best[0] if best else None,
                duration_s=time.monotonic() - step_t0,
            )
        )

        if best is not None:
            record, cand_seq, window_score = best
            mutations.append(
                {"position": position, "from": record.from_residue, "to": record.to_residue}
            )
            current_sequence = cand_seq
            current_score = window_score

    return SteeringTrace(
        input_sequence=parent_sequence,
        output_sequence=current_sequence,
        target_window_start=target_start,
        target_window_end=target_end,
        hla_allele=config.hla_allele,
        mutations=mutations,
        initial_score=orig_score.to_dict(),
        final_score=current_score.to_dict(),
        steps=steps,
        model_version=heads.model_version,
        esm3_model_id=heads.esm3_model_id,
        dataset_version_hash=heads.dataset_version_hash,
        constraints=config.to_dict(),
        total_duration_s=time.monotonic() - t0,
    )
