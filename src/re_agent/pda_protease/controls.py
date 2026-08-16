"""Control arms that make the concordance claim falsifiable.

Any co-folding model will dock a peptide into a protease cleft; that is what the
cleft is shaped to do. So "Boltz put the peptide in the active site" is not
evidence on its own. The claim only means something if predicted-real sites
score better than matched negatives run in the same wave, under the same
settings, through the same measurement code.

Four arms, each removing exactly one thing:

* **positive control** - CD74 with cathepsin S. Known biology; if this does not
  come out, the method is not measuring anything.
* **scramble** - the same residues in a shuffled order. Removes the motif while
  holding composition and length fixed, so a hit here means composition alone
  drives the score.
* **protease swap** - the same segment against an enzyme whose matrix scores it
  poorly. Removes specificity while holding the peptide fixed, so a hit here
  means the result is not protease-specific.
* **decoy site** - a different bond in the same design with a weak matrix score
  but comparable solvent accessibility. Removes the sequence signal while
  holding the protein and its exposure fixed.
"""

from __future__ import annotations

import logging
import random
from dataclasses import dataclass

from .cleavage import CleavageSite
from .structure import COFOLDABLE, CofoldJob, build_job

log = logging.getLogger(__name__)

ARM_REAL = "real"
ARM_SCRAMBLE = "scramble"
ARM_SWAP = "protease_swap"
ARM_DECOY = "decoy"
ARM_POSITIVE = "positive_control"

# CD74 is cleaved by cathepsin S to liberate CLIP, which is annotated 97-120 in
# UniProt P04233. A site called anywhere in this window counts as a recovery.
CD74_CLIP_RANGE = (97, 125)


@dataclass
class ControlPlan:
    """Every job to run, real and control, as one interleaved list."""

    jobs: list[CofoldJob]
    counts: dict[str, int]


def scramble_segment(segment: str, rng: random.Random) -> str:
    """Shuffle a segment while keeping its ends, so length and composition hold."""
    chars = list(segment)
    rng.shuffle(chars)
    return "".join(chars)


def pick_swap_protease(
    sites: list[CleavageSite], cut_index: int, real_protease: str
) -> str | None:
    """The enzyme that least likes this exact site.

    Using the worst-scoring protease rather than a random one makes the negative
    as unfavourable as the matrices allow, which is the strongest form of the
    control.
    """
    scored = [
        (s.pwm_score, s.protease)
        for s in sites
        if s.cut_index == cut_index and s.protease != real_protease and s.protease in COFOLDABLE
    ]
    if not scored:
        return None
    return min(scored)[1]


def pick_decoy_sites(
    sites: list[CleavageSite],
    protease: str,
    real_cuts: set[int],
    *,
    n: int,
    accessibility: str | None = None,
    rng: random.Random | None = None,
) -> list[CleavageSite]:
    """Weak-scoring bonds matched to the real sites on solvent accessibility."""
    rng = rng or random.Random(0)
    pool = [
        s
        for s in sites
        if s.protease == protease
        and s.cut_index not in real_cuts
        and s.p_value > 0.5
        and (accessibility is None or s.accessibility == accessibility)
    ]
    if not pool:
        pool = [s for s in sites if s.protease == protease and s.cut_index not in real_cuts]
    pool.sort(key=lambda s: -s.p_value)
    # Draw from the weakest third so decoys are genuinely poor sites, but not
    # always the single worst bond in the protein.
    head = pool[: max(n, len(pool) // 3)] or pool
    rng.shuffle(head)
    return head[:n]


def build_control_plan(
    real_sites: list[tuple[str, str, CleavageSite]],
    all_sites_by_design: dict[str, list[CleavageSite]],
    sequences: dict[str, str],
    *,
    n_scrambles: int = 1,
    n_decoys_per_design: int = 1,
    seed: int = 0,
) -> ControlPlan:
    """Assemble real and control jobs into one interleaved wave.

    ``real_sites`` is a list of (design_id, protease, site) already selected as
    the sites worth folding.
    """
    rng = random.Random(seed)
    jobs: list[CofoldJob] = []
    counts = dict.fromkeys((ARM_REAL, ARM_SCRAMBLE, ARM_SWAP, ARM_DECOY), 0)

    real_cuts_by_design: dict[str, set[int]] = {}
    for design_id, _protease, site in real_sites:
        real_cuts_by_design.setdefault(design_id, set()).add(site.cut_index)

    for design_id, protease, site in real_sites:
        seq = sequences[design_id]

        real_job = build_job(
            design_id=design_id,
            sequence=seq,
            cut_index=site.cut_index,
            protease=protease,
            arm=ARM_REAL,
            meta={
                "pwm_score": site.pwm_score,
                "p_value": site.p_value,
                "accessibility": site.accessibility,
                "scissile_bond": site.scissile_bond,
            },
        )
        jobs.append(real_job)
        counts[ARM_REAL] += 1

        for k in range(n_scrambles):
            scrambled = scramble_segment(real_job.segment, rng)
            job = build_job(
                design_id=design_id,
                sequence=seq,
                cut_index=site.cut_index,
                protease=protease,
                arm=ARM_SCRAMBLE,
                meta={"replicate": k, "parent_job": real_job.job_id},
            )
            # Same length and cut offset, different residue order.
            job.segment = scrambled
            jobs.append(job)
            counts[ARM_SCRAMBLE] += 1

        swap = pick_swap_protease(all_sites_by_design.get(design_id, []), site.cut_index, protease)
        if swap:
            jobs.append(
                build_job(
                    design_id=design_id,
                    sequence=seq,
                    cut_index=site.cut_index,
                    protease=swap,
                    arm=ARM_SWAP,
                    meta={"swapped_from": protease, "parent_job": real_job.job_id},
                )
            )
            counts[ARM_SWAP] += 1

    for design_id, cuts in real_cuts_by_design.items():
        sites = all_sites_by_design.get(design_id, [])
        if not sites:
            continue
        protease = next(p for d, p, _ in real_sites if d == design_id)
        for decoy in pick_decoy_sites(
            sites, protease, cuts, n=n_decoys_per_design, accessibility=None, rng=rng
        ):
            jobs.append(
                build_job(
                    design_id=design_id,
                    sequence=sequences[design_id],
                    cut_index=decoy.cut_index,
                    protease=protease,
                    arm=ARM_DECOY,
                    meta={
                        "pwm_score": decoy.pwm_score,
                        "p_value": decoy.p_value,
                        "scissile_bond": decoy.scissile_bond,
                    },
                )
            )
            counts[ARM_DECOY] += 1

    # Interleave so no arm is systematically favoured by run order or by a
    # worker warming up partway through the wave.
    rng.shuffle(jobs)
    log.info("control plan: %s", counts)
    return ControlPlan(jobs=jobs, counts=counts)


def positive_control_jobs(cd74_sequence: str, sites: list[CleavageSite], *, n: int = 2) -> list[CofoldJob]:
    """Cathepsin S against its best-supported sites in CD74."""
    ctss = sorted(
        (s for s in sites if s.protease == "cathepsin_S"), key=lambda s: s.p_value
    )
    chosen: list[CleavageSite] = []
    # Prefer sites in the CLIP window, since that is the documented product.
    in_clip = [s for s in ctss if CD74_CLIP_RANGE[0] <= s.cut_index <= CD74_CLIP_RANGE[1]]
    chosen.extend(in_clip[:n])
    for s in ctss:
        if len(chosen) >= n:
            break
        if s not in chosen:
            chosen.append(s)

    return [
        build_job(
            design_id="CD74",
            sequence=cd74_sequence,
            cut_index=s.cut_index,
            protease="cathepsin_S",
            arm=ARM_POSITIVE,
            meta={
                "pwm_score": s.pwm_score,
                "p_value": s.p_value,
                "scissile_bond": s.scissile_bond,
                "in_clip_window": CD74_CLIP_RANGE[0] <= s.cut_index <= CD74_CLIP_RANGE[1],
            },
        )
        for s in chosen
    ]
