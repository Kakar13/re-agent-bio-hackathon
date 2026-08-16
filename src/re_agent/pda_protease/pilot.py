"""Pilot-set selection.

The pilot deliberately spans the range rather than sampling it. A set of designs
that all look alike cannot show that the method discriminates, so it is built
from three poles plus the reference proteins:

* **high-signal** - designs with the most significant cathepsin sites per
  residue. If the pipeline works anywhere, it works here.
* **predicted-clean** - designs with the fewest. These are the negative pole: if
  they light up downstream anyway, the method is measuring length or composition
  rather than proteolysis.
* **mid** - designs from the middle of the distribution, so the pilot is not only
  its own extremes.
* **natural anchors** - human serum albumin, clinically tolerated, and CD74,
  which carries the cathepsin S positive control.

Selection uses sequence only. Structures are fetched afterwards for the chosen
designs, so the pilot does not silently exclude anything on structural grounds.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from .cleavage import DEFAULT_P_THRESHOLD, scan_design
from .merops import SpecificityMatrix
from .pda import Design

log = logging.getLogger(__name__)

POLE_HIGH = "high-signal"
POLE_CLEAN = "predicted-clean"
POLE_MID = "mid"
POLE_ANCHOR = "natural-anchor"


@dataclass
class PilotEntry:
    design: Design
    pole: str
    site_density: float
    n_significant: int


def site_density(
    design: Design,
    matrices: dict[str, SpecificityMatrix],
    *,
    p_threshold: float = DEFAULT_P_THRESHOLD,
    n_shuffles: int = 50,
) -> tuple[float, int]:
    """Significant cathepsin sites per 100 residues, sequence only.

    A cheap screen: fewer shuffles than the full scan, and no structure, since
    this only has to rank designs rather than produce final numbers.
    """
    sites = scan_design(
        design.design_id, design.sequence, matrices, n_shuffles=n_shuffles, seed=1
    )
    bonds = {s.cut_index for s in sites if s.p_value <= p_threshold}
    n = len(bonds)
    return (100.0 * n / max(design.length, 1)), n


def select_pilot(
    pool: list[Design],
    matrices: dict[str, SpecificityMatrix],
    anchors: list[Design],
    *,
    n_high: int = 6,
    n_clean: int = 4,
    n_mid: int = 3,
    min_length: int = 60,
    max_length: int = 250,
    n_shuffles: int = 50,
) -> list[PilotEntry]:
    """Rank the pool by cleavage-site density and take from both ends plus the middle."""
    candidates = [d for d in pool if min_length <= d.length <= max_length]
    log.info("ranking %d candidate designs by cleavage-site density", len(candidates))

    scored: list[tuple[float, int, Design]] = []
    for d in candidates:
        density, n = site_density(d, matrices, n_shuffles=n_shuffles)
        scored.append((density, n, d))
    scored.sort(key=lambda x: -x[0])

    entries: list[PilotEntry] = []
    taken: set[str] = set()

    def take(items, pole: str, k: int) -> None:
        for density, n, d in items:
            if len(
                [e for e in entries if e.pole == pole]
            ) >= k:
                break
            if d.design_id in taken:
                continue
            taken.add(d.design_id)
            entries.append(PilotEntry(design=d, pole=pole, site_density=density, n_significant=n))

    take(scored, POLE_HIGH, n_high)
    take(list(reversed(scored)), POLE_CLEAN, n_clean)
    mid = len(scored) // 2
    take(scored[mid : mid + n_mid * 3], POLE_MID, n_mid)

    for a in anchors:
        density, n = site_density(a, matrices, n_shuffles=n_shuffles)
        entries.append(PilotEntry(design=a, pole=POLE_ANCHOR, site_density=density, n_significant=n))

    for e in entries:
        log.info(
            "pilot %-16s %-14s len=%-4d sites=%-3d density=%.2f/100aa",
            e.design.design_id,
            e.pole,
            e.design.length,
            e.n_significant,
            e.site_density,
        )
    return entries
