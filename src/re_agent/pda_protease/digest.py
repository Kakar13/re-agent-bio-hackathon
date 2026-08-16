"""Turn resolved cleavage sites into the peptides that could actually exist.

This is the step that makes protease attribution meaningful. A conventional
scan asks "which 15-mers in this sequence bind MHC class II"; every window is
assumed to exist. Here a peptide only exists if proteolysis could have produced
it, and it carries the identity of the protease that produced each terminus.

The two pathways differ, so they are digested differently:

* **Class I** - the proteasome sets the C-terminus, then ERAP1 trims the
  N-terminus in the endoplasmic reticulum. So the C-terminus must be a NetChop
  site while the N-terminus is left ragged.
* **Class II** - endolysosomal cathepsins cut the antigen into fragments and the
  open-ended class II groove binds within them. So a 15-mer may not cross a
  high-confidence internal cut, because that peptide would already be severed.

``unconstrained_peptides`` produces the naive every-window set, kept as the
diagnostic baseline the constrained set is measured against.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field

from .cleavage import CleavageSite

log = logging.getLogger(__name__)

MHCI_LENGTHS = (8, 9, 10, 11)
MHCII_LENGTH = 15
MHCII_MIN_FRAGMENT = 11  # shortest fragment still able to fill the groove

ARM_MHCI = "mhc_i"
ARM_MHCII = "mhc_ii"


@dataclass
class Peptide:
    """A peptide with the proteolytic history that produced it."""

    design_id: str
    peptide: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive
    arm: str
    n_term_source: str  # protease or mechanism setting the N-terminus
    c_term_source: str
    constrained: bool = True
    fragment_id: str = ""
    supporting_sites: list[int] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.peptide)

    @property
    def key(self) -> str:
        return f"{self.design_id}:{self.start}-{self.end}:{self.arm}"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["length"] = self.length
        return d


@dataclass
class Fragment:
    """A stretch of sequence bounded by two cathepsin cuts."""

    design_id: str
    start: int
    end: int
    sequence: str
    n_term_protease: str
    c_term_protease: str

    @property
    def fragment_id(self) -> str:
        return f"{self.design_id}:{self.start}-{self.end}"


def _protease_at(sites: list[CleavageSite], cut_index: int) -> str:
    """Best-supported protease for a bond, for provenance labelling."""
    best, best_p = "unknown", 2.0
    for s in sites:
        if s.cut_index == cut_index and s.p_value < best_p:
            best, best_p = s.protease, s.p_value
    return best


# ------------------------------------------------------------------ class I


def digest_mhci(
    design_id: str,
    sequence: str,
    netchop_cuts: set[int],
    *,
    lengths: tuple[int, ...] = MHCI_LENGTHS,
) -> list[Peptide]:
    """Peptides whose C-terminus is a predicted proteasome cut.

    ``netchop_cuts`` holds 0-based indices of residues whose C-terminal bond is
    cleaved, so a peptide ending at that residue is inclusive of it.
    """
    out: list[Peptide] = []
    for c in sorted(netchop_cuts):
        end = c + 1
        for L in lengths:
            start = end - L
            if start < 0:
                continue
            out.append(
                Peptide(
                    design_id=design_id,
                    peptide=sequence[start:end],
                    start=start,
                    end=end,
                    arm=ARM_MHCI,
                    n_term_source="ERAP1-ragged",
                    c_term_source="proteasome",
                    supporting_sites=[c],
                )
            )
    log.debug("%s: %d class I peptides from %d cuts", design_id, len(out), len(netchop_cuts))
    return out


# ----------------------------------------------------------------- class II


def build_fragments(
    design_id: str,
    sequence: str,
    cut_indices: set[int],
    sites: list[CleavageSite],
) -> list[Fragment]:
    """Split a sequence at high-confidence cathepsin cuts."""
    bounds = sorted(c for c in cut_indices if 0 < c < len(sequence))
    edges = [0, *bounds, len(sequence)]
    frags: list[Fragment] = []
    for i in range(len(edges) - 1):
        s, e = edges[i], edges[i + 1]
        if e - s <= 0:
            continue
        frags.append(
            Fragment(
                design_id=design_id,
                start=s,
                end=e,
                sequence=sequence[s:e],
                n_term_protease=("N-terminus" if s == 0 else _protease_at(sites, s)),
                c_term_protease=("C-terminus" if e == len(sequence) else _protease_at(sites, e)),
            )
        )
    return frags


def digest_mhcii(
    design_id: str,
    sequence: str,
    cut_indices: set[int],
    sites: list[CleavageSite],
    *,
    length: int = MHCII_LENGTH,
    min_fragment: int = MHCII_MIN_FRAGMENT,
) -> list[Peptide]:
    """15-mers confined within cathepsin fragments.

    A fragment shorter than the nominal length is still presentable if it can
    fill the groove, so it is emitted whole rather than discarded.
    """
    out: list[Peptide] = []
    for frag in build_fragments(design_id, sequence, cut_indices, sites):
        flen = frag.end - frag.start
        if flen < min_fragment:
            continue
        if flen < length:
            out.append(
                Peptide(
                    design_id=design_id,
                    peptide=frag.sequence,
                    start=frag.start,
                    end=frag.end,
                    arm=ARM_MHCII,
                    n_term_source=frag.n_term_protease,
                    c_term_source=frag.c_term_protease,
                    fragment_id=frag.fragment_id,
                    supporting_sites=[frag.start, frag.end],
                )
            )
            continue
        for off in range(flen - length + 1):
            s = frag.start + off
            e = s + length
            out.append(
                Peptide(
                    design_id=design_id,
                    peptide=sequence[s:e],
                    start=s,
                    end=e,
                    arm=ARM_MHCII,
                    # Only a terminus that coincides with a fragment edge was set
                    # by a protease; an interior offset is just a register shift.
                    n_term_source=(frag.n_term_protease if s == frag.start else "internal"),
                    c_term_source=(frag.c_term_protease if e == frag.end else "internal"),
                    fragment_id=frag.fragment_id,
                    supporting_sites=[frag.start, frag.end],
                )
            )
    log.debug("%s: %d class II peptides from %d cuts", design_id, len(out), len(cut_indices))
    return out


def unconstrained_peptides(
    design_id: str, sequence: str, *, length: int = MHCII_LENGTH, arm: str = ARM_MHCII
) -> list[Peptide]:
    """Every window, ignoring proteolysis. The baseline, not a result."""
    return [
        Peptide(
            design_id=design_id,
            peptide=sequence[i : i + length],
            start=i,
            end=i + length,
            arm=arm,
            n_term_source="none",
            c_term_source="none",
            constrained=False,
        )
        for i in range(len(sequence) - length + 1)
    ]


def dedupe(peptides: list[Peptide]) -> list[Peptide]:
    """Collapse peptides identical in sequence, arm and position."""
    seen: dict[str, Peptide] = {}
    for p in peptides:
        seen.setdefault(p.key, p)
    return list(seen.values())
