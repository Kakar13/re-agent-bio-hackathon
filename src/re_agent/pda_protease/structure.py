"""Boltz-2 co-folding of a protease with a design segment, on Modal.

This stage is deliberately downstream and deliberately ignorant. Cleavage sites
were resolved from MEROPS statistics and the design's own crystal structure;
Boltz-2 sees only two sequences and never sees the position weight matrix score.
If the co-fold places the predicted scissile bond in the catalytic cleft, that is
an independent line of evidence rather than a restatement of the input.

Protease constructs are the mature catalytic domains, sliced from UniProt by
their annotated chain boundaries: the proenzyme's propeptide blocks the active
site, so folding the full-length precursor would bury the very cleft being
measured.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .paths import BOLTZ_CACHE, ensure_caches
from .pda import fetch_uniprot_sequence

log = logging.getLogger(__name__)

SEGMENT_FLANK = 14  # residues either side of the scissile bond


@dataclass(frozen=True)
class ProteaseConstruct:
    """A mature protease domain and its catalytic machinery.

    Positions are 1-based in full-length UniProt numbering, matching how the
    catalytic residues are annotated, and converted to mature-domain indices on
    demand.
    """

    name: str
    uniprot: str
    mature_start: int
    mature_end: int
    nucleophile: int  # Cys SG, or the catalytic Asp of an aspartic peptidase
    nucleophile_atom: str
    catalytic_other: dict[str, int] = field(default_factory=dict)
    family: str = ""
    note: str = ""

    def sequence(self) -> str:
        full = fetch_uniprot_sequence(self.uniprot)
        return full[self.mature_start - 1 : self.mature_end]

    def mature_index(self, full_position: int) -> int:
        """0-based index of a full-length position within the mature chain."""
        return full_position - self.mature_start

    @property
    def nucleophile_index(self) -> int:
        return self.mature_index(self.nucleophile)


# Cysteine cathepsins share the papain fold: Cys-His-Asn triad, nucleophile is
# the Cys thiolate, which attacks the P1 carbonyl carbon. Cathepsin D is an
# aspartic peptidase and is excluded from co-folding because its mature form is
# a two-chain heterodimer that a single-chain co-fold would misrepresent.
PROTEASES: dict[str, ProteaseConstruct] = {
    "cathepsin_S": ProteaseConstruct(
        name="cathepsin_S",
        uniprot="P25774",
        mature_start=115,
        mature_end=331,
        nucleophile=139,
        nucleophile_atom="SG",
        catalytic_other={"His": 278, "Asn": 298},
        family="papain (Cys)",
        note="endosomal; generates CLIP from CD74",
    ),
    "cathepsin_L": ProteaseConstruct(
        name="cathepsin_L",
        uniprot="P07711",
        mature_start=114,
        mature_end=333,
        nucleophile=138,
        nucleophile_atom="SG",
        catalytic_other={"His": 276, "Asn": 300},
        family="papain (Cys)",
        note="mature form is nicked into heavy 114-288 and light 292-333",
    ),
    "cathepsin_K": ProteaseConstruct(
        name="cathepsin_K",
        uniprot="P43235",
        mature_start=115,
        mature_end=329,
        nucleophile=139,
        nucleophile_atom="SG",
        catalytic_other={"His": 276, "Asn": 296},
        family="papain (Cys)",
    ),
    "cathepsin_B": ProteaseConstruct(
        name="cathepsin_B",
        uniprot="P07858",
        mature_start=80,
        mature_end=333,
        nucleophile=108,
        nucleophile_atom="SG",
        catalytic_other={"His": 278, "Asn": 298},
        family="papain (Cys)",
        note="occluding loop gives it carboxydipeptidase character",
    ),
    "legumain": ProteaseConstruct(
        name="legumain",
        uniprot="Q99538",
        mature_start=18,
        mature_end=323,
        nucleophile=189,
        nucleophile_atom="SG",
        catalytic_other={"His": 148},
        family="legumain (Cys)",
        note="strict Asn/Asp P1 specificity",
    ),
}

COFOLDABLE = tuple(PROTEASES)


@dataclass
class CofoldJob:
    """One protease-plus-segment complex to predict."""

    job_id: str
    protease: str
    segment: str
    segment_start: int  # 0-based offset of the segment in the parent sequence
    cut_offset: int  # index within the segment of the residue after the bond
    design_id: str
    arm: str  # real | scramble | protease_swap | decoy | positive_control
    source_cut_index: int = -1
    meta: dict = field(default_factory=dict)

    @property
    def p1_offset(self) -> int:
        """Index within the segment of the P1 residue, which donates the
        carbonyl carbon the nucleophile attacks."""
        return self.cut_offset - 1

    def to_dict(self) -> dict:
        return asdict(self)


def make_segment(sequence: str, cut_index: int, flank: int = SEGMENT_FLANK) -> tuple[str, int, int]:
    """Cut out a segment centred on a scissile bond.

    Returns the segment, its start offset in the parent, and the index of the
    first residue after the bond within the segment.
    """
    start = max(0, cut_index - flank)
    end = min(len(sequence), cut_index + flank)
    return sequence[start:end], start, cut_index - start


def build_job(
    *,
    design_id: str,
    sequence: str,
    cut_index: int,
    protease: str,
    arm: str = "real",
    flank: int = SEGMENT_FLANK,
    meta: dict | None = None,
) -> CofoldJob:
    segment, seg_start, cut_off = make_segment(sequence, cut_index, flank)
    raw = f"{design_id}|{protease}|{cut_index}|{arm}|{segment}"
    jid = f"{arm}_{protease}_{design_id}_{cut_index}_{hashlib.sha1(raw.encode()).hexdigest()[:6]}"
    return CofoldJob(
        job_id=jid,
        protease=protease,
        segment=segment,
        segment_start=seg_start,
        cut_offset=cut_off,
        design_id=design_id,
        arm=arm,
        source_cut_index=cut_index,
        meta=meta or {},
    )


def _cache_path(job: CofoldJob) -> Path:
    ensure_caches()
    return BOLTZ_CACHE / f"{job.job_id}.json"


def run_cofold(
    job: CofoldJob,
    *,
    use_msa: bool = True,
    recycling_steps: int | None = None,
    sampling_steps: int | None = None,
    timeout: int = 3600,
    force: bool = False,
) -> dict | None:
    """Predict one complex, caching the coordinates and confidence metrics."""
    cache = _cache_path(job)
    if cache.exists() and not force:
        log.debug("boltz cache hit %s", job.job_id)
        return json.loads(cache.read_text())

    from proto_tools.modal.client import dispatch_to_modal
    from proto_tools.tools.structure_prediction.boltz2.boltz2 import Boltz2Config, Boltz2Input

    protease_seq = PROTEASES[job.protease].sequence()
    chains = [
        {"sequence": protease_seq, "entity_type": "protein"},
        {"sequence": job.segment, "entity_type": "protein"},
    ]
    cfg: dict = {"use_msa": use_msa, "timeout": timeout, "verbose": 1}
    if recycling_steps is not None:
        cfg["recycling_steps"] = recycling_steps
    if sampling_steps is not None:
        cfg["sampling_steps"] = sampling_steps

    log.info("boltz %s (%d + %d aa)", job.job_id, len(protease_seq), len(job.segment))
    # Dispatch to the deployed Modal app. Calling run_boltz2 directly would build
    # a local conda environment and run on this laptop's CPU instead.
    out = dispatch_to_modal(
        "boltz2-prediction",
        Boltz2Input(complexes=[{"chains": chains}]),
        Boltz2Config(**cfg),
    )

    if not out.structures:
        log.warning("boltz returned no structure for %s: %s", job.job_id, out.errors)
        return None

    st = out.structures[0]
    m = st.metrics
    record = {
        "job": job.to_dict(),
        "protease_length": len(protease_seq),
        "structure_format": getattr(st, "structure_format", "cif"),
        "structure": st.structure,
        "metrics": {
            k: getattr(m, k, None)
            for k in (
                "confidence_score",
                "ptm",
                "iptm",
                "complex_plddt",
                "complex_iplddt",
                "avg_pae",
                "pair_chains_iptm",
            )
        },
        "success": bool(out.success),
        "errors": list(out.errors or []),
    }
    cache.write_text(json.dumps(record))
    return record


def run_wave(jobs: list[CofoldJob], **kwargs) -> dict[str, dict]:
    """Run a batch against the deployed Modal app.

    The container stays warm between calls, so the model loads once for the wave
    rather than once per job.
    """
    results: dict[str, dict] = {}
    pending = [j for j in jobs if not _cache_path(j).exists()]
    log.info("boltz wave: %d jobs, %d already cached", len(jobs), len(jobs) - len(pending))

    for i, job in enumerate(jobs, 1):
        try:
            rec = run_cofold(job, **kwargs)
            if rec:
                results[job.job_id] = rec
        except Exception as exc:  # noqa: BLE001 - one failure must not sink the wave
            log.warning("boltz job %s failed: %s", job.job_id, exc)
        log.info("boltz wave progress %d/%d", i, len(jobs))
    return results
