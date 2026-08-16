"""Measure whether a co-fold actually poses a bond for catalysis.

A high interface confidence only says two chains stick together; it does not say
they are arranged as an enzyme and its substrate. The discriminating question is
geometric: does the catalytic nucleophile sit against the carbonyl carbon of the
bond that the sequence model predicted would be cut, with the flanking residues
seated in the subsites that select them.

Three measurements, in decreasing order of how hard they are to satisfy by
accident:

* distance from the nucleophile to the scissile carbonyl carbon,
* burial of the P1 and P2 side chains against the enzyme (the S1/S2 subsites),
* interface confidence, ipTM.

Thresholds are deliberately loose. A co-folded model is not a crystal structure,
and the analysis compares real sites against controls rather than against an
absolute standard.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from io import StringIO

import numpy as np

from .structure import PROTEASES, CofoldJob

log = logging.getLogger(__name__)

# A thiolate poised to attack sits ~3 A away; a co-fold that places it within
# this range has the bond in the cleft rather than merely nearby.
ENGAGED_DISTANCE = 8.0
POSED_DISTANCE = 5.0
CONTACT_CUTOFF = 4.5

BACKBONE = {"N", "CA", "C", "O", "OXT"}


@dataclass
class Geometry:
    """Active-site geometry for one predicted complex."""

    job_id: str
    protease: str
    design_id: str
    arm: str
    nucleophile_distance: float | None
    p1_residue: str
    p2_residue: str
    p1_contacts: int
    p2_contacts: int
    min_interface_distance: float | None
    iptm: float | None
    complex_plddt: float | None
    interface_plddt: float | None
    engaged: bool
    posed: bool

    def to_dict(self) -> dict:
        return asdict(self)


def _parse(text: str, fmt: str):
    import biotite.structure.io.pdb as pdb
    import biotite.structure.io.pdbx as pdbx

    handle = StringIO(text)
    if (fmt or "cif").lower() == "cif":
        return pdbx.get_structure(pdbx.CIFFile.read(handle), model=1, extra_fields=["b_factor"])
    return pdb.get_structure(pdb.PDBFile.read(handle), model=1, extra_fields=["b_factor"])


def _split_chains(arr, protease_length: int):
    """Identify which chain is the enzyme by length, not by chain letter."""
    import biotite.structure as struc

    best = None
    for cid in struc.get_chains(arr):
        sub = arr[arr.chain_id == cid]
        n_res = struc.get_residue_count(sub)
        delta = abs(n_res - protease_length)
        if best is None or delta < best[0]:
            best = (delta, cid, n_res)
    if best is None:
        return None, None
    _, prot_chain, _ = best
    protease = arr[arr.chain_id == prot_chain]
    substrate = arr[arr.chain_id != prot_chain]
    return protease, substrate


def _residue_atoms(chain, residue_offset: int):
    """Atoms of the nth residue in a chain, by order rather than by number."""
    import biotite.structure as struc

    res_ids = struc.get_residues(chain)[0]
    if residue_offset < 0 or residue_offset >= len(res_ids):
        return None
    return chain[chain.res_id == res_ids[residue_offset]]


def analyse(record: dict) -> Geometry | None:
    """Measure one cached Boltz record."""
    job = CofoldJob(**record["job"])
    construct = PROTEASES.get(job.protease)
    if construct is None:
        return None

    try:
        arr = _parse(record["structure"], record.get("structure_format", "cif"))
    except Exception as exc:  # noqa: BLE001
        log.warning("could not parse structure for %s: %s", job.job_id, exc)
        return None

    import biotite.structure as struc

    arr = arr[struc.filter_amino_acids(arr)]
    protease, substrate = _split_chains(arr, record.get("protease_length", 0))
    if protease is None or substrate is None or substrate.array_length() == 0:
        log.warning("could not split chains for %s", job.job_id)
        return None

    metrics = record.get("metrics", {})

    # Catalytic nucleophile, and the carbonyl carbon it would attack.
    nuc_atoms = _residue_atoms(protease, construct.nucleophile_index)
    p1_atoms = _residue_atoms(substrate, job.p1_offset)
    p2_atoms = _residue_atoms(substrate, job.p1_offset - 1)

    distance = None
    if nuc_atoms is not None and p1_atoms is not None:
        sel = nuc_atoms[nuc_atoms.atom_name == construct.nucleophile_atom]
        carbonyl = p1_atoms[p1_atoms.atom_name == "C"]
        if sel.array_length() and carbonyl.array_length():
            distance = float(np.linalg.norm(sel.coord[0] - carbonyl.coord[0]))

    def side_chain_contacts(atoms) -> int:
        """Protease heavy atoms packed against a substrate side chain."""
        if atoms is None or protease.array_length() == 0:
            return 0
        sc = atoms[~np.isin(atoms.atom_name, list(BACKBONE))]
        sc = sc[sc.element != "H"]
        if sc.array_length() == 0:
            return 0
        heavy = protease[protease.element != "H"]
        d = np.linalg.norm(heavy.coord[:, None, :] - sc.coord[None, :, :], axis=-1)
        return int((d < CONTACT_CUTOFF).sum())

    p1_res = str(p1_atoms.res_name[0]) if p1_atoms is not None and p1_atoms.array_length() else "?"
    p2_res = str(p2_atoms.res_name[0]) if p2_atoms is not None and p2_atoms.array_length() else "?"

    heavy_p = protease[protease.element != "H"]
    heavy_s = substrate[substrate.element != "H"]
    min_iface = None
    if heavy_p.array_length() and heavy_s.array_length():
        d = np.linalg.norm(heavy_p.coord[:, None, :] - heavy_s.coord[None, :, :], axis=-1)
        min_iface = float(d.min())

    iface_plddt = None
    if hasattr(substrate, "b_factor") and substrate.array_length():
        # Boltz writes per-atom pLDDT into the B-factor column.
        iface_plddt = float(np.nanmean(substrate.b_factor))

    return Geometry(
        job_id=job.job_id,
        protease=job.protease,
        design_id=job.design_id,
        arm=job.arm,
        nucleophile_distance=distance,
        p1_residue=p1_res,
        p2_residue=p2_res,
        p1_contacts=side_chain_contacts(p1_atoms),
        p2_contacts=side_chain_contacts(p2_atoms),
        min_interface_distance=min_iface,
        iptm=metrics.get("iptm"),
        complex_plddt=metrics.get("complex_plddt"),
        interface_plddt=iface_plddt,
        engaged=bool(distance is not None and distance <= ENGAGED_DISTANCE),
        posed=bool(distance is not None and distance <= POSED_DISTANCE),
    )


def analyse_all(records: dict[str, dict]) -> list[Geometry]:
    out = []
    for rec in records.values():
        g = analyse(rec)
        if g is not None:
            out.append(g)
    return out
