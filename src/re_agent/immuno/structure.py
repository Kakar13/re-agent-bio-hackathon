"""Strict sequence-to-PDB mapping for residue-level evidence visualization."""

from __future__ import annotations

import hashlib
from pathlib import Path

from re_agent.immuno.contracts import StructureReference

_THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _pdb_chain_sequence(pdb_text: str, chain_id: str) -> tuple[str, list[str]]:
    residues: list[str] = []
    residue_ids: list[str] = []
    seen: set[tuple[str, str]] = set()

    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            continue
        if line[21].strip() != chain_id:
            continue
        residue_name = line[17:20].strip().upper()
        amino_acid = _THREE_TO_ONE.get(residue_name)
        if amino_acid is None:
            continue
        sequence_number = line[22:26].strip()
        insertion_code = line[26].strip()
        key = (sequence_number, insertion_code)
        if key in seen:
            continue
        seen.add(key)
        residues.append(amino_acid)
        residue_ids.append(f"{sequence_number}{insertion_code}")

    return "".join(residues), residue_ids


def structure_reference_from_pdb(
    path: Path,
    *,
    sequence: str,
    chain_id: str,
    repository_root: Path,
) -> StructureReference:
    """Create a reference only when one PDB chain exactly matches the screened sequence."""

    resolved = path.resolve()
    root = repository_root.resolve()
    if resolved != root and root not in resolved.parents:
        raise ValueError("structure path must be inside the repository")
    if resolved.suffix.lower() != ".pdb":
        raise ValueError("3D evidence overlays currently require a PDB file")
    if not resolved.is_file():
        raise FileNotFoundError(f"structure file does not exist: {resolved}")

    pdb_bytes = resolved.read_bytes()
    pdb_text = pdb_bytes.decode("utf-8")
    structure_sequence, residue_ids = _pdb_chain_sequence(pdb_text, chain_id)
    if not structure_sequence:
        raise ValueError(f"structure contains no canonical residues for chain {chain_id!r}")
    mapping_status = "verified_exact_sequence"
    unresolved_positions: list[int] = []
    if structure_sequence != sequence:
        start = sequence.find(structure_sequence)
        terminal_trim_is_safe = (
            start >= 0
            and len(structure_sequence) >= 15
            and len(structure_sequence) / len(sequence) >= 0.9
            and start <= 5
            and len(sequence) - start - len(structure_sequence) <= 5
        )
        if not terminal_trim_is_safe:
            raise ValueError(
                "structure chain sequence does not exactly match the screened sequence: "
                f"chain {chain_id!r} has {len(structure_sequence)} residues, "
                f"screened sequence has {len(sequence)}"
            )
        end = start + len(structure_sequence)
        unresolved_positions = [
            *range(1, start + 1),
            *range(end + 1, len(sequence) + 1),
        ]
        residue_ids = [
            *([""] * start),
            *residue_ids,
            *([""] * (len(sequence) - end)),
        ]
        mapping_status = "verified_terminal_trim"

    return StructureReference(
        path=str(resolved.relative_to(root)),
        format="pdb",
        chain_id=chain_id,
        residue_ids=residue_ids,
        unresolved_sequence_positions=unresolved_positions,
        sequence_sha256=_sha256_bytes(sequence.encode()),
        structure_sha256=_sha256_bytes(pdb_bytes),
        mapping_status=mapping_status,
    )
