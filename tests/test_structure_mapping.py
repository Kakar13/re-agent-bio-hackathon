from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from re_agent.immuno.structure import structure_reference_from_pdb


def _atom(serial: int, residue: str, chain: str, number: int) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{number:4d}    "
        "   0.000   0.000   0.000  1.00 20.00           C"
    )


class StructureMappingTests(unittest.TestCase):
    def test_exact_chain_sequence_produces_residue_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "candidate.pdb"
            path.write_text(
                "\n".join(
                    [
                        _atom(1, "ALA", "A", 7),
                        _atom(2, "CYS", "A", 8),
                        _atom(3, "ASP", "A", 9),
                        "END",
                    ]
                )
                + "\n"
            )

            reference = structure_reference_from_pdb(
                path,
                sequence="ACD",
                chain_id="A",
                repository_root=root,
            )

        self.assertEqual(reference.residue_ids, ["7", "8", "9"])
        self.assertEqual(reference.mapping_status, "verified_exact_sequence")
        self.assertEqual(reference.path, "candidate.pdb")

    def test_unresolved_terminal_residues_preserve_full_sequence_alignment(self):
        residues = [
            "ALA",
            "CYS",
            "ASP",
            "GLU",
            "PHE",
            "GLY",
            "HIS",
            "ILE",
            "LYS",
            "LEU",
            "MET",
            "ASN",
            "PRO",
            "GLN",
            "ARG",
            "SER",
            "THR",
            "VAL",
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "terminal-trim.pdb"
            path.write_text(
                "\n".join(
                    _atom(index, residue, "A", index)
                    for index, residue in enumerate(residues, start=1)
                )
                + "\n"
            )

            reference = structure_reference_from_pdb(
                path,
                sequence="GACDEFGHIKLMNPQRSTVP",
                chain_id="A",
                repository_root=root,
            )

        self.assertEqual(reference.mapping_status, "verified_terminal_trim")
        self.assertEqual(reference.unresolved_sequence_positions, [1, 20])
        self.assertEqual(reference.residue_ids, ["", *map(str, range(1, 19)), ""])

    def test_sequence_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "candidate.pdb"
            path.write_text(_atom(1, "ALA", "A", 1) + "\n")

            with self.assertRaisesRegex(ValueError, "does not exactly match"):
                structure_reference_from_pdb(
                    path,
                    sequence="AC",
                    chain_id="A",
                    repository_root=root,
                )


if __name__ == "__main__":
    unittest.main()
