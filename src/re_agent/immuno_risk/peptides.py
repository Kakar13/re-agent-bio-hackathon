"""Peptide window helpers and simple physicochemical descriptors."""

from __future__ import annotations

AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = set("AILMFVWY")
POSITIVE = set("KRH")
NEGATIVE = set("DE")


def clean_sequence(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in AA)


def sliding_windows(sequence: str, lengths: range | list[int]) -> list[tuple[int, int, str]]:
    seq = clean_sequence(sequence)
    out: list[tuple[int, int, str]] = []
    for L in lengths:
        if L > len(seq):
            continue
        for i in range(0, len(seq) - L + 1):
            out.append((i, i + L, seq[i : i + L]))
    return out


def peptide_descriptors(peptide: str) -> dict[str, float]:
    p = clean_sequence(peptide)
    n = max(len(p), 1)
    return {
        "length": float(len(p)),
        "hydrophobic_fraction": sum(1 for c in p if c in HYDROPHOBIC) / n,
        "charge_proxy": (sum(1 for c in p if c in POSITIVE) - sum(1 for c in p if c in NEGATIVE))
        / n,
        "aromatic_fraction": sum(1 for c in p if c in "FWY") / n,
        "cysteine_fraction": sum(1 for c in p if c == "C") / n,
    }


def net_charge_ph74(sequence: str) -> float:
    seq = clean_sequence(sequence)
    return float(sum(1 for c in seq if c in POSITIVE) - sum(1 for c in seq if c in NEGATIVE))
