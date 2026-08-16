"""Structural gate on cleavage sites, from experimental coordinates.

A protease cannot cut a bond it cannot reach. Threading a substrate through a
cathepsin cleft requires roughly eight residues in an extended conformation, so
accessibility is scored over the whole P4-P4' span rather than at one residue.

Two things keep this honest:

* Coordinates come from the design's own solved structure, never from a
  prediction, which is what lets Boltz stay an independent downstream check.
* Buried sites are classified, not discarded. In the endosome GILT reduces
  disulfides and the antigen partially unfolds, so a site buried in the crystal
  can still be cleaved. Reporting the split between reachable-as-folded and
  unfolding-required is more informative than a hard filter, especially for
  hyperstable de novo designs.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

# Theoretical maximum solvent accessibility per residue, Tien et al. 2013
# (PLoS ONE 8:e80635). Used to turn absolute SASA into a 0-1 relative value.
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}

# A span is treated as reachable without unfolding when it is meaningfully
# solvent-exposed and not locked inside regular secondary structure.
RSA_ACCESSIBLE = 0.25
COIL_FRACTION_ACCESSIBLE = 0.5

NATIVE_ACCESSIBLE = "native-accessible"
UNFOLDING_REQUIRED = "unfolding-required"
UNOBSERVED = "unobserved"


@dataclass
class ResidueStructure:
    """Per-residue structural features, indexed against the design sequence."""

    seq_index: int
    residue: str
    observed: bool
    rsa: float | None = None
    bfactor: float | None = None
    sse: str = "-"  # a=helix, b=strand, c=coil, -=unobserved


@dataclass
class SpanFeatures:
    """Structural summary of one P4-P4' span."""

    mean_rsa: float | None
    mean_bfactor: float | None
    coil_fraction: float
    observed_fraction: float
    classification: str


def _load_chain(cif_path: Path, chain_id: str):
    import biotite.structure as struc
    import biotite.structure.io.pdbx as pdbx

    f = pdbx.CIFFile.read(str(cif_path))
    # Keep B-factors; drop waters/ligands and alternate models.
    arr = pdbx.get_structure(f, model=1, extra_fields=["b_factor"])
    arr = arr[struc.filter_amino_acids(arr)]
    arr = arr[~arr.hetero | np.isin(arr.res_name, ["MSE"])]

    for field in ("label_asym_id", "chain_id"):
        if field in arr.get_annotation_categories():
            sel = arr[getattr(arr, field) == chain_id]
            if sel.array_length() > 0:
                return sel
    # Fall back to the first chain present rather than failing outright.
    chains = struc.get_chains(arr)
    if len(chains) == 0:
        return None
    log.debug("chain %s not found in %s; using %s", chain_id, cif_path.name, chains[0])
    return arr[arr.chain_id == chains[0]]


def _observed_sequence(chain) -> tuple[str, list[int]]:
    import biotite.structure as struc

    res_ids = struc.get_residues(chain)[0]
    res_names = struc.get_residues(chain)[1]
    seq = "".join(THREE_TO_ONE.get(rn, "X") for rn in res_names)
    return seq, list(res_ids)


def _align_to_design(observed: str, design: str) -> dict[int, int]:
    """Map observed-residue order onto design-sequence indices.

    Crystal structures routinely miss disordered loops and may carry expression
    tags, so positions cannot be assumed to line up. Exact substring match covers
    the common case cheaply; otherwise fall back to a global alignment.
    """
    if not observed:
        return {}
    pos = design.find(observed)
    if pos >= 0:
        return {i: pos + i for i in range(len(observed))}

    from Bio import Align

    aligner = Align.PairwiseAligner(
        mode="global", open_gap_score=-10, extend_gap_score=-0.5, match_score=2, mismatch_score=-1
    )
    try:
        aln = aligner.align(design, observed)[0]
    except Exception:  # noqa: BLE001
        return {}
    mapping: dict[int, int] = {}
    for (d0, d1), (o0, o1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(d1 - d0):
            mapping[o0 + k] = d0 + k
    return mapping


def compute_structure_features(
    cif_path: Path, chain_id: str, design_sequence: str
) -> list[ResidueStructure]:
    """Per-residue RSA, B-factor and secondary structure, aligned to the sequence."""
    import biotite.structure as struc

    features = [
        ResidueStructure(seq_index=i, residue=aa, observed=False)
        for i, aa in enumerate(design_sequence)
    ]

    chain = _load_chain(cif_path, chain_id)
    if chain is None or chain.array_length() == 0:
        log.warning("no usable chain in %s", cif_path.name)
        return features

    observed_seq, _ = _observed_sequence(chain)
    mapping = _align_to_design(observed_seq, design_sequence)
    if not mapping:
        log.warning("could not align %s chain %s to its sequence", cif_path.name, chain_id)
        return features

    try:
        atom_sasa = struc.sasa(chain, vdw_radii="Single")
        res_sasa = struc.apply_residue_wise(chain, atom_sasa, np.nansum)
    except Exception as exc:  # noqa: BLE001
        log.warning("SASA failed for %s: %s", cif_path.name, exc)
        return features

    try:
        sse = struc.annotate_sse(chain)
    except Exception:  # noqa: BLE001
        sse = np.array(["c"] * len(res_sasa))

    b_atoms = getattr(chain, "b_factor", None)
    res_b = (
        struc.apply_residue_wise(chain, b_atoms, np.nanmean)
        if b_atoms is not None
        else np.full(len(res_sasa), np.nan)
    )

    for obs_i, design_i in mapping.items():
        if obs_i >= len(res_sasa) or design_i >= len(features):
            continue
        aa = design_sequence[design_i]
        max_asa = MAX_ASA.get(aa)
        rsa = float(res_sasa[obs_i]) / max_asa if max_asa else None
        f = features[design_i]
        f.observed = True
        f.rsa = min(rsa, 1.5) if rsa is not None else None
        b = float(res_b[obs_i]) if obs_i < len(res_b) else float("nan")
        f.bfactor = None if np.isnan(b) else b
        f.sse = str(sse[obs_i]) if obs_i < len(sse) else "c"

    n_obs = sum(1 for f in features if f.observed)
    log.debug("%s: %d/%d residues observed", cif_path.name, n_obs, len(features))
    return features


def span_features(
    features: list[ResidueStructure], cut_index: int, *, flank: int = 4
) -> SpanFeatures:
    """Summarise the P4-P4' span around a bond before ``cut_index``."""
    start, end = cut_index - flank, cut_index + flank
    if start < 0 or end > len(features):
        return SpanFeatures(None, None, 0.0, 0.0, UNOBSERVED)

    span = features[start:end]
    observed = [f for f in span if f.observed]
    observed_fraction = len(observed) / len(span)
    if not observed:
        # Absent from the density map means disordered, and disordered loops are
        # the classic protease-labile sites. Treat it as reachable rather than
        # unknown, but record that no coordinates backed the call.
        return SpanFeatures(None, None, 1.0, 0.0, NATIVE_ACCESSIBLE)

    rsas = [f.rsa for f in observed if f.rsa is not None]
    bs = [f.bfactor for f in observed if f.bfactor is not None]
    coil_fraction = sum(1 for f in observed if f.sse not in ("a", "b")) / len(observed)

    mean_rsa = float(np.mean(rsas)) if rsas else None
    mean_b = float(np.mean(bs)) if bs else None

    # Unobserved residues are disordered, which means flexible and reachable.
    if observed_fraction < 0.5:
        classification = NATIVE_ACCESSIBLE
    elif mean_rsa is None:
        classification = UNOBSERVED
    elif mean_rsa >= RSA_ACCESSIBLE or coil_fraction >= COIL_FRACTION_ACCESSIBLE:
        classification = NATIVE_ACCESSIBLE
    else:
        classification = UNFOLDING_REQUIRED

    return SpanFeatures(mean_rsa, mean_b, coil_fraction, observed_fraction, classification)
