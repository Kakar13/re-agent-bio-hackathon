#!/usr/bin/env python3
"""Join per-residue solvent accessibility onto the PDA MHC-I epitope profile.

Answers a design question, not an immunological one: of the 9-mers flagged as
MHC-I risks, how many sit on the surface where ProteinMPNN can resample them
freely, and how many are buried in the core where changing them threatens the
fold?

Structures come from RCSB for the PDA parents themselves, which are deposited
entries, so no design-time structure prediction is needed.
"""

from __future__ import annotations

import argparse
import gzip
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.e2e_pls.accessibility import (
    BURIED_RSA_THRESHOLD,
    align_chain_to_design,
    project_onto_design,
    relative_accessibility,
    window_accessibility,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROFILE = (
    REPO_ROOT / "data/processed/profiles/a0201-pda-mhci-profile-v4/pda_mhci_profile.parquet"
)
DEFAULT_DESIGNS = REPO_ROOT / "data/processed/pda_designs.parquet"
DEFAULT_CACHE = REPO_ROOT / "data/raw/rcsb_cif"
DEFAULT_OUTPUT = REPO_ROOT / "results/benchmarks/epitope_accessibility"
RCSB_TEMPLATE = "https://files.rcsb.org/download/{code}.cif.gz"
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def fetch_structure(code: str, cache_dir: Path) -> Path | None:
    """Download one mmCIF, caching so reruns cost nothing."""

    import httpx

    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{code.lower()}.cif.gz"
    if path.exists() and path.stat().st_size > 0:
        return path
    try:
        response = httpx.get(RCSB_TEMPLATE.format(code=code.lower()), timeout=60.0)
        response.raise_for_status()
    except Exception as error:  # noqa: BLE001 - a missing entry must not abort the sweep
        print(f"  fetch failed for {code}: {type(error).__name__}", flush=True)
        return None
    path.write_bytes(response.content)
    return path


def chain_accessibility(path: Path) -> dict[str, tuple[str, np.ndarray]]:
    """Return observed sequence and relative SASA for every protein chain."""

    from Bio.PDB import MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley

    parser = MMCIFParser(QUIET=True)
    with gzip.open(path, "rt") as handle, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure(path.stem, handle)

    model = next(iter(structure))
    # SASA is computed on the whole deposited model, so residues buried at a
    # crystallographic interface read as buried. That is the right answer for
    # an assembly and a caveat for a monomeric design.
    ShrakeRupley().compute(model, level="R")

    chains: dict[str, tuple[str, np.ndarray]] = {}
    for chain in model:
        codes: list[str] = []
        values: list[float] = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            letter = THREE_TO_ONE.get(residue.get_resname().upper())
            if letter is None:
                continue
            codes.append(letter)
            values.append(float(getattr(residue, "sasa", np.nan)))
        if len(codes) >= 20:
            chains[chain.id] = ("".join(codes), relative_accessibility(codes, np.array(values)))
    return chains


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--designs", type=Path, default=DEFAULT_DESIGNS)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--limit-parents", type=int, default=0)
    args = parser.parse_args()

    designs = pd.read_parquet(args.designs)
    if args.limit_parents:
        designs = designs.head(args.limit_parents)
    profile = pd.read_parquet(args.profile)
    profile = profile[profile["parent_sequence_id"].isin(set(designs["parent"]))]

    per_parent_rsa: dict[str, np.ndarray] = {}
    unmatched: list[str] = []
    missing_structure: list[str] = []
    structure_cache: dict[str, dict[str, tuple[str, np.ndarray]]] = {}

    for position, design in enumerate(designs.itertuples(index=False)):
        code = str(design.pdb).lower()
        parent = str(design.parent)
        if position % 200 == 0:
            print(f"parents {position:,}/{len(designs):,}", flush=True)
        if code not in structure_cache:
            path = fetch_structure(code, args.cache_dir)
            structure_cache[code] = chain_accessibility(path) if path is not None else {}
        chains = structure_cache[code]
        if not chains:
            missing_structure.append(parent)
            continue

        sequence = str(design.seq)
        best: tuple[int, np.ndarray] | None = None
        for chain_sequence, chain_rsa in chains.values():
            offset = align_chain_to_design(sequence, chain_sequence)
            if offset is None:
                continue
            if best is None or len(chain_rsa) > len(best[1]):
                best = (offset, chain_rsa)
        if best is None:
            unmatched.append(parent)
            continue
        per_parent_rsa[parent] = project_onto_design(len(sequence), best[0], best[1])

    covered = profile["parent_sequence_id"].isin(per_parent_rsa)
    scored = profile.loc[covered].copy()
    mean_rsa = np.full(len(scored), np.nan)
    buried_fraction = np.full(len(scored), np.nan)
    observed = np.zeros(len(scored), dtype=np.int64)
    for parent, block in scored.groupby("parent_sequence_id", sort=False):
        rows = scored.index.get_indexer(block.index)
        window_mean, window_buried, window_observed = window_accessibility(
            per_parent_rsa[str(parent)],
            block["start"].to_numpy(),
            block["end"].to_numpy(),
        )
        mean_rsa[rows] = window_mean
        buried_fraction[rows] = window_buried
        observed[rows] = window_observed
    scored["window_mean_rsa"] = mean_rsa
    scored["window_buried_fraction"] = buried_fraction
    scored["window_observed_residues"] = observed
    resolved = scored.dropna(subset=["window_mean_rsa"]).copy()
    resolved["window_is_buried"] = resolved["window_mean_rsa"] < BURIED_RSA_THRESHOLD

    by_class = (
        resolved.groupby("binder_class_student_oof")
        .agg(
            n=("window_mean_rsa", "size"),
            mean_rsa=("window_mean_rsa", "mean"),
            median_rsa=("window_mean_rsa", "median"),
            buried_share=("window_is_buried", "mean"),
        )
        .reindex(["strong", "weak", "nonbinder"])
        .dropna(how="all")
    )
    flagged = resolved[resolved["binder_class_student_oof"].isin(["strong", "weak"])]

    report = {
        "profile": str(args.profile.relative_to(REPO_ROOT)),
        "coverage": {
            "parents_requested": int(len(designs)),
            "parents_with_structure_match": len(per_parent_rsa),
            "parents_missing_structure": len(missing_structure),
            "parents_sequence_unmatched": len(unmatched),
            "occurrences_scored": int(len(resolved)),
            "occurrences_in_covered_parents": int(len(scored)),
        },
        "buried_definition": f"window mean relative SASA < {BURIED_RSA_THRESHOLD}",
        "by_binder_class": {
            str(index): {
                "n": int(row["n"]),
                "mean_rsa": float(row["mean_rsa"]),
                "median_rsa": float(row["median_rsa"]),
                "buried_share": float(row["buried_share"]),
            }
            for index, row in by_class.iterrows()
        },
        "flagged_epitopes": {
            "n": int(len(flagged)),
            "buried_share": float(flagged["window_is_buried"].mean()) if len(flagged) else None,
            "surface_share": (
                float(1.0 - flagged["window_is_buried"].mean()) if len(flagged) else None
            ),
        },
        "risk_vs_exposure_spearman": float(
            resolved["overall_mhci_risk"].corr(
                resolved["window_mean_rsa"], method="spearman"
            )
        ),
        "boundaries": [
            "Burial does not gate MHC-I presentation; peptides derive from unfolded "
            "protein. This measures design actionability, not immunogenicity.",
            "SASA is computed on the full deposited model, so interface residues in "
            "multi-chain entries read as buried.",
            "Windows with fewer than five observed residues are excluded rather than "
            "averaged over partial coverage.",
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    resolved[
        [
            "parent_sequence_id",
            "start",
            "end",
            "peptide",
            "binder_class_student_oof",
            "overall_mhci_risk",
            "window_mean_rsa",
            "window_buried_fraction",
            "window_is_buried",
        ]
    ].to_parquet(args.output_dir / "epitope_accessibility.parquet", index=False)

    coverage = report["coverage"]
    lines = [
        "# Are flagged MHC-I epitopes reachable by redesign?",
        "",
        f"{coverage['occurrences_scored']:,} 9-mer windows across "
        f"{coverage['parents_with_structure_match']:,} PDA parents with a matched "
        "deposited structure.",
        f"Buried is defined as {report['buried_definition']}.",
        "",
        "| Student binder class | Windows | Mean RSA | Buried share |",
        "| --- | ---: | ---: | ---: |",
    ]
    for name, block in report["by_binder_class"].items():
        lines.append(
            f"| {name} | {block['n']:,} | {block['mean_rsa']:.3f} "
            f"| {block['buried_share']:.1%} |"
        )
    flagged_block = report["flagged_epitopes"]
    if flagged_block["buried_share"] is not None:
        lines += [
            "",
            f"Of {flagged_block['n']:,} windows flagged strong or weak, "
            f"{flagged_block['surface_share']:.1%} are surface-exposed and can be "
            f"resampled directly; {flagged_block['buried_share']:.1%} are buried and "
            "need a backbone-aware fix.",
        ]
    lines += [
        "",
        f"Spearman between MHC-I risk and window exposure: "
        f"{report['risk_vs_exposure_spearman']:.3f}.",
        "",
        "## Boundaries",
        "",
    ]
    lines += [f"- {boundary}" for boundary in report["boundaries"]]
    lines.append("")
    (args.output_dir / "REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
