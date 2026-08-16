#!/usr/bin/env python3
"""Run the frozen team MHC-I surrogate on one traceable PDA design."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.dataset as ds

from re_agent.agent.tools import screen_candidate

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PARENT = "pda:9s14:0"
DEFAULT_POOL = ROOT / "data" / "processed" / "pda_designs.parquet"
DEFAULT_TRAINING_DATA = ROOT / "data" / "dataset.parquet 2"
DEFAULT_CHECKPOINT = ROOT / "models" / "chao1" / "cv5_heads.pkl 2"


def _sample_metadata(row: pd.Series, training_data: Path) -> dict[str, Any]:
    sequence = str(row["seq"])
    parent_hash = hashlib.sha256(sequence.encode()).hexdigest()[:16]
    sample_windows = {sequence[index : index + 9] for index in range(len(sequence) - 8)}
    overlap = {"parent_rows": None, "overlapping_9mers": None}
    if training_data.exists():
        dataset = ds.dataset(training_data, format="parquet")
        parent_rows = dataset.to_table(
            columns=["parent_sequence_hash"],
            filter=ds.field("parent_sequence_hash") == parent_hash,
        ).num_rows
        peptide_rows = dataset.to_table(
            columns=["peptide"],
            filter=ds.field("peptide").isin(sorted(sample_windows)),
        ).to_pylist()
        overlap = {
            "parent_rows": parent_rows,
            "overlapping_9mers": len({item["peptide"] for item in peptide_rows}),
            "sample_9mers": len(sample_windows),
        }
    return {
        "source": "Protein Design Archive",
        "snapshot": "20260802_data_curated.json",
        "parent": row["parent"],
        "pdb": row["pdb"],
        "length": int(row["length"]),
        "tm_score_to_nearest_natural_fold": float(row["tm_natural"]),
        "nearest_natural_structure": row["tm_partner"],
        "novelty_bin": row["novelty_bin"],
        "release_date": row["release_date"],
        "in_rcsb_training_pool": bool(row["in_rcsb_pool"]),
        "training_overlap_audit": overlap,
        "smoke_test_only": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", default=DEFAULT_PARENT)
    parser.add_argument("--pool", type=Path, default=DEFAULT_POOL)
    parser.add_argument("--training-data", type=Path, default=DEFAULT_TRAINING_DATA)
    parser.add_argument(
        "--checkpoint",
        type=Path,
        action="append",
        dest="checkpoints",
        help="checkpoint to attach; repeat to compare models",
    )
    args = parser.parse_args()

    designs = pd.read_parquet(args.pool)
    matches = designs.loc[designs["parent"] == args.parent]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one PDA row for {args.parent!r}, found {len(matches)}")
    row = matches.iloc[0]
    metadata = _sample_metadata(row, args.training_data)
    checkpoints = [
        (path if path.is_absolute() else ROOT / path).resolve()
        for path in (args.checkpoints or [DEFAULT_CHECKPOINT])
    ]
    artifact = screen_candidate.invoke(
        {
            "sequence": row["seq"],
            "candidate_id": row["parent"],
            "mhci_surrogate_checkpoints": [
                str(checkpoint.relative_to(ROOT)) for checkpoint in checkpoints
            ],
            "source_metadata": metadata,
        }
    )
    assessment = artifact["payload"]["assessment"]
    surrogates = assessment["mhc_i_surrogate_results"]
    prediction_sets = [row["predictions"] for row in surrogates]
    output = {
        "artifact_path": artifact["path"],
        "artifact_sha256": artifact["sha256"],
        "candidate_id": assessment["candidate_id"],
        "sequence_length": len(assessment["sequence"]),
        "training_overlap_audit": metadata["training_overlap_audit"],
        "mhci_surrogates": [
            {
                "adapter_id": surrogate["adapter_id"],
                "status": surrogate["status"],
                "checkpoint_sha256": surrogate["provenance"]["parameters"][
                    "checkpoint_sha256"
                ],
                "summary": surrogate["protein_summary"],
            }
            for surrogate in surrogates
        ],
        "custom_predictions_identical": (
            all(predictions == prediction_sets[0] for predictions in prediction_sets[1:])
            if len(prediction_sets) > 1
            else None
        ),
        "mhc_ii_provider_statuses": {
            result["provider_id"]: result["status"] for result in assessment["mhc_results"]
        },
        "combined_rank_score": assessment["combined_rank_score"],
        "claim_boundary": artifact["payload"]["claim_boundary"],
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
