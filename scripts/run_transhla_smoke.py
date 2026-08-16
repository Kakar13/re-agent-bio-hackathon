#!/usr/bin/env python3
"""Run pretrained ESM2-derived TransHLA-I/II on a full protein sequence."""

from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
from pathlib import Path

os.environ.setdefault("HF_HOME", "/tmp/huggingface")
os.environ.setdefault("HF_MODULES_CACHE", "/tmp/huggingface/modules")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import esm
import torch
from esm.data import Alphabet
from esm.model.esm2 import ESM2
from safetensors.torch import load_model
from transformers import AutoTokenizer

TOP7 = (
    "DIQVQVNIDDNGKNFDYTYTVTTESELQKVLNELKDYIKKQGAKRVRISITARTKKEAEKFAAILI"
    "KVFAELGYNDINVTFDGDTVTVEGQLE"
)


def _empty_esm2_650m() -> tuple[ESM2, Alphabet]:
    """Build the ESM2 architecture without downloading duplicate base weights."""
    alphabet = Alphabet.from_architecture("ESM-1b")
    model = ESM2(
        num_layers=33,
        embed_dim=1280,
        attention_heads=20,
        alphabet=alphabet,
        token_dropout=True,
    )
    return model, alphabet


def windows(sequence: str, lengths: range) -> list[tuple[int, str]]:
    return [
        (start + 1, sequence[start : start + length])
        for length in lengths
        for start in range(len(sequence) - length + 1)
    ]


def score_arm(
    sequence: str,
    *,
    model_path: Path,
    lengths: range,
    padded_length: int,
    device: torch.device,
) -> list[dict[str, object]]:
    peptides = windows(sequence, lengths)
    tokenizer = AutoTokenizer.from_pretrained("/tmp/esm2-tokenizer", local_files_only=True)
    esm.pretrained.esm2_t33_650M_UR50D = _empty_esm2_650m
    sys.path.insert(0, str(model_path))
    try:
        if model_path.name == "transhla-i":
            from configuration_TransHLA_I import TransHLA_I_Config
            from modeling_TransHLA_I import TransHLA_I_Model

            model = TransHLA_I_Model(TransHLA_I_Config())
        else:
            from configuration_TransHLA_II import TransHLA_II_Config
            from modeling_TransHLA_II import TransHLA_II_Model

            model = TransHLA_II_Model(TransHLA_II_Config())
    finally:
        sys.path.pop(0)
    load_model(model, model_path / "model.safetensors", strict=False)
    model = model.to(device).eval()
    encoded = tokenizer([peptide for _, peptide in peptides])["input_ids"]
    encoded = [row + ([1] * (padded_length - len(row))) for row in encoded]
    tensor = torch.tensor(encoded, dtype=torch.long, device=device)

    scored: list[dict[str, object]] = []
    with torch.inference_mode():
        for offset in range(0, len(tensor), 64):
            outputs, _representations = model(tensor[offset : offset + 64])
            probabilities = outputs[:, 1].detach().cpu().tolist()
            for (start, peptide), probability in zip(
                peptides[offset : offset + 64], probabilities, strict=True
            ):
                scored.append(
                    {
                        "start": start,
                        "end": start + len(peptide) - 1,
                        "peptide": peptide,
                        "probability": round(float(probability), 6),
                        "predicted_epitope": bool(probability >= 0.5),
                    }
                )
    del model
    if device.type == "mps":
        torch.mps.empty_cache()
    return sorted(scored, key=lambda row: float(row["probability"]), reverse=True)


def summarize(scored: list[dict[str, object]]) -> dict[str, object]:
    probabilities = [float(row["probability"]) for row in scored]
    positive_count = sum(bool(row["predicted_epitope"]) for row in scored)
    return {
        "window_count": len(scored),
        "positive_count": positive_count,
        "positive_fraction": round(positive_count / max(len(scored), 1), 6),
        "mean_probability": round(statistics.fmean(probabilities), 6),
        "median_probability": round(statistics.median(probabilities), 6),
        "top_hits": scored[:20],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", default=TOP7)
    parser.add_argument("--sequence-id", default="Top7_1QYS")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/immuno_risk/transhla_top7_smoke.json"),
    )
    args = parser.parse_args()

    sequence = "".join(args.sequence.split()).upper()
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    mhc_i = score_arm(
        sequence,
        model_path=Path("/tmp/transhla-i"),
        lengths=range(8, 15),
        padded_length=16,
        device=device,
    )
    mhc_ii = score_arm(
        sequence,
        model_path=Path("/tmp/transhla-ii"),
        lengths=range(13, 22),
        padded_length=23,
        device=device,
    )
    result = {
        "sequence_id": args.sequence_id,
        "sequence": sequence,
        "device": str(device),
        "interpretation": (
            "Out-of-distribution smoke test on a de novo protein. "
            "TransHLA probabilities are not validation or calibrated de novo risk."
        ),
        "mhc_i": summarize(mhc_i),
        "mhc_ii": summarize(mhc_ii),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
