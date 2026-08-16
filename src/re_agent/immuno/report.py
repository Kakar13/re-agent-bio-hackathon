"""`predict` entry point: binder FASTA in, risk + confidence + heatmap out."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from re_agent.immuno.config import PATHS, ensure_dirs
from re_agent.immuno.data import clean_sequence, read_fasta
from re_agent.immuno.explain import plot_heatmap, top_hotspots
from re_agent.immuno.score import (
    Bundle,
    build_reference_risks,
    load_bundle,
    score_sequence,
    score_to_dict,
)


def ensure_reference(bundle: Bundle, model_name: str) -> None:
    if len(bundle.reference_risks) == 0:
        print("building natural reference cohort risks (one-off)...")
        build_reference_risks(bundle, model_name)


def assess(
    sequences: dict[str, str], model_name: str = "mean_teacher", mc_passes: int = 20
) -> list[dict]:
    ensure_dirs()
    bundle = load_bundle(model_name)
    ensure_reference(bundle, model_name)

    reports = []
    for name, seq in sequences.items():
        score = score_sequence(seq, bundle, name=name, mc_passes=mc_passes)
        payload = score_to_dict(score)
        payload["hotspots"] = top_hotspots(score.per_residue, score.sequence)
        payload["model"] = model_name

        fig_path = PATHS.figures / f"{name}_heatmap.png"
        plot_heatmap(score, fig_path)
        payload["heatmap"] = str(fig_path)

        json_path = PATHS.reports / f"{name}.json"
        json_path.write_text(json.dumps(payload, indent=2))
        payload["report"] = str(json_path)
        reports.append(payload)

        print(
            f"\n{name}: risk {payload['risk']:.3f} "
            f"({payload['risk_percentile_vs_natural']:.0f}th pct vs length-matched natural, "
            f"peak window {payload['peak_window_risk']:.3f}) | "
            f"confidence {payload['confidence']:.2f} "
            f"(stability {payload['confidence_breakdown']['stability']:.2f}, "
            f"agreement {payload['confidence_breakdown']['agreement']:.2f}, "
            f"familiarity {payload['confidence_breakdown']['familiarity']:.2f})"
        )
        print(f"  top regions ({len(payload['top_regions'])} shown):")
        for region in payload["top_regions"][:5]:
            print(
                f"    {region['start']:>4}-{region['end']:<4} {region['peptide']}  "
                f"risk {region['risk']:.3f}  conf {region['confidence']:.2f}"
            )
        print(f"  heatmap -> {fig_path}")
    return reports


def main() -> None:
    parser = argparse.ArgumentParser(description="Assess immunogenicity risk of a protein/binder")
    parser.add_argument("input", help="FASTA path, or a raw amino-acid sequence")
    parser.add_argument("--model", default="mean_teacher", choices=["mean_teacher", "baseline"])
    parser.add_argument("--name", default=None, help="name when passing a raw sequence")
    parser.add_argument("--mc-passes", type=int, default=20)
    args = parser.parse_args()

    path = Path(args.input)
    if path.exists():
        sequences = read_fasta(path)
    else:
        seq = clean_sequence(args.input)
        if not seq:
            raise SystemExit("input is neither an existing FASTA path nor a valid sequence")
        sequences = {args.name or "query": seq}

    if not sequences:
        raise SystemExit("no usable sequences found")
    assess(sequences, model_name=args.model, mc_passes=args.mc_passes)


if __name__ == "__main__":
    main()
