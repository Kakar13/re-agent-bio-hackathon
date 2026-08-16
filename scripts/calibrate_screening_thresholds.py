#!/usr/bin/env python3
"""Pick student cutoffs by recall instead of by inverted rank.

The deployed binder classes threshold a predicted percentile rank, which makes
them inherit the student's conservative bias: it calls fewer strong binders than
the teacher does. For a screen that is the wrong direction of error, so this
reports what recall the default buys and what each alternative cutoff costs.

Everything is measured on out-of-fold predictions, so the operating points are
not fitted on rows the scoring head trained on.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from re_agent.e2e_pls.netmhcpan_student import (
    STRONG_BINDER_RANK,
    WEAK_BINDER_RANK,
    screening_thresholds,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROFILE = (
    REPO_ROOT / "data/processed/profiles/a0201-pda-mhci-profile-v4/pda_mhci_profile.parquet"
)
DEFAULT_OUTPUT = REPO_ROOT / "results/benchmarks/screening_calibration"


def render(report: dict) -> str:
    lines = [
        "# Recall-driven screening cutoffs for the NetMHCpan student",
        "",
        f"Measured on {report['n_unique_peptides']:,} out-of-fold PDA peptides.",
        "",
    ]
    for name, block in report["channels"].items():
        default = block["default_cutoff"]
        lines += [
            f"## {name} (teacher EL rank <= {block['teacher_binder_rank']})",
            "",
            f"{block['n_teacher_binders']:,} teacher binders. The default "
            f"`{default['rule']}` flags {default['n_flagged']:,} peptides at "
            f"{default['recall']:.1%} recall and {default['precision']:.1%} precision.",
            "",
            "| Target recall | Propensity cutoff | Predicted rank | Recall | Precision | Flagged |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        for point in block["operating_points"]:
            lines.append(
                f"| {point['target_recall']:.0%} "
                f"| {point['propensity_threshold']:.4f} "
                f"| {point['implied_predicted_rank']:.3f} "
                f"| {point['recall']:.1%} "
                f"| {point['precision']:.1%} "
                f"| {point['n_flagged']:,} |"
            )
        lines.append("")
    lines += [
        "## How to read this",
        "",
        "Precision falls as recall rises because the extra flagged peptides are",
        "mostly teacher non-binders. The right operating point depends on what a",
        "false flag costs: during generation it costs one resampled window, which",
        "is cheap, so high recall is usually correct.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    frame = pd.read_parquet(
        args.profile,
        columns=["peptide", "student_el_propensity_oof", "netmhcpan_el_rank"],
    ).drop_duplicates(subset="peptide")

    propensities = frame["student_el_propensity_oof"].to_numpy()
    teacher_ranks = frame["netmhcpan_el_rank"].to_numpy()
    report = {
        "profile": str(args.profile.relative_to(REPO_ROOT)),
        "n_unique_peptides": int(len(frame)),
        "scoring": "out-of-fold student EL propensity",
        "channels": {
            "strong_binder": screening_thresholds(
                propensities, teacher_ranks, binder_rank=STRONG_BINDER_RANK
            ),
            "weak_or_better": screening_thresholds(
                propensities, teacher_ranks, binder_rank=WEAK_BINDER_RANK
            ),
        },
        "boundaries": [
            "Recall is measured against the NetMHCpan teacher, not experimental binding.",
            "Cutoffs are fitted on PDA de novo peptides and may shift on other domains.",
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    markdown = render(report)
    (args.output_dir / "REPORT.md").write_text(markdown)
    print(markdown)


if __name__ == "__main__":
    main()
