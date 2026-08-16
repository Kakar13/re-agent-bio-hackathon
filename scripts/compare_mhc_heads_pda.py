#!/usr/bin/env python3
"""Compare the legacy chao1 MHCflurry head against the NetMHCpan student on PDA.

The external NY-ESO cohort is pre-screened on NetMHCpan and therefore cannot
separate binding models. The full PDA profile can: it spans the whole rank range
with 96% non-binders, so agreement with the NetMHCpan teacher is measurable.

Both scores are evaluated on the same rows against the same teacher. The student
is read from its out-of-fold column, so neither model has seen these labels for
the rows it is scored on.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.e2e_pls.benchmark import teacher_agreement
from re_agent.e2e_pls.netmhcpan_student import (
    STRONG_BINDER_RANK,
    WEAK_BINDER_RANK,
    rank_to_propensity,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROFILE = (
    REPO_ROOT / "data/processed/profiles/a0201-pda-mhci-profile-v4/pda_mhci_profile.parquet"
)
DEFAULT_OUTPUT = REPO_ROOT / "results/benchmarks/pda_mhc_head_comparison"
SCORES = {
    "chao1_mhcflurry_presentation": "chao1_mhcflurry_presentation",
    "netmhcpan_student_el_oof": "student_el_propensity_oof",
    "netmhcpan_student_ba_oof": "student_ba_propensity_oof",
    "netmhcpan_student_affinity_oof": "student_ba_affinity_score_oof",
}


def score_against_teacher(frame: pd.DataFrame, column: str) -> dict:
    teacher_rank = frame["netmhcpan_el_rank"].to_numpy()
    return teacher_agreement(
        frame[column].to_numpy(),
        teacher_rank,
        rank_to_propensity(teacher_rank),
        strong_rank=STRONG_BINDER_RANK,
        weak_rank=WEAK_BINDER_RANK,
    )


def _affinity_recovery(frame: pd.DataFrame) -> dict | None:
    """How closely the student recovers the teacher's nanomolar IC50.

    chao1 has no comparable lane: its MHC head emits a calibrated cosine
    similarity, not an affinity, so only the student can be scored here.
    """

    column = "student_predicted_ba_ic50_nm_oof"
    if column not in frame or "netmhcpan_ba_ic50_nm" not in frame:
        return None
    predicted = frame[column].to_numpy(dtype=float)
    actual = frame["netmhcpan_ba_ic50_nm"].to_numpy(dtype=float)
    fold_error = np.exp(np.abs(np.log(predicted) - np.log(actual)))
    binders = actual <= 500.0
    return {
        "n_rows": int(len(frame)),
        "n_binders_ic50_500nm": int(binders.sum()),
        "median_fold_error": float(np.median(fold_error)),
        "median_fold_error_binders": float(np.median(fold_error[binders])),
        "within_2x": float((fold_error <= 2.0).mean()),
        "within_2x_binders": float((fold_error[binders] <= 2.0).mean()),
        "within_5x_binders": float((fold_error[binders] <= 5.0).mean()),
        "spearman_vs_teacher_ic50": float(
            pd.Series(predicted).corr(pd.Series(actual), method="spearman")
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    frame = pd.read_parquet(args.profile)
    # Peptides repeat across parents; deduplicate so shared 9-mers are not
    # counted once per occurrence and allowed to dominate the correlation.
    unique = frame.drop_duplicates(subset="peptide", keep="first").reset_index(drop=True)

    report = {
        "profile": str(args.profile.relative_to(REPO_ROOT)),
        "n_occurrences": int(len(frame)),
        "n_unique_peptides": int(len(unique)),
        "teacher": "NetMHCpan 4.1 EL, HLA-A*02:01",
        "prevalence": {
            "strong_binders": int((unique["netmhcpan_el_rank"] <= STRONG_BINDER_RANK).sum()),
            "weak_or_better": int((unique["netmhcpan_el_rank"] <= WEAK_BINDER_RANK).sum()),
        },
        "agreement": {
            name: score_against_teacher(unique, column)
            for name, column in SCORES.items()
            if column in unique
        },
        "affinity_recovery": _affinity_recovery(unique),
        "boundaries": [
            "The chao1 MHC head was calibrated against MHCflurry, so this measures "
            "cross-teacher disagreement, not an error rate against experimental data.",
            "Student values are out-of-fold; chao1 never saw these labels at all, so "
            "neither model is scored in sample.",
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# chao1 MHCflurry head vs NetMHCpan student on the PDA corpus",
        "",
        f"{report['n_unique_peptides']:,} unique 9-mers from "
        f"{report['n_occurrences']:,} occurrences, scored against NetMHCpan 4.1 EL.",
        f"{report['prevalence']['strong_binders']:,} strong binders "
        f"(EL rank <= {STRONG_BINDER_RANK}), "
        f"{report['prevalence']['weak_or_better']:,} at rank <= {WEAK_BINDER_RANK}.",
        "",
        "| Score | Spearman vs teacher | Strong AUROC | Strong AP | Weak AUROC | Weak AP |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for name, metrics in report["agreement"].items():
        lines.append(
            f"| {name} | {metrics['spearman_vs_teacher_el_rank']:.3f} "
            f"| {metrics['strong_binder_roc_auc']:.3f} "
            f"| {metrics['strong_binder_average_precision']:.3f} "
            f"| {metrics['weak_binder_roc_auc']:.3f} "
            f"| {metrics['weak_binder_average_precision']:.3f} |"
        )
    recovery = report.get("affinity_recovery")
    if recovery:
        lines += [
            "",
            "## Recovering the teacher's nanomolar affinity",
            "",
            "Only the student has an affinity lane; the chao1 MHC head emits a "
            "calibrated cosine similarity with no nanomolar interpretation.",
            "",
            f"- Spearman against teacher IC50: "
            f"{recovery['spearman_vs_teacher_ic50']:.3f}",
            f"- Median fold error: {recovery['median_fold_error']:.2f}x overall, "
            f"{recovery['median_fold_error_binders']:.2f}x on the "
            f"{recovery['n_binders_ic50_500nm']:,} binders at IC50 <= 500 nM",
            f"- Within 2x of the teacher: {recovery['within_2x']:.1%} overall, "
            f"{recovery['within_2x_binders']:.1%} on binders",
            f"- Within 5x on binders: {recovery['within_5x_binders']:.1%}",
        ]
    lines += ["", "## Boundaries", ""]
    lines += [f"- {boundary}" for boundary in report["boundaries"]]
    lines.append("")
    (args.output_dir / "REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
