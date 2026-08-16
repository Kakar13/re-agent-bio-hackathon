#!/usr/bin/env python3
"""Compare the chao1 MHCflurry binding head against NetMHCpan BA/EL on an
external de novo HLA-A*02:01 cohort (Visani et al. NY-ESO-1 designs).

Two endpoints are reported separately because only the first is unbiased:

1. Teacher imitation. The student was distilled from NetMHCpan on PDA de novo
   proteins; these peptides come from an unrelated design campaign, so agreement
   with freshly called NetMHCpan measures out-of-corpus generalization.
2. Experimental T-cell activation. The published cohort was pre-screened with
   NetMHCpan, so presentation scores are range-restricted here and this endpoint
   can only bound discrimination, never rank the two binding models fairly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.e2e_pls.benchmark import (
    binder_agreement,
    evaluate_external_scores,
    imitation_metrics,
)
from re_agent.e2e_pls.netmhcpan_student import (
    load_student_ensemble,
    predict_student_ensemble,
    propensity_to_rank,
    rank_to_propensity,
)
from re_agent.immuno.e2e_pls_pickle import (
    _ESM2SegmentEncoder,
    _interp,
    load_validated_bundle,
)
from re_agent.immuno.netmhcpan import NetMHCpanTeacher, NetMHCpanTeacherConfig

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COHORT = REPO_ROOT / "data/processed/benchmarks/visani_nyeso_hla_a0201.csv"
DEFAULT_CHAO1 = REPO_ROOT / "models/chao1/cv5_heads.pkl 2"
DEFAULT_STUDENT = REPO_ROOT / "models/a0201-netmhcpan-pda-cv5-v4/checkpoint"
DEFAULT_CACHE = REPO_ROOT / "data/raw/netmhcpan-nyeso"
DEFAULT_OUTPUT = REPO_ROOT / "results/benchmarks/nyeso_a0201"
ALLELE = "HLA-A*02:01"
PEPTIDE_LENGTH = 9


def load_cohort(path: Path) -> pd.DataFrame:
    """Deduplicate the published cohort to one row per designed peptide."""

    frame = pd.read_csv(path)
    required = {"sequence", "resp", "r_mean"}
    missing = required - set(frame)
    if missing:
        raise ValueError(f"cohort is missing columns: {sorted(missing)}")

    # The published table carries unstimulated assay controls with no peptide.
    controls = frame["sequence"].isna()
    if controls.any():
        print(f"dropping {int(controls.sum())} assay control rows without a peptide")
        frame = frame.loc[~controls].copy()

    frame["peptide"] = frame["sequence"].astype(str).str.strip().str.upper()
    bad_length = frame["peptide"].str.len().ne(PEPTIDE_LENGTH)
    if bad_length.any():
        raise ValueError(f"cohort contains {int(bad_length.sum())} non-9-mer peptides")

    # One peptide appears twice under different design protocols. Averaging the
    # replicate activation keeps every peptide independent for the bootstrap.
    collapsed = (
        frame.groupby("peptide", as_index=False)
        .agg(
            r_mean=("r_mean", "mean"),
            resp=("resp", "max"),
            netmhcpan_published_binder_el_0p5=(
                "is_binder_by_netmhcpan_at_el_rank_0p5",
                "max",
            ),
            netmhcpan_published_binder_el_2p0=(
                "is_binder_by_netmhcpan_at_el_rank_2p0",
                "max",
            ),
            n_design_rows=("peptide", "size"),
        )
        .sort_values("peptide", ignore_index=True)
    )
    collapsed["resp"] = collapsed["resp"].astype(bool)
    return collapsed


def embed_peptides(peptides: list[str], *, batch_size: int) -> np.ndarray:
    """Mean-pool frozen ESM-2 over the nine peptide residues.

    Synthetic peptides are pulsed onto cells without flanking context, so both
    models receive the bare 9-mer. This matches the student's training pooling
    and is identical input for chao1, keeping the head-to-head fair.
    """

    encoder = _ESM2SegmentEncoder(batch_size=batch_size)
    embeddings = encoder.embed(peptides)
    pooled = np.stack([embedding[:PEPTIDE_LENGTH].mean(axis=0) for embedding in embeddings])
    return pooled.astype("float32")


def score_chao1(bundle: dict, pooled: np.ndarray) -> np.ndarray:
    """Reproduce the chao1 MHCflurry-derived presentation head exactly."""

    mhc = bundle["mhc"]
    projected = pooled @ mhc["projection_weight"].T + mhc["projection_bias"]
    projected /= np.clip(np.linalg.norm(projected, axis=1, keepdims=True), 1e-6, None)
    cosine = projected @ mhc["centroids"][ALLELE]
    return _interp(mhc["calibrators"][ALLELE], cosine)


def _auc_row(name: str, metrics: dict) -> str:
    low, high = metrics["bootstrap_95_ci"]["roc_auc"]
    spearman_low, spearman_high = metrics["bootstrap_95_ci"]["spearman_vs_activity"]
    return (
        f"| {name} | {metrics['roc_auc']:.3f} [{low:.3f}, {high:.3f}] "
        f"| {metrics['average_precision']:.3f} "
        f"| {metrics['spearman_vs_activity']:.3f} [{spearman_low:.3f}, {spearman_high:.3f}] |"
    )


def write_report(report: dict, path: Path) -> None:
    """Render the two endpoints as a reviewable markdown summary."""

    cohort = report["cohort"]
    imitation = report["endpoint_1_teacher_imitation"]
    activation = report["endpoint_2_experimental_activation"]
    restriction = report["range_restriction"]

    lines = [
        "# External benchmark: de novo HLA-A*02:01 NY-ESO-1 designs",
        "",
        f"Cohort: {cohort['n_peptides']} unique 9-mers, "
        f"{cohort['n_responders']} experimental responders, allele {cohort['allele']}.",
        f"Origin: {cohort['design_origin']}.",
        "",
        "## Endpoint 1 - out-of-corpus teacher imitation (unbiased)",
        "",
        "The student was distilled from NetMHCpan on PDA de novo proteins. These peptides",
        "come from an unrelated design campaign, so agreement with freshly called",
        "NetMHCpan measures generalization beyond the training corpus.",
        "",
        "| Channel | Spearman vs teacher | MAE | RMSE |",
        "| --- | --- | --- | --- |",
    ]
    for channel in ("el", "ba"):
        entry = imitation[channel]
        lines.append(
            f"| {channel.upper()} | {entry['spearman_vs_teacher_propensity']:.3f} "
            f"| {entry['propensity_mae']:.4f} | {entry['propensity_rmse']:.4f} |"
        )

    strong = imitation["binder_agreement"]["el_strong"]
    weak = imitation["binder_agreement"]["el_weak"]
    lines += [
        "",
        f"Binder-call agreement: {weak['agreement']:.1%} at EL rank <= 2.0 "
        f"({weak['student_binders']}/{weak['teacher_binders']} student/teacher calls), "
        f"{strong['agreement']:.1%} at EL rank <= 0.5 "
        f"({strong['student_binders']}/{strong['teacher_binders']}).",
        "",
        "## Endpoint 2 - experimental T-cell activation (limits result)",
        "",
        f"{restriction['teacher_el_rank_le_2p0']}/{restriction['n_peptides']} peptides are "
        f"already NetMHCpan EL rank <= 2.0 and "
        f"{restriction['teacher_el_rank_le_0p5']}/{restriction['n_peptides']} are <= 0.5.",
        restriction["note"],
        "",
        f"n = {activation['n_rows']} peptides, {activation['n_positive']} responders, "
        f"{activation['n_negative']} non-responders.",
        "",
        "| Score | ROC AUC [95% CI] | Avg precision | Spearman vs activation [95% CI] |",
        "| --- | --- | --- | --- |",
    ]
    lines += [_auc_row(name, metrics) for name, metrics in activation["metrics"].items()]
    lines += [
        "",
        "All confidence intervals overlap, so this cohort does not distinguish the chao1",
        "MHCflurry head from the NetMHCpan student or teacher on activation. The real",
        "NetMHCpan teacher itself reaches only "
        f"{activation['metrics']['netmhcpan_teacher_el']['roc_auc']:.3f} AUC, which bounds",
        "what any faithful surrogate of it can achieve here.",
        "",
        "## Boundaries",
        "",
    ]
    lines += [f"- {boundary}" for boundary in report["boundaries"]]
    lines.append("")
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--chao1-checkpoint", type=Path, default=DEFAULT_CHAO1)
    parser.add_argument("--student-checkpoint", type=Path, default=DEFAULT_STUDENT)
    parser.add_argument("--teacher-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--embedding-batch-size", type=int, default=32)
    parser.add_argument("--bootstrap-samples", type=int, default=2_000)
    args = parser.parse_args()

    cohort = load_cohort(args.cohort)
    peptides = cohort["peptide"].tolist()
    print(f"cohort: {len(peptides)} unique 9-mers, {int(cohort['resp'].sum())} responders")

    print("calling the NetMHCpan teacher through IEDB ...")
    teacher = NetMHCpanTeacher(
        cache_dir=args.teacher_cache,
        config=NetMHCpanTeacherConfig(allele=ALLELE, peptide_length=PEPTIDE_LENGTH),
    )
    labeled = teacher.label(cohort, peptide_column="peptide")
    if labeled["netmhcpan_el_rank"].isna().any() or labeled["netmhcpan_ba_rank"].isna().any():
        raise RuntimeError("NetMHCpan teacher returned incomplete labels")

    print("embedding peptides with frozen ESM-2 ...")
    pooled = embed_peptides(peptides, batch_size=args.embedding_batch_size)

    print("scoring chao1 and the NetMHCpan student ...")
    bundle, chao1_sha256 = load_validated_bundle(args.chao1_checkpoint)
    students, student_sha256, student_manifest = load_student_ensemble(args.student_checkpoint)
    student_predictions = predict_student_ensemble(students, pooled)

    scores = labeled.copy()
    scores["chao1_mhcflurry_presentation"] = score_chao1(bundle, pooled)
    scores["student_el_propensity"] = student_predictions[:, 0].astype(np.float64)
    scores["student_ba_propensity"] = student_predictions[:, 1].astype(np.float64)
    scores["student_predicted_el_rank"] = propensity_to_rank(scores["student_el_propensity"])
    scores["student_predicted_ba_rank"] = propensity_to_rank(scores["student_ba_propensity"])
    scores["teacher_el_propensity"] = rank_to_propensity(scores["netmhcpan_el_rank"].to_numpy())
    scores["teacher_ba_propensity"] = rank_to_propensity(scores["netmhcpan_ba_rank"].to_numpy())

    imitation = {
        "el": imitation_metrics(
            scores["student_el_propensity"].to_numpy(),
            scores["teacher_el_propensity"].to_numpy(),
        ),
        "ba": imitation_metrics(
            scores["student_ba_propensity"].to_numpy(),
            scores["teacher_ba_propensity"].to_numpy(),
        ),
        "binder_agreement": {
            "el_strong": binder_agreement(
                scores["student_predicted_el_rank"].to_numpy(),
                scores["netmhcpan_el_rank"].to_numpy(),
                0.5,
            ),
            "el_weak": binder_agreement(
                scores["student_predicted_el_rank"].to_numpy(),
                scores["netmhcpan_el_rank"].to_numpy(),
                2.0,
            ),
        },
    }

    activation = evaluate_external_scores(
        scores,
        {
            "chao1_mhcflurry_presentation": "chao1_mhcflurry_presentation",
            "netmhcpan_student_el": "student_el_propensity",
            "netmhcpan_student_ba": "student_ba_propensity",
            "netmhcpan_teacher_el": "teacher_el_propensity",
            "netmhcpan_teacher_ba": "teacher_ba_propensity",
        },
        label_column="resp",
        activity_column="r_mean",
        bootstrap_samples=args.bootstrap_samples,
    )

    range_restriction = {
        "teacher_el_rank_le_0p5": int((scores["netmhcpan_el_rank"] <= 0.5).sum()),
        "teacher_el_rank_le_2p0": int((scores["netmhcpan_el_rank"] <= 2.0).sum()),
        "n_peptides": int(len(scores)),
        "note": (
            "The published cohort was pre-screened with NetMHCpan, so nearly every "
            "peptide is already a predicted strong binder. Activation metrics below "
            "are bounded by that restriction and do not rank the binding models."
        ),
    }

    report = {
        "cohort": {
            "source": str(args.cohort.relative_to(REPO_ROOT)),
            "citation": (
                "Visani et al. (2025) T cell receptor specificity landscape revealed "
                "through de novo peptide design. PNAS. doi:10.1073/pnas.2504783122"
            ),
            "allele": ALLELE,
            "n_peptides": int(len(scores)),
            "n_responders": int(scores["resp"].sum()),
            "design_origin": (
                "HERMES-designed NY-ESO-1 variants for the 1G4 TCR, selected on "
                "TCRdock/AlphaFold3 PAE; unrelated to the PDA training corpus"
            ),
        },
        "models": {
            "chao1_mhcflurry_presentation": {
                "checkpoint_sha256": chao1_sha256,
                "description": "frozen chao1 MHC head, MHCflurry-derived calibration",
            },
            "netmhcpan_student": {
                "checkpoint_sha256": student_sha256,
                "model_version": student_manifest["model_version"],
                "description": "five-fold PDA-trained NetMHCpan EL/BA student ensemble",
            },
            "netmhcpan_teacher": {
                "version": str(scores["netmhcpan_version"].iloc[0]),
                "source": str(scores["netmhcpan_source"].iloc[0]),
            },
        },
        "endpoint_1_teacher_imitation": imitation,
        "endpoint_2_experimental_activation": activation,
        "range_restriction": range_restriction,
        "boundaries": [
            "Peptides are pulsed synthetically, so no flanking context exists; the "
            "chao1 cleavage and TAP heads are out of scope for this cohort.",
            "T-cell activation is downstream of presentation and TCR recognition; "
            "no model here is trained on activation labels.",
            "Endpoint 2 is a limits result, not a ranking of the binding models.",
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.output_dir / "scores.csv", index=False)
    (args.output_dir / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    write_report(report, args.output_dir / "REPORT.md")
    print(f"wrote {args.output_dir}")


if __name__ == "__main__":
    main()
