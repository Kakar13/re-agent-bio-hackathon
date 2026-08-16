#!/usr/bin/env python3
"""Build the occurrence-level PDA MHC-I processing profile on a detached Modal GPU.

One row per 9-mer occurrence in the PDA de novo corpus, carrying every axis of
the agent-facing profile:

* cleavage  - calibrated N/C probabilities from the frozen chao1 heads
* transport - TAP score and bootstrap uncertainty from chao1
* binding   - NetMHCpan student BA propensity
* presentation - NetMHCpan student EL propensity
* binder class - strong/weak/nonbinder from the predicted EL rank

Student values are out-of-fold: each peptide is scored by the one fold head that
held its parent component out, so the table carries no in-sample optimism. The
deployment ensemble mean is emitted alongside it for comparison.

Two encoders' worth of pooling are needed and they are not interchangeable. The
student was trained on bare 9-mers, while chao1's cleavage and TAP heads were
trained on flanked segments, so each lane is embedded the way its model expects.

    uv run --extra proto modal run --detach scripts/build_mhci_profile_modal.py \
      --source-corpus-run a0201-netmhcpan-pda-full-v1 \
      --student-run a0201-netmhcpan-pda-cv5-v1 \
      --run-name a0201-pda-mhci-profile-v1
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import modal

APP_NAME = "mhci-pda-profile"
ARTIFACT_VOLUME_NAME = "mhci-netmhcpan-artifacts"
WEIGHTS_VOLUME_NAME = "mhci-esm2-weights"
CHAO1_REMOTE_PATH = "/opt/chao1.pkl"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.4,<2.8",
        "fair-esm==2.0.0",
        "httpx>=0.28",
        "numpy>=1.26",
        "pandas>=2.2",
        "pyarrow>=17",
        "pydantic>=2.8",
        "scikit-learn>=1.5",
        "scipy>=1.13",
    )
    .env({"TORCH_HOME": "/weights/torch"})
    .add_local_file("models/chao1/cv5_heads.pkl 2", CHAO1_REMOTE_PATH)
    .add_local_python_source("re_agent")
)

app = modal.App(APP_NAME)
artifact_volume = modal.Volume.from_name(ARTIFACT_VOLUME_NAME, create_if_missing=True)
weights_volume = modal.Volume.from_name(WEIGHTS_VOLUME_NAME, create_if_missing=True)


@app.function(
    image=image,
    gpu="A10G",
    cpu=4,
    memory=32768,
    timeout=6 * 60 * 60,
    volumes={"/artifacts": artifact_volume, "/weights": weights_volume},
)
def build_profile(
    *,
    source_corpus_run: str,
    student_run: str,
    run_name: str,
    pda_designs_bytes: bytes,
    embedding_batch_size: int = 256,
) -> dict:
    import io

    import esm
    import numpy as np
    import pandas as pd
    import torch

    from re_agent.e2e_pls import data
    from re_agent.e2e_pls.netmhcpan_student import (
        AFFINITY_INDEX,
        PROFILE_CAVEATS,
        STRONG_BINDER_RANK,
        WEAK_BINDER_RANK,
        affinity_to_ic50_nm,
        load_student_ensemble,
        propensity_to_rank,
    )
    from re_agent.immuno.e2e_pls_pickle import (
        PEPTIDE_LENGTH,
        _CleavageMLP,
        _interp,
        load_validated_bundle,
    )

    artifact_volume.reload()
    device = torch.device("cuda")

    # ---- peptide lane: teacher labels and out-of-fold student predictions ----
    corpus_path = Path("/artifacts/corpora") / source_corpus_run / "pda_training.parquet"
    if not corpus_path.exists():
        raise FileNotFoundError(corpus_path)
    corpus = pd.read_parquet(corpus_path)
    if corpus["peptide"].duplicated().any():
        raise ValueError("source corpus must hold one row per unique peptide")
    if "cv_fold" not in corpus:
        raise ValueError("source corpus must carry cv_fold for out-of-fold scoring")

    student_dir = Path("/artifacts/models") / student_run
    rows = pd.read_parquet(student_dir / "rows.parquet", columns=["peptide"])
    embeddings = np.load(student_dir / "embeddings.float16.npy")
    if len(rows) != len(embeddings):
        raise ValueError("stored embeddings and rows.parquet disagree on length")
    if not rows["peptide"].astype(str).equals(corpus["peptide"].astype(str)):
        raise ValueError("stored embedding rows are not aligned with the source corpus")

    students, student_sha256, student_manifest = load_student_ensemble(student_dir / "checkpoint")
    fold_models = [model.to(device).eval() for model in students]
    n_outputs = fold_models[0].n_outputs
    has_affinity = n_outputs > AFFINITY_INDEX
    n_peptides = len(corpus)
    fold_predictions = np.empty((len(fold_models), n_peptides, n_outputs), dtype=np.float32)
    embedding_tensor = torch.from_numpy(embeddings.astype(np.float32))
    with torch.inference_mode():
        for offset in range(0, n_peptides, 4096):
            batch = embedding_tensor[offset : offset + 4096].to(device)
            for fold_index, model in enumerate(fold_models):
                fold_predictions[fold_index, offset : offset + len(batch)] = (
                    model(batch).cpu().numpy()
                )

    folds = corpus["cv_fold"].to_numpy(dtype=np.int64)
    if folds.min() < 0 or folds.max() >= len(fold_models):
        raise ValueError("cv_fold values fall outside the five trained folds")
    peptide_index = np.arange(n_peptides)
    out_of_fold = fold_predictions[folds, peptide_index, :]
    ensemble = fold_predictions.mean(axis=0)

    peptide_lane = pd.DataFrame(
        {
            "peptide": corpus["peptide"].astype(str),
            "cv_fold": folds,
            "student_el_propensity_oof": out_of_fold[:, 0].astype(np.float64),
            "student_ba_propensity_oof": out_of_fold[:, 1].astype(np.float64),
            "student_el_propensity_ensemble": ensemble[:, 0].astype(np.float64),
            "student_ba_propensity_ensemble": ensemble[:, 1].astype(np.float64),
        }
    )
    for channel in ("el", "ba"):
        peptide_lane[f"student_predicted_{channel}_rank_oof"] = propensity_to_rank(
            peptide_lane[f"student_{channel}_propensity_oof"].to_numpy()
        )
    if has_affinity:
        peptide_lane["student_ba_affinity_score_oof"] = out_of_fold[
            :, AFFINITY_INDEX
        ].astype(np.float64)
        peptide_lane["student_ba_affinity_score_ensemble"] = ensemble[
            :, AFFINITY_INDEX
        ].astype(np.float64)
        peptide_lane["student_predicted_ba_ic50_nm_oof"] = affinity_to_ic50_nm(
            peptide_lane["student_ba_affinity_score_oof"].to_numpy()
        )
    teacher_columns = [
        column
        for column in (
            "netmhcpan_el_score",
            "netmhcpan_el_rank",
            "netmhcpan_ba_score",
            "netmhcpan_ba_rank",
            "netmhcpan_ba_ic50_nm",
            "binder_class",
            "netmhcpan_allele",
            "netmhcpan_version",
            "parent_component_id",
        )
        if column in corpus
    ]
    peptide_lane = pd.concat([peptide_lane, corpus[teacher_columns].reset_index(drop=True)], axis=1)

    # ---- occurrence lane: every 9-mer window with its flanking context ----
    designs = pd.read_parquet(io.BytesIO(pda_designs_bytes))
    tiled_rows: list[dict] = []
    for design in designs.itertuples(index=False):
        windows = data.tile_protein(
            parent_id=str(design.parent),
            sequence=str(design.seq),
            source_domain="de_novo",
        )
        for window in windows:
            window["pda_release_date"] = str(design.release_date)
            window["pda_novelty_bin"] = str(design.novelty_bin)
        tiled_rows.extend(windows)
    occurrences = pd.DataFrame(tiled_rows)
    occurrences = occurrences.merge(peptide_lane, on="peptide", how="left", validate="many_to_one")
    if occurrences["netmhcpan_el_rank"].isna().any():
        missing = int(occurrences["netmhcpan_el_rank"].isna().sum())
        raise RuntimeError(f"{missing} PDA occurrences have no NetMHCpan label")

    # ---- chao1 lane: flanked-context cleavage, TAP, and MHCflurry presentation ----
    bundle, chao1_sha256 = load_validated_bundle(Path(CHAO1_REMOTE_PATH))
    cleavage_model = _CleavageMLP()
    cleavage_model.load_state_dict(bundle["cleavage"]["state_dict"], strict=True)
    cleavage_model = cleavage_model.eval().to(device)

    tap = bundle["tap"]
    tap_coef = torch.tensor(np.asarray(tap["coef"], dtype=np.float32), device=device)
    tap_intercept = float(tap["intercept"])
    tap_boot_coef = torch.tensor(np.asarray(tap["bootstrap_coef"], dtype=np.float32), device=device)
    tap_boot_intercept = torch.tensor(
        np.asarray(tap["bootstrap_intercept"], dtype=np.float32), device=device
    )
    mhc = bundle["mhc"]
    mhc_weight = torch.tensor(np.asarray(mhc["projection_weight"], dtype=np.float32), device=device)
    mhc_bias = torch.tensor(np.asarray(mhc["projection_bias"], dtype=np.float32), device=device)
    mhc_centroid = torch.tensor(
        np.asarray(mhc["centroids"]["HLA-A*02:01"], dtype=np.float32), device=device
    )

    context_frame = occurrences[["n_flank", "peptide", "c_flank"]].astype(str)
    contexts = pd.MultiIndex.from_frame(context_frame)
    unique_contexts = contexts.unique()
    occurrence_context = unique_contexts.get_indexer(contexts)
    if (occurrence_context < 0).any():
        raise RuntimeError("failed to index every occurrence context")
    unique_frame = unique_contexts.to_frame(index=False)
    unique_frame.columns = ["n_flank", "peptide", "c_flank"]
    segments = (
        unique_frame["n_flank"] + unique_frame["peptide"] + unique_frame["c_flank"]
    ).tolist()
    n_lengths = unique_frame["n_flank"].str.len().to_numpy()
    print(f"embedding {len(segments):,} unique flanked contexts", flush=True)

    encoder, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    encoder = encoder.eval().half().to(device)
    batch_converter = alphabet.get_batch_converter()

    n_contexts = len(segments)
    cleavage_logits = np.empty((n_contexts, 2), dtype=np.float32)
    tap_mean = np.empty(n_contexts, dtype=np.float32)
    tap_uncertainty = np.empty(n_contexts, dtype=np.float32)
    chao1_cosine = np.empty(n_contexts, dtype=np.float32)

    for offset in range(0, n_contexts, embedding_batch_size):
        batch = segments[offset : offset + embedding_batch_size]
        batch_n_lengths = n_lengths[offset : offset + embedding_batch_size]
        _, _, tokens = batch_converter(
            [(f"segment-{offset + index}", segment) for index, segment in enumerate(batch)]
        )
        with torch.inference_mode():
            result = encoder(tokens.to(device), repr_layers=[33], return_contacts=False)
        representations = result["representations"][33].float()

        n_vectors = []
        c_vectors = []
        mer_vectors = []
        for index, n_length in enumerate(batch_n_lengths):
            # Offset by one to skip BOS, matching the reference adapter's slicing.
            residues = representations[index, 1 : 1 + len(batch[index])]
            c_site = int(n_length) + PEPTIDE_LENGTH
            n_vectors.append(
                residues[max(0, int(n_length) - 3) : min(len(residues), int(n_length) + 3)].mean(
                    dim=0
                )
            )
            c_vectors.append(
                residues[max(0, c_site - 3) : min(len(residues), c_site + 3)].mean(dim=0)
            )
            mer_vectors.append(residues[int(n_length) : c_site].mean(dim=0))

        n_array = torch.stack(n_vectors)
        c_array = torch.stack(c_vectors)
        mer_array = torch.stack(mer_vectors)
        with torch.inference_mode():
            logits = cleavage_model(n_array, c_array)
            tap_values = mer_array @ tap_coef + tap_intercept
            tap_boot = mer_array @ tap_boot_coef.T + tap_boot_intercept
            projected = mer_array @ mhc_weight.T + mhc_bias
            projected = projected / projected.norm(dim=1, keepdim=True).clamp(min=1e-6)
            cosine = projected @ mhc_centroid

        span = slice(offset, offset + len(batch))
        cleavage_logits[span] = logits.cpu().numpy()
        tap_mean[span] = tap_values.cpu().numpy()
        tap_uncertainty[span] = tap_boot.std(dim=1).cpu().numpy()
        chao1_cosine[span] = cosine.cpu().numpy()

        if offset % (embedding_batch_size * 200) == 0:
            print(f"  contexts {offset:,}/{n_contexts:,}", flush=True)

    probabilities = 1.0 / (1.0 + np.exp(-cleavage_logits))
    cleavage_n = _interp(bundle["cleavage"]["calibrator_n"], probabilities[:, 0])
    cleavage_c = _interp(bundle["cleavage"]["calibrator_c"], probabilities[:, 1])
    chao1_presentation = _interp(mhc["calibrators"]["HLA-A*02:01"], chao1_cosine)

    occurrences["cleavage_n_probability"] = cleavage_n[occurrence_context]
    occurrences["cleavage_c_probability"] = cleavage_c[occurrence_context]
    occurrences["tap_log_ic50_relative"] = tap_mean[occurrence_context].astype(np.float64)
    occurrences["tap_uncertainty"] = tap_uncertainty[occurrence_context].astype(np.float64)
    occurrences["chao1_mhcflurry_presentation"] = chao1_presentation[occurrence_context]

    # ---- derived profile fields ----
    predicted_el_rank = occurrences["student_predicted_el_rank_oof"].to_numpy()
    occurrences["binder_class_student_oof"] = np.select(
        [predicted_el_rank <= STRONG_BINDER_RANK, predicted_el_rank <= WEAK_BINDER_RANK],
        ["strong", "weak"],
        default="nonbinder",
    )
    occurrences["risk_band_student_oof"] = np.select(
        [predicted_el_rank <= STRONG_BINDER_RANK, predicted_el_rank <= WEAK_BINDER_RANK],
        ["high", "moderate"],
        default="low",
    )
    occurrences["overall_mhci_risk"] = occurrences["student_el_propensity_oof"]
    composite_inputs = np.clip(
        np.stack(
            [
                occurrences["cleavage_n_probability"].to_numpy(),
                occurrences["cleavage_c_probability"].to_numpy(),
                occurrences["student_ba_propensity_oof"].to_numpy(),
            ]
        ),
        1e-6,
        1.0,
    )
    occurrences["composite_processing_risk"] = np.exp(np.log(composite_inputs).mean(axis=0))

    run_dir = Path("/artifacts/profiles") / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    profile_path = run_dir / "pda_mhci_profile.parquet"
    occurrences.to_parquet(profile_path, index=False)

    manifest = {
        "schema_version": "1.0.0",
        "purpose": (
            "occurrence-level MHC-I processing profile over the PDA de novo corpus; "
            "teacher-imitation and processing surrogate outputs, not experimental labels"
        ),
        "grain": "one row per 9-mer occurrence (parent_sequence_id, start)",
        "n_rows": int(len(occurrences)),
        "n_unique_peptides": int(occurrences["peptide"].nunique()),
        "n_unique_contexts": int(n_contexts),
        "n_parent_proteins": int(occurrences["parent_sequence_id"].nunique()),
        "allele": "HLA-A*02:01",
        "profile_axes": {
            "cleavage": ["cleavage_n_probability", "cleavage_c_probability"],
            "transport": ["tap_log_ic50_relative", "tap_uncertainty"],
            "binding": (
                ["student_ba_propensity_oof", "student_predicted_ba_ic50_nm_oof"]
                if has_affinity
                else ["student_ba_propensity_oof"]
            ),
            "presentation": ["student_el_propensity_oof"],
            "binder_class": ["binder_class_student_oof"],
            "overall": ["overall_mhci_risk", "composite_processing_risk"],
            "immunogenicity": None,
        },
        "student_scoring": {
            "mode": "out_of_fold",
            "rule": "each peptide scored by the fold head that held its component out",
            "ensemble_columns": [
                "student_el_propensity_ensemble",
                "student_ba_propensity_ensemble",
            ],
            "checkpoint_sha256": student_sha256,
            "model_version": student_manifest["model_version"],
        },
        "chao1_scoring": {
            "checkpoint_sha256": chao1_sha256,
            "heads": ["cleavage", "tap", "mhc"],
            "pooling": "flanked segment; peptide-residue mean for TAP and MHC",
            "note": (
                "the chao1 MHC head is retained only as an MHCflurry-derived comparison "
                "lane and is not part of the deployed profile"
            ),
        },
        "binder_class_thresholds": {
            "strong": f"predicted EL rank <= {STRONG_BINDER_RANK}",
            "weak": f"predicted EL rank <= {WEAK_BINDER_RANK}",
        },
        "source_corpus_run": source_corpus_run,
        "student_run": student_run,
        "pda_designs_sha256": hashlib.sha256(pda_designs_bytes).hexdigest(),
        "artifact_sha256": hashlib.sha256(profile_path.read_bytes()).hexdigest(),
        "caveats": [
            *PROFILE_CAVEATS,
            "Cleavage and TAP come from the legacy chao1 heads and inherit their "
            "MHCflurry-era training distribution.",
            "No immunogenicity axis is provided; that needs TCR activation or ADA labels.",
        ],
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    artifact_volume.commit()

    summary = {
        "binder_class_counts": occurrences["binder_class_student_oof"].value_counts().to_dict(),
        "teacher_binder_class_counts": occurrences["binder_class"].value_counts().to_dict(),
        "mean_composite_processing_risk": float(occurrences["composite_processing_risk"].mean()),
    }
    return {"manifest": manifest, "summary": summary, "artifact_path": str(profile_path)}


@app.local_entrypoint()
def main(
    source_corpus_run: str = "a0201-netmhcpan-pda-full-v1",
    student_run: str = "a0201-netmhcpan-pda-cv5-v1",
    run_name: str = "a0201-pda-mhci-profile-v1",
    pda_designs: str = "data/processed/pda_designs.parquet",
    embedding_batch_size: int = 256,
) -> None:
    designs_bytes = Path(pda_designs).read_bytes()
    result = build_profile.remote(
        source_corpus_run=source_corpus_run,
        student_run=student_run,
        run_name=run_name,
        pda_designs_bytes=designs_bytes,
        embedding_batch_size=embedding_batch_size,
    )
    print(json.dumps(result, indent=2))
    print(
        "\nInspect cloud artifacts with:\n"
        f"  uv run --extra proto modal volume ls {ARTIFACT_VOLUME_NAME} profiles/{run_name}"
    )
