#!/usr/bin/env python3
"""Embed a cloud NetMHCpan corpus and train its student head on Modal.

Run this only after ``build_netmhcpan_corpus_modal.py`` completes:

    uv run --extra proto modal run --detach scripts/train_mhci_student_modal.py \
      --corpus-run-name a0201-netmhcpan-10k \
      --run-name a0201-netmhcpan-student-10k
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import modal

APP_NAME = "mhci-netmhcpan-student"
ARTIFACT_VOLUME_NAME = "mhci-netmhcpan-artifacts"
WEIGHTS_VOLUME_NAME = "mhci-esm2-weights"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.4,<2.8",
        "fair-esm==2.0.0",
        "numpy>=1.26",
        "pandas>=2.2",
        "pyarrow>=17",
        "scikit-learn>=1.5",
        "scipy>=1.13",
    )
    .env({"TORCH_HOME": "/weights/torch"})
    .add_local_python_source("re_agent")
)

app = modal.App(APP_NAME)
weights_volume = modal.Volume.from_name(WEIGHTS_VOLUME_NAME, create_if_missing=True)
artifact_volume = modal.Volume.from_name(ARTIFACT_VOLUME_NAME, create_if_missing=True)


@app.function(
    image=image,
    gpu="A10G",
    cpu=4,
    memory=32768,
    timeout=6 * 60 * 60,
    volumes={"/weights": weights_volume, "/artifacts": artifact_volume},
)
def embed_and_train(
    *,
    corpus_run_name: str,
    run_name: str,
    embedding_batch_size: int = 256,
    train_batch_size: int = 512,
    epochs: int = 60,
    hidden_dim: int = 256,
    seed: int = 0,
) -> dict:
    import esm
    import numpy as np
    import pandas as pd
    import torch

    from re_agent.e2e_pls.netmhcpan_student import (
        StudentTrainConfig,
        save_student_checkpoint,
        save_student_cv_checkpoints,
        train_student,
        train_student_cv,
    )

    corpus_dir = Path("/artifacts/corpora") / corpus_run_name
    training_path = corpus_dir / "training.parquet"
    pda_training_path = corpus_dir / "pda_training.parquet"
    challenge_path = corpus_dir / "pda_challenge.parquet"
    artifact_volume.reload()
    if not training_path.exists():
        training_path = pda_training_path
    if not training_path.exists():
        raise FileNotFoundError(
            f"cloud corpus has neither training.parquet nor pda_training.parquet: {corpus_dir}"
        )
    training = pd.read_parquet(training_path)
    frames = [training]
    if challenge_path.exists():
        frames.append(pd.read_parquet(challenge_path))
    frame = pd.concat(frames, ignore_index=True)
    corpus_bytes = training_path.read_bytes()
    if challenge_path.exists():
        corpus_bytes += challenge_path.read_bytes()
    corpus_sha256 = hashlib.sha256(corpus_bytes).hexdigest()

    device = torch.device("cuda")
    encoder, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    encoder = encoder.eval().half().to(device)
    batch_converter = alphabet.get_batch_converter()

    unique_peptides = list(dict.fromkeys(frame["peptide"].astype(str)))
    unique_embeddings = np.empty((len(unique_peptides), 1280), dtype=np.float16)
    for offset in range(0, len(unique_peptides), embedding_batch_size):
        peptides = unique_peptides[offset : offset + embedding_batch_size]
        _, _, tokens = batch_converter(
            [(f"peptide-{offset + index}", peptide) for index, peptide in enumerate(peptides)]
        )
        with torch.inference_mode():
            result = encoder(
                tokens.to(device),
                repr_layers=[33],
                return_contacts=False,
            )
        pooled = result["representations"][33][:, 1:10].mean(dim=1)
        unique_embeddings[offset : offset + len(peptides)] = pooled.cpu().numpy()

    lookup = {peptide: index for index, peptide in enumerate(unique_peptides)}
    row_indices = np.fromiter(
        (lookup[peptide] for peptide in frame["peptide"].astype(str)),
        dtype=np.int64,
        count=len(frame),
    )
    embeddings = unique_embeddings[row_indices]

    run_dir = Path("/artifacts/models") / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    np.save(run_dir / "embeddings.float16.npy", embeddings)
    frame.to_parquet(run_dir / "rows.parquet", index=False)
    embedding_manifest = {
        "encoder_model_id": "esm2_t33_650M_UR50D",
        "layer": 33,
        "pooling": "mean over 9 peptide residues; BOS/EOS excluded",
        "dtype": "float16",
        "shape": list(embeddings.shape),
        "n_unique_peptides": len(unique_peptides),
        "corpus_sha256": corpus_sha256,
        "source_corpus_run": corpus_run_name,
    }
    (run_dir / "embedding_manifest.json").write_text(
        json.dumps(embedding_manifest, indent=2) + "\n"
    )

    train_config = StudentTrainConfig(
        batch_size=train_batch_size,
        epochs=epochs,
        hidden_dim=hidden_dim,
        seed=seed,
    )
    if "cv_fold" in frame.columns:
        students, metrics = train_student_cv(
            embeddings.astype(np.float32),
            frame,
            train_config,
            device=device,
        )
        save_student_cv_checkpoints(
            run_dir / "checkpoint",
            students,
            metrics,
            corpus_sha256=corpus_sha256,
        )
        training_mode = "five_fold_cv"
    else:
        student, metrics = train_student(
            embeddings.astype(np.float32),
            frame,
            train_config,
            device=device,
        )
        save_student_checkpoint(
            run_dir / "checkpoint",
            student,
            metrics,
            corpus_sha256=corpus_sha256,
        )
        training_mode = "split"
    weights_volume.commit()
    artifact_volume.commit()
    return {
        "run_name": run_name,
        "corpus_run_name": corpus_run_name,
        "training_mode": training_mode,
        "artifact_volume": ARTIFACT_VOLUME_NAME,
        "artifact_path": f"/models/{run_name}",
        "embedding_manifest": embedding_manifest,
        "metrics": metrics,
    }


@app.local_entrypoint()
def main(
    corpus_run_name: str = "a0201-netmhcpan-10k",
    run_name: str = "a0201-netmhcpan-student-10k",
    embedding_batch_size: int = 256,
    train_batch_size: int = 512,
    epochs: int = 60,
    hidden_dim: int = 256,
    seed: int = 0,
) -> None:
    result = embed_and_train.remote(
        corpus_run_name=corpus_run_name,
        run_name=run_name,
        embedding_batch_size=embedding_batch_size,
        train_batch_size=train_batch_size,
        epochs=epochs,
        hidden_dim=hidden_dim,
        seed=seed,
    )
    print(json.dumps(result, indent=2))
    print(
        "\nInspect cloud artifacts with:\n"
        f"  uv run --extra proto modal volume ls {ARTIFACT_VOLUME_NAME} models/{run_name}"
    )
