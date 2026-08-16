"""MHCflurry-independent HLA-A*02:01 corpus construction.

Parent proteins are selected and tiled before NetMHCpan is called. The resulting
training rows retain separate EL and BA channels; BA is also projected into the
legacy E2E-PLS MHC columns so the existing schema and embedding code can inspect
the dataset without treating EL as a second binding score.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.e2e_pls import data, schema
from re_agent.e2e_pls.build_dataset import _load_parent_sources
from re_agent.immuno.netmhcpan import NetMHCpanTeacher, NetMHCpanTeacherConfig

DEFAULT_RAW_DIR = Path("data/raw/netmhcpan-corpus")
DEFAULT_PDA_DESIGNS = Path("data/processed/pda_designs.parquet")
PDA_STUDY_ID = "pda_designs_v1"
PDA_CORPUS_MODE = "pda_only_full"
N_CV_FOLDS = 5


class _UnionFind:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            value, self.parent[value] = self.parent[value], root
        return root

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if left_root < right_root:
            self.parent[right_root] = left_root
        else:
            self.parent[left_root] = right_root


def _sample_by_split(frame: pd.DataFrame, target_rows: int, seed: int) -> pd.DataFrame:
    """Randomly sample before prediction while preserving grouped split ratios."""

    if target_rows <= 0 or len(frame) <= target_rows:
        return frame.reset_index(drop=True)
    ratios = {"train": 0.8, "val": 0.1, "test": 0.1}
    selected: list[int] = []
    remaining_target = target_rows
    for offset, split in enumerate(("train", "val", "test")):
        available = frame.loc[frame["split"] == split]
        requested = (
            remaining_target
            if split == "test"
            else min(len(available), int(target_rows * ratios[split]))
        )
        chosen = available.sample(n=min(requested, len(available)), random_state=seed + offset)
        selected.extend(chosen.index.tolist())
        remaining_target -= len(chosen)
    if remaining_target > 0:
        remaining = frame.drop(index=selected)
        extra = remaining.sample(
            n=min(remaining_target, len(remaining)),
            random_state=seed + 10,
        )
        selected.extend(extra.index.tolist())
    return frame.loc[selected].sample(frac=1, random_state=seed).reset_index(drop=True)


def build_training_pool(
    *,
    target_rows: int,
    n_human_proteins: int,
    n_viral_proteins: int,
    n_bacterial_proteins: int,
    n_de_novo_proteins: int,
    raw_dir: Path = DEFAULT_RAW_DIR,
    seed: int = 0,
) -> tuple[pd.DataFrame, dict]:
    """Fetch, tile, group-split, and neutrally sample proteins without MHC labels."""

    sources, parent_sequences = _load_parent_sources(
        n_human_proteins,
        n_viral_proteins,
        n_bacterial_proteins,
        n_de_novo_proteins,
        data.DEFAULT_UNIPROT_MIN_LENGTH,
        data.DEFAULT_UNIPROT_MAX_LENGTH,
        data.DEFAULT_RCSB_MIN_LENGTH,
        data.DEFAULT_RCSB_MAX_LENGTH,
        raw_dir,
        seed,
    )
    pool = data.build_candidate_pool(sources)
    pool = data.finalize_candidate_columns(pool)
    pool = data.assign_clusters_and_splits(pool, seed=seed)
    pool = _sample_by_split(pool, target_rows, seed)
    metadata = {
        "n_parent_proteins": len(parent_sequences),
        "source_parent_counts": {source: len(rows) for source, rows in sources.items()},
        "candidate_rows_before_sampling": sum(
            max(0, len(sequence) - data.DEFAULT_PEPTIDE_LEN)
            for sequence in parent_sequences.values()
        ),
    }
    return pool, metadata


def attach_netmhcpan_labels(
    pool: pd.DataFrame,
    teacher: NetMHCpanTeacher,
) -> pd.DataFrame:
    """Attach EL/BA labels and populate legacy binding-only MHC columns from BA."""

    labeled = teacher.label(pool)
    labeled["mhc_affinity_nm"] = labeled["netmhcpan_ba_ic50_nm"].astype(np.float32)
    labeled["mhc_percentile"] = labeled["netmhcpan_ba_rank"].astype(np.float32)
    labeled["label_origin"] = "teacher"
    labeled["label_model_version"] = "netmhcpan=4.1;channels=el,ba"
    labeled["generated_at"] = datetime.now(UTC).isoformat()
    validation = schema.validate_dataframe(labeled)
    if not validation.ok:
        raise ValueError(f"NetMHCpan corpus failed E2E-PLS schema: {validation.errors}")
    return labeled


def build_pda_challenge_pool(
    *,
    pda_designs_path: Path,
    target_rows: int,
    seed: int,
    exclude_peptides: set[str] | None = None,
) -> pd.DataFrame:
    """Retile PDA parents into 9-mers and sample without using predictor outputs."""

    if target_rows <= 0:
        return pd.DataFrame()
    if not pda_designs_path.exists():
        raise FileNotFoundError(
            f"PDA designs not found at {pda_designs_path}; run scripts/build_pda_pool.py first"
        )
    designs = pd.read_parquet(pda_designs_path)
    rows = []
    for design in designs.itertuples(index=False):
        tiled = data.tile_protein(
            parent_id=str(design.parent),
            sequence=str(design.seq),
            source_domain="de_novo",
        )
        for row in tiled:
            row.update(
                {
                    "pda_release_date": str(design.release_date),
                    "pda_tm_natural": design.tm_natural,
                    "pda_novelty_bin": str(design.novelty_bin),
                }
            )
        rows.extend(tiled)
    challenge = pd.DataFrame(rows).drop_duplicates("peptide", keep="first")
    if exclude_peptides:
        challenge = challenge.loc[~challenge["peptide"].isin(exclude_peptides)]
    if len(challenge) > target_rows:
        challenge = challenge.sample(n=target_rows, random_state=seed)
    challenge = challenge.reset_index(drop=True)
    challenge["split"] = "challenge"
    challenge["row_id"] = (
        challenge["parent_sequence_id"]
        + ":"
        + challenge["start"].astype(str)
        + "-"
        + challenge["end"].astype(str)
    )
    return challenge


def label_pda_challenge(
    challenge: pd.DataFrame,
    teacher: NetMHCpanTeacher,
) -> pd.DataFrame:
    if challenge.empty:
        return challenge
    return teacher.label(challenge)


def _load_pda_parents(pda_designs_path: Path) -> pd.DataFrame:
    """Load all PDA parents and retain their source metadata for every prediction."""

    if not pda_designs_path.exists():
        raise FileNotFoundError(
            f"PDA designs not found at {pda_designs_path}; run scripts/build_pda_pool.py first"
        )
    designs = pd.read_parquet(pda_designs_path)
    required = {"parent", "seq"}
    missing = required - set(designs)
    if missing:
        raise ValueError(f"PDA designs are missing columns: {sorted(missing)}")
    if designs["parent"].duplicated().any():
        raise ValueError("PDA parent identifiers must be unique")

    parents = pd.DataFrame(
        {
            "parent_sequence_id": designs["parent"].astype(str),
            "sequence": designs["seq"].astype(str).str.upper(),
        }
    )
    invalid = ~parents["sequence"].map(lambda value: set(value) <= schema.CANONICAL_RESIDUES)
    if invalid.any():
        sample = parents.loc[invalid, "parent_sequence_id"].head(5).tolist()
        raise ValueError(f"PDA parents contain noncanonical residues (e.g. {sample})")

    parents["parent_sequence_hash"] = parents["sequence"].map(
        lambda value: hashlib.sha256(value.encode()).hexdigest()[:16]
    )
    for source_column, output_column in (
        ("pdb", "pda_pdb"),
        ("release_date", "pda_release_date"),
        ("tm_natural", "pda_tm_natural"),
        ("tm_partner", "pda_tm_partner"),
        ("seq_bitscore_natural", "pda_seq_bitscore_natural"),
        ("classification", "pda_classification"),
        ("novelty_bin", "pda_novelty_bin"),
        ("in_rcsb_pool", "pda_in_rcsb_pool"),
    ):
        if source_column in designs:
            parents[output_column] = designs[source_column].to_numpy()
    return parents


def _component_fold_map(
    occurrences: pd.DataFrame,
    *,
    n_folds: int = N_CV_FOLDS,
    seed: int = 0,
) -> tuple[dict[str, int], dict[str, str], dict[str, int]]:
    """Connect parents sharing a peptide and balance components by unique peptides."""

    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    parent_ids = sorted(occurrences["parent_sequence_id"].astype(str).unique())
    union_find = _UnionFind(parent_ids)
    first_parent: dict[str, str] = {}
    for peptide, parent_id in zip(
        occurrences["peptide"].astype(str),
        occurrences["parent_sequence_id"].astype(str),
        strict=True,
    ):
        previous = first_parent.setdefault(peptide, parent_id)
        union_find.union(previous, parent_id)

    component_parents: dict[str, list[str]] = {}
    for parent_id in parent_ids:
        component_parents.setdefault(union_find.find(parent_id), []).append(parent_id)
    parent_component: dict[str, str] = {}
    for members in component_parents.values():
        component_id = "pda_component:" + hashlib.sha256(
            "\n".join(members).encode()
        ).hexdigest()[:16]
        parent_component.update(dict.fromkeys(members, component_id))

    peptide_components = (
        occurrences[["peptide", "parent_sequence_id"]]
        .drop_duplicates("peptide")
        .assign(
            component_id=lambda frame: frame["parent_sequence_id"].map(parent_component)
        )
    )
    component_weights = (
        peptide_components.groupby("component_id")["peptide"].nunique().astype(int).to_dict()
    )
    tie_break = lambda component_id: hashlib.sha256(  # noqa: E731
        f"{seed}:{component_id}".encode()
    ).hexdigest()
    ordered_components = sorted(
        component_weights,
        key=lambda component_id: (-component_weights[component_id], tie_break(component_id)),
    )
    fold_tie_order = sorted(
        range(n_folds),
        key=lambda fold: hashlib.sha256(f"{seed}:fold:{fold}".encode()).hexdigest(),
    )
    fold_tie_rank = {fold: rank for rank, fold in enumerate(fold_tie_order)}
    fold_weights = [0] * n_folds
    component_folds: dict[str, int] = {}
    for component_id in ordered_components:
        fold = min(range(n_folds), key=lambda value: (fold_weights[value], fold_tie_rank[value]))
        component_folds[component_id] = fold
        fold_weights[fold] += component_weights[component_id]
    return component_folds, parent_component, component_weights


def _binder_class(rank: pd.Series) -> pd.Series:
    values = pd.to_numeric(rank, errors="raise")
    return pd.Series(
        np.select(
            [values <= 0.5, values <= 2.0],
            ["strong", "weak"],
            default="nonbinder",
        ),
        index=rank.index,
        dtype="string",
    )


def build_full_pda_corpus(
    *,
    pda_designs_path: Path,
    teacher: NetMHCpanTeacher,
    parent_batch_size: int = 20,
    seed: int = 0,
) -> tuple[pd.DataFrame, dict]:
    """Label every PDA parent, deduplicate peptides, and assign leakage-safe folds."""

    if parent_batch_size <= 0:
        raise ValueError("parent_batch_size must be positive")
    parents = _load_pda_parents(pda_designs_path)
    occurrences = teacher.label_parent_sequences(
        parents,
        id_column="parent_sequence_id",
        sequence_column="sequence",
        parents_per_batch=parent_batch_size,
    )
    expected_windows = int((parents["sequence"].str.len() - data.DEFAULT_PEPTIDE_LEN + 1).sum())
    if len(occurrences) != expected_windows:
        raise RuntimeError(
            f"NetMHCpan returned {len(occurrences):,} PDA windows; expected {expected_windows:,}"
        )
    if occurrences["peptide"].str.len().ne(data.DEFAULT_PEPTIDE_LEN).any():
        raise RuntimeError("NetMHCpan returned a non-9-mer for the PDA corpus")

    component_folds, parent_components, component_weights = _component_fold_map(
        occurrences,
        seed=seed,
    )
    occurrences = occurrences.copy()
    occurrences["parent_component_id"] = occurrences["parent_sequence_id"].map(parent_components)
    occurrences["cv_fold"] = occurrences["parent_component_id"].map(component_folds).astype("int8")
    occurrences = occurrences.sort_values(
        ["peptide", "parent_sequence_id", "start", "end"],
        kind="mergesort",
    )

    grouped = occurrences.groupby("peptide", sort=True, observed=True)
    corpus = grouped.first().reset_index()
    corpus["occurrence_count"] = grouped.size().to_numpy(dtype=np.int32)
    corpus["parent_count"] = grouped["parent_sequence_id"].nunique().to_numpy(dtype=np.int32)
    corpus["occurrence_parent_ids"] = grouped["parent_sequence_id"].apply(
        lambda values: json.dumps(sorted(set(values)), separators=(",", ":"))
    ).to_numpy()
    corpus["occurrence_locations"] = grouped.apply(
        lambda frame: json.dumps(
            [
                {
                    "parent_sequence_id": str(row.parent_sequence_id),
                    "start": int(row.start),
                    "end": int(row.end),
                }
                for row in frame.itertuples(index=False)
            ],
            separators=(",", ":"),
        ),
        include_groups=False,
    ).to_numpy()

    now = datetime.now(UTC).isoformat()
    corpus["row_id"] = corpus["peptide"].map(
        lambda peptide: f"pda9:{hashlib.sha256(peptide.encode()).hexdigest()[:16]}"
    )
    corpus["source_domain"] = "de_novo"
    corpus["mhc_affinity_nm"] = corpus["netmhcpan_ba_ic50_nm"].astype(np.float32)
    corpus["mhc_percentile"] = corpus["netmhcpan_ba_rank"].astype(np.float32)
    corpus["hla_allele"] = data.DEFAULT_HLA
    corpus["label_origin"] = "teacher"
    corpus["label_model_version"] = "netmhcpan=4.1;channels=el,ba"
    corpus["source_uri"] = corpus.get(
        "pda_pdb",
        pd.Series("", index=corpus.index, dtype="string"),
    ).fillna("").astype(str).map(
        lambda pdb_id: f"https://www.rcsb.org/structure/{pdb_id}" if pdb_id else ""
    )
    corpus["license"] = data.LICENSE_DE_NOVO
    corpus["generated_at"] = now
    corpus["study_id"] = PDA_STUDY_ID
    corpus["protein_cluster_id"] = corpus["parent_component_id"]
    corpus["peptide_cluster_id"] = corpus["parent_component_id"]
    corpus["split"] = np.select(
        [corpus["cv_fold"].eq(0), corpus["cv_fold"].eq(1)],
        ["test", "val"],
        default="train",
    )
    corpus["encoder_model_id"] = schema.ENCODER_MODEL_ID
    corpus["pooling_recipe"] = "mean_9mer"
    corpus["embedding_cache_key"] = [
        schema.embedding_cache_key(peptide, n_flank or "", c_flank or "")
        for peptide, n_flank, c_flank in zip(
            corpus["peptide"],
            corpus["n_flank"],
            corpus["c_flank"],
            strict=True,
        )
    ]
    corpus["has_structure_track"] = False
    for column in ("cleave_n_prob", "cleave_c_prob", "tap_log_ic50_relative"):
        corpus[column] = np.nan
    corpus["binder_class"] = _binder_class(corpus["netmhcpan_el_rank"])
    corpus["netmhcpan_el_binder_class"] = corpus["binder_class"]
    corpus["netmhcpan_ba_binder_class"] = _binder_class(corpus["netmhcpan_ba_rank"])

    validation = schema.validate_dataframe(corpus)
    if not validation.ok:
        raise ValueError(f"PDA NetMHCpan corpus failed E2E-PLS schema: {validation.errors}")
    observed_folds = sorted(corpus["cv_fold"].unique().tolist())
    if len(component_weights) >= N_CV_FOLDS and observed_folds != list(range(N_CV_FOLDS)):
        raise RuntimeError(f"expected CV folds 0..4, observed {observed_folds}")

    metadata = {
        "n_parent_proteins": len(parents),
        "n_parent_components": len(component_weights),
        "n_windows_before_deduplication": len(occurrences),
        "n_unique_peptides": len(corpus),
        "n_duplicate_occurrences": len(occurrences) - len(corpus),
        "n_shared_parent_peptides": int(corpus["parent_count"].gt(1).sum()),
        "component_unique_peptide_weights": {
            key: int(value) for key, value in sorted(component_weights.items())
        },
    }
    return corpus, metadata


def build_pda_corpus_artifacts(
    *,
    output_dir: Path,
    cache_dir: Path,
    pda_designs_path: Path = DEFAULT_PDA_DESIGNS,
    api_batch_size: int = 500,
    parent_batch_size: int = 20,
    seed: int = 0,
    teacher: NetMHCpanTeacher | None = None,
) -> dict:
    """Write the distinct, full-PDA teacher corpus and its provenance manifest."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if teacher is None:
        teacher = NetMHCpanTeacher(
            Path(cache_dir) / "pda_only_predictions",
            NetMHCpanTeacherConfig(batch_size=api_batch_size),
        )
    corpus, metadata = build_full_pda_corpus(
        pda_designs_path=Path(pda_designs_path),
        teacher=teacher,
        parent_batch_size=parent_batch_size,
        seed=seed,
    )
    corpus_path = output_dir / "pda_training.parquet"
    manifest_path = output_dir / "manifest.json"
    corpus.to_parquet(corpus_path, index=False)
    source_hash = hashlib.sha256(Path(pda_designs_path).read_bytes()).hexdigest()
    manifest = {
        "schema_version": "1.0.0",
        "mode": PDA_CORPUS_MODE,
        "purpose": (
            "NetMHCpan teacher distillation on all PDA parents; "
            "not experimental ground truth"
        ),
        "teacher": {
            "provider": "NetMHCpan",
            "version": "4.1",
            "channels": ["EL", "BA"],
            "transport": "IEDB Tools API",
            "parent_batch_size": parent_batch_size,
        },
        "allele": data.DEFAULT_HLA,
        "peptide_length": data.DEFAULT_PEPTIDE_LEN,
        "selection": (
            "all PDA parents and all 9-mer windows labeled before "
            "exact-peptide deduplication"
        ),
        "deduplication": {
            "key": "peptide",
            "retained_counts": ["occurrence_count", "parent_count"],
            "retained_provenance": ["occurrence_parent_ids", "occurrence_locations"],
        },
        "cross_validation": {
            "n_folds": N_CV_FOLDS,
            "column": "cv_fold",
            "grouping": "connected components of parents sharing any exact 9-mer",
            "balancing_weight": "unique peptide rows",
            "compatibility_split": {"test": [0], "val": [1], "train": [2, 3, 4]},
        },
        "binder_class": {
            "default_column": "binder_class",
            "rank_channel": "EL",
            "thresholds": {"strong": "<=0.5", "weak": "<=2.0", "nonbinder": ">2.0"},
            "experimental_labels": False,
        },
        "n_rows": len(corpus),
        "cv_fold_counts": {
            str(key): int(value)
            for key, value in corpus["cv_fold"].value_counts().sort_index().items()
        },
        "binder_class_counts": {
            str(key): int(value) for key, value in corpus["binder_class"].value_counts().items()
        },
        "source": {
            "path": str(pda_designs_path),
            "sha256": source_hash,
        },
        "dataset_version_hash": schema.dataset_version_hash(corpus),
        "artifact_sha256": _frame_sha256(corpus),
        "seed": seed,
        "created_at": datetime.now(UTC).isoformat(),
        "artifacts": {"pda_training": str(corpus_path)},
        **metadata,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Full PDA NetMHCpan corpus complete -> {output_dir}", flush=True)
    return manifest


def corpus_manifest(
    training: pd.DataFrame,
    challenge: pd.DataFrame,
    source_metadata: dict,
    *,
    requested_training_rows: int,
    requested_challenge_rows: int,
    seed: int,
) -> dict:
    return {
        "schema_version": "1.0.0",
        "purpose": "NetMHCpan teacher distillation; not experimental ground truth",
        "teacher": {
            "provider": "NetMHCpan",
            "version": "4.1",
            "channels": ["EL", "BA"],
            "transport": "IEDB Tools API",
        },
        "allele": data.DEFAULT_HLA,
        "peptide_length": data.DEFAULT_PEPTIDE_LEN,
        "selection": "parent proteins and peptide rows selected before MHC prediction",
        "requested_training_rows": requested_training_rows,
        "requested_challenge_rows": requested_challenge_rows,
        "n_training_rows": len(training),
        "n_challenge_rows": len(challenge),
        "split_counts": {
            key: int(value) for key, value in training["split"].value_counts().items()
        },
        "source_domain_counts": {
            key: int(value) for key, value in training["source_domain"].value_counts().items()
        },
        "dataset_version_hash": schema.dataset_version_hash(training),
        "challenge_sha256": _frame_sha256(challenge),
        "seed": seed,
        "created_at": datetime.now(UTC).isoformat(),
        **source_metadata,
    }


def _frame_sha256(frame: pd.DataFrame) -> str | None:
    if frame.empty:
        return None
    hashed = pd.util.hash_pandas_object(frame.sort_values("row_id"), index=False)
    return hashlib.sha256(hashed.to_numpy(dtype="uint64", copy=False).tobytes()).hexdigest()


def build_corpus_artifacts(
    *,
    output_dir: Path,
    cache_dir: Path,
    target_rows: int,
    pda_challenge_rows: int,
    n_human_proteins: int,
    n_viral_proteins: int,
    n_bacterial_proteins: int,
    n_de_novo_proteins: int,
    pda_designs_path: Path = DEFAULT_PDA_DESIGNS,
    api_batch_size: int = 500,
    seed: int = 0,
) -> dict:
    """Execute acquisition and labeling, then write inspectable artifacts."""

    output_dir = Path(output_dir)
    cache_dir = Path(cache_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pool, source_metadata = build_training_pool(
        target_rows=target_rows,
        n_human_proteins=n_human_proteins,
        n_viral_proteins=n_viral_proteins,
        n_bacterial_proteins=n_bacterial_proteins,
        n_de_novo_proteins=n_de_novo_proteins,
        raw_dir=cache_dir / "sources",
        seed=seed,
    )
    print(
        f"neutral training pool: {len(pool):,} rows; "
        f"splits={pool['split'].value_counts().to_dict()}",
        flush=True,
    )
    teacher = NetMHCpanTeacher(
        cache_dir / "predictions",
        NetMHCpanTeacherConfig(batch_size=api_batch_size),
    )
    training = attach_netmhcpan_labels(pool, teacher)
    challenge_pool = build_pda_challenge_pool(
        pda_designs_path=pda_designs_path,
        target_rows=pda_challenge_rows,
        seed=seed,
        exclude_peptides=set(training["peptide"]),
    )
    print(f"neutral PDA challenge pool: {len(challenge_pool):,} rows", flush=True)
    challenge = label_pda_challenge(challenge_pool, teacher)

    training_path = output_dir / "training.parquet"
    challenge_path = output_dir / "pda_challenge.parquet"
    manifest_path = output_dir / "manifest.json"
    training.to_parquet(training_path, index=False)
    if not challenge.empty:
        challenge.to_parquet(challenge_path, index=False)
    manifest = corpus_manifest(
        training,
        challenge,
        source_metadata,
        requested_training_rows=target_rows,
        requested_challenge_rows=pda_challenge_rows,
        seed=seed,
    )
    manifest["artifacts"] = {
        "training": str(training_path),
        "pda_challenge": str(challenge_path) if not challenge.empty else None,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"NetMHCpan corpus complete -> {output_dir}", flush=True)
    return manifest
