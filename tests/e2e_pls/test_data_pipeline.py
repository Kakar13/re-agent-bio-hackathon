"""Tests for the offline, network-free parts of the Track 1 pipeline.

Live network calls (UniProt, RCSB, DS613) are exercised by actually running
`build_dataset.py`, not here -- a unit suite shouldn't depend on external
services being reachable.
"""

import numpy as np
import pandas as pd

from re_agent.e2e_pls import data, dataset_card, schema


def test_tile_protein_excludes_first_position_only():
    seq = "MKTAYIAKQRQISFVKSHFS"  # length 20
    windows = data.tile_protein("p1", seq, "natural_human")
    # start runs 1..(n-9), so n-9 windows, not n-9+1
    assert len(windows) == len(seq) - 9
    assert windows[0]["start"] == 1
    for w in windows:
        assert seq[w["start"] : w["end"]] == w["peptide"]
        assert w["n_flank"] == seq[max(0, w["start"] - 4) : w["start"]]
        assert w["c_flank"] == seq[w["end"] : w["end"] + 4]


def test_tile_protein_skips_noncanonical_residues():
    seq = "MKTAYIAXQRQISFVKSHFS"  # X at position 7
    windows = data.tile_protein("p1", seq, "natural_human")
    assert all("X" not in w["peptide"] for w in windows)


def test_tile_protein_skips_windows_with_noncanonical_flank():
    # X sits outside every peptide but inside some windows' flanks
    seq = "MKTAYIAXKQRQISFVKSHFS"
    windows = data.tile_protein("p1", seq, "natural_human")
    assert all("X" not in w["n_flank"] and "X" not in w["c_flank"] for w in windows)


def test_build_candidate_pool_dedupes_by_peptide():
    seqs_a = {"p1": "MKTAYIAKQRQISFVKSHFS"}
    seqs_b = {"p2": "MKTAYIAKQRQISFVKSHFS"}  # identical sequence -> identical peptides
    pool = data.build_candidate_pool({"natural_human": seqs_a, "de_novo": seqs_b})
    assert pool["peptide"].is_unique


def test_build_candidate_pool_chunked_matches_unchunked():
    seqs = {
        "p1": "MKTAYIAKQRQISFVKSHFS",
        "p2": "MSTAVLENPGLGRKLSDFGQ",
        "p3": "GAVLIFYWKRHNDQECMSTP",
    }
    full = data.build_candidate_pool({"natural_human": seqs}, chunk_rows=100_000)
    chunked = data.build_candidate_pool({"natural_human": seqs}, chunk_rows=5)
    assert list(full["peptide"]) == list(chunked["peptide"])


def test_quantile_stratified_sample_respects_target_and_spans_range():
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"score": rng.uniform(0, 100, size=1000)})
    sampled = data.quantile_stratified_sample(
        df, target_n=200, score_col="score", n_bins=10, seed=0
    )
    assert 180 <= len(sampled) <= 200
    # spans low and high scores, not just one cluster
    assert sampled["score"].min() < 20
    assert sampled["score"].max() > 80


def test_quantile_stratified_sample_noop_when_pool_smaller_than_target():
    df = pd.DataFrame({"score": [1.0, 2.0, 3.0]})
    sampled = data.quantile_stratified_sample(df, target_n=100, score_col="score")
    assert len(sampled) == 3


def test_inherit_or_assign_splits_keeps_existing_protein_split():
    existing = pd.DataFrame(
        {
            "parent_sequence_id": ["p1", "p1", "p2"],
            "split": ["train", "train", "val"],
            "start": [1, 6, 1],
        }
    )
    new = pd.DataFrame(
        {
            "parent_sequence_id": ["p1", "p3", "p3"],
            "start": [11, 1, 6],
        }
    )
    out = data.inherit_or_assign_splits(new, existing, seed=0)
    assert (out.loc[out["parent_sequence_id"] == "p1", "split"] == "train").all()
    assert out.loc[out["parent_sequence_id"] == "p3", "split"].nunique() == 1
    assert set(out["split"]) <= set(schema.SPLIT_NAMES)


def test_assign_clusters_and_splits_keeps_protein_whole():
    seqs_a = {"p1": "MKTAYIAKQRQISFVKSHFS", "p2": "MSTAVLENPGLGRKLSDFGQ"}
    pool = data.build_candidate_pool({"natural_human": seqs_a})
    pool = data.assign_clusters_and_splits(pool, seed=0)
    splits_per_protein = pool.groupby("protein_cluster_id")["split"].nunique()
    assert (splits_per_protein == 1).all()
    assert set(pool["split"]) <= set(schema.SPLIT_NAMES)


def test_finalize_candidate_columns_produces_valid_schema_subset():
    seqs = {"uniprot_TEST1": "MKTAYIAKQRQISFVKSHFS"}
    pool = data.build_candidate_pool({"natural_human": seqs})
    pool = data.assign_clusters_and_splits(pool, seed=0)
    pool = data.finalize_candidate_columns(pool)
    for col in schema.REQUIRED_COLUMNS:
        assert col in pool.columns, col
    assert (pool["license"] == data.LICENSE_NATURAL).all()
    assert (pool["study_id"] == data.STUDY_BY_DOMAIN["natural_human"]).all()
    assert pool["row_id"].is_unique


def test_build_ds613_rows_shape_and_splits():
    residues = "ACDEFGHIKLMNPQRSTVWY"
    sequences = [f"AAAAAAA{residues[i % len(residues)]}A" for i in range(50)]  # unique 9-mers
    ds613_df = pd.DataFrame({"Sequence": sequences, "log(IC50_relative)": np.linspace(-4, 4, 50)})
    rows = data.build_ds613_rows(ds613_df, seed=0)
    assert len(rows) == 50
    assert rows["row_id"].is_unique
    assert set(rows["split"]) <= set(schema.SPLIT_NAMES)
    assert rows["cleave_n_prob"].isna().all()
    assert rows["label_origin"].eq("measured").all()


def test_impute_tap_labels_fills_missing_and_leaves_measured_untouched():
    from re_agent.e2e_pls import label

    residues = "ACDEFGHIKLMNPQRSTVWY"
    n_measured = 30
    n_missing = 20
    measured_peptides = [f"AAAAAAA{residues[i % 20]}A" for i in range(n_measured)]
    missing_peptides = [f"CCCCCCC{residues[i % 20]}C" for i in range(n_missing)]
    df = pd.DataFrame(
        {
            "peptide": measured_peptides + missing_peptides,
            "tap_log_ic50_relative": [float(i) for i in range(n_measured)] + [np.nan] * n_missing,
            "label_model_version": ["mhcflurry=x"] * (n_measured + n_missing),
            "study_id": ["ds613"] * n_measured + ["uniprot_x"] * n_missing,
        }
    )
    imputed = label.impute_tap_labels(df)

    # measured values are unchanged
    np.testing.assert_array_equal(
        imputed.loc[: n_measured - 1, "tap_log_ic50_relative"].values,
        df.loc[: n_measured - 1, "tap_log_ic50_relative"].values,
    )
    # every missing value is now filled
    assert imputed["tap_log_ic50_relative"].notna().all()
    # imputed rows carry the ridge tag; measured rows do not
    assert imputed.loc[n_measured:, "label_model_version"].str.contains(label.TAP_RIDGE_TAG).all()
    assert (
        not imputed.loc[: n_measured - 1, "label_model_version"]
        .str.contains(label.TAP_RIDGE_TAG)
        .any()
    )


def test_has_measured_tap_flags_only_non_imputed_rows():
    from re_agent.e2e_pls import label

    df = pd.DataFrame(
        {
            "label_model_version": [
                "mhcflurry=2.2.1",
                f"mhcflurry=2.2.1;{label.TAP_RIDGE_TAG}",
                f"pepsickle=x;mhcflurry=y;{label.TAP_RIDGE_TAG}",
                "",
            ]
        }
    )
    mask = label.has_measured_tap(df)
    np.testing.assert_array_equal(mask.values, [True, False, False, True])


def test_dataset_card_hard_gate_rejects_invalid_dataframe(tmp_path):
    from re_agent.e2e_pls import fixtures

    df = fixtures.generate_fixture()
    df.loc[df.index[0], "end"] = df.loc[df.index[0], "end"] + 5  # break coordinate consistency
    try:
        dataset_card.build_dataset_card(df, tmp_path, sources=[], build_config={})
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_dataset_card_writes_manifest_and_card_for_valid_dataframe(tmp_path):
    from re_agent.e2e_pls import fixtures

    df = fixtures.generate_fixture()
    provenance = data.SourceProvenance(
        name="test-source",
        url="https://example.com",
        retrieved_at=data.now_iso(),
        sha256="abc",
        license="CC0",
        n_records=len(df),
    )
    manifest = dataset_card.build_dataset_card(
        df, tmp_path, sources=[provenance], build_config={"seed": 0}
    )

    assert (tmp_path / "dataset.parquet").exists()
    assert (tmp_path / "dataset_manifest.json").exists()
    assert (tmp_path / "dataset_card.md").exists()
    assert manifest["n_rows"] == len(df)
    assert manifest["dataset_version_hash"] == schema.dataset_version_hash(df)

    reloaded = pd.read_parquet(tmp_path / "dataset.parquet")
    assert len(reloaded) == len(df)
    assert schema.validate_dataframe(reloaded).ok
