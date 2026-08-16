import pandas as pd

from re_agent.e2e_pls import schema


def test_dev_fixture_passes_validation(dev_df):
    result = schema.validate_dataframe(dev_df)
    assert result.ok, result.errors


def test_duplicate_row_id_detected(dev_df):
    dupes = pd.concat([dev_df, dev_df.iloc[:2]], ignore_index=True)
    result = schema.validate_dataframe(dupes)
    assert not result.ok
    assert any("duplicate row_id" in e for e in result.errors)


def test_coordinate_mismatch_detected(dev_df):
    bad = dev_df.copy()
    bad.loc[bad.index[0], "end"] = bad.loc[bad.index[0], "end"] + 1
    result = schema.validate_dataframe(bad)
    assert not result.ok
    assert any("end - start" in e for e in result.errors)


def test_length_mismatch_detected(dev_df):
    bad = dev_df.copy()
    bad.loc[bad.index[0], "length"] = 99
    result = schema.validate_dataframe(bad)
    assert not result.ok
    assert any("length != len(peptide)" in e for e in result.errors)


def test_noncanonical_residue_detected(dev_df):
    bad = dev_df.copy()
    peptide = bad.loc[bad.index[0], "peptide"]
    bad.loc[bad.index[0], "peptide"] = "X" + peptide[1:]
    result = schema.validate_dataframe(bad)
    assert not result.ok
    assert any("noncanonical" in e for e in result.errors)


def test_invalid_enum_detected(dev_df):
    bad = dev_df.copy()
    bad.loc[bad.index[0], "split"] = "holdout"
    result = schema.validate_dataframe(bad)
    assert not result.ok
    assert any("split" in e for e in result.errors)


def test_split_leakage_detected(dev_df):
    bad = dev_df.copy()
    # force one row's protein_cluster_id to appear in two different splits
    other_split = next(s for s in schema.SPLIT_NAMES if s != bad.loc[bad.index[0], "split"])
    target_cluster = bad.loc[bad.index[0], "protein_cluster_id"]
    same_cluster_mask = bad["protein_cluster_id"] == target_cluster
    bad.loc[bad[same_cluster_mask].index[:1], "split"] = other_split
    result = schema.validate_dataframe(bad)
    assert not result.ok
    assert any("split-group overlap" in e for e in result.errors)


def test_missing_target_is_warning_not_error(dev_df):
    partial = dev_df.copy()
    partial.loc[partial.index[0], "cleave_n_prob"] = None
    result = schema.validate_dataframe(partial)
    assert result.ok
    assert any("cleave_n_prob" in w for w in result.warnings)


def test_dataset_version_hash_deterministic_and_order_independent(dev_df):
    h1 = schema.dataset_version_hash(dev_df)
    h2 = schema.dataset_version_hash(dev_df.sample(frac=1, random_state=1).reset_index(drop=True))
    assert h1 == h2


def test_dataset_version_hash_changes_with_content(dev_df):
    h1 = schema.dataset_version_hash(dev_df)
    mutated = dev_df.copy()
    mutated.loc[mutated.index[0], "cleave_n_prob"] = 0.999999
    h2 = schema.dataset_version_hash(mutated)
    assert h1 != h2


def test_embedding_cache_key_deterministic():
    k1 = schema.embedding_cache_key("AAIKLMNPQ", "GSHM", "STAV")
    k2 = schema.embedding_cache_key("AAIKLMNPQ", "GSHM", "STAV")
    k3 = schema.embedding_cache_key("AAIKLMNPQ", "GSHM", "STAX")
    assert k1 == k2
    assert k1 != k3
