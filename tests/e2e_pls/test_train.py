import numpy as np

from re_agent.e2e_pls import fixtures, schema
from re_agent.e2e_pls.model import ThreeHeadModel
from re_agent.e2e_pls.train import train


def test_train_writes_loadable_checkpoint_and_metrics(tmp_path):
    output_dir = tmp_path / "checkpoint"
    cache_path = tmp_path / "cache" / "embeddings.dat"

    metrics = train(
        data_path=None,  # dev fixture
        output_dir=output_dir,
        esm3_mode="mock",
        cache_path=cache_path,
        hidden_dim=16,
        cleavage_epochs=20,
        tap_bootstrap=10,
        mhc_output_dim=16,
    )

    assert (output_dir / "manifest.json").exists()
    assert (output_dir / "metrics.json").exists()
    assert "cleavage" in metrics
    assert "tap" in metrics
    assert "mhc" in metrics
    assert metrics["n_train_rows"] > 0

    loaded = ThreeHeadModel.load(output_dir)
    assert loaded.dataset_version_hash  # non-empty: computed from the fixture


def test_train_tolerates_partial_label_coverage(tmp_path):
    """Real Track 1 data has partial coverage per target (e.g. DS613 rows
    have no cleavage label, non-DS613 rows have no measured TAP label) --
    train() must filter per-head rather than crash on NaN.
    """
    df = fixtures.generate_fixture()
    rng = np.random.default_rng(0)
    no_cleavage = rng.choice(df.index, size=len(df) // 3, replace=False)
    no_tap = rng.choice(df.index.difference(no_cleavage), size=len(df) // 3, replace=False)
    df.loc[no_cleavage, ["cleave_n_prob", "cleave_c_prob"]] = np.nan
    df.loc[no_tap, "tap_log_ic50_relative"] = np.nan
    assert schema.validate_dataframe(df).ok

    data_path = tmp_path / "partial.parquet"
    df.to_parquet(data_path, index=False)

    metrics = train(
        data_path=str(data_path),
        output_dir=tmp_path / "checkpoint",
        esm3_mode="mock",
        cache_path=tmp_path / "cache" / "embeddings.dat",
        hidden_dim=16,
        cleavage_epochs=20,
        tap_bootstrap=10,
        mhc_output_dim=16,
    )

    assert metrics["cleavage"]["n_rows"] < metrics["n_test_rows"]
    assert metrics["tap"]["n_rows"] < metrics["n_test_rows"]
    assert metrics["mhc"]["n_rows"] == metrics["n_test_rows"]  # untouched, full coverage
