import numpy as np
import pytest

from re_agent.e2e_pls.model import (
    CleavageHead,
    IsotonicCalibrator,
    MhcHead,
    TapHead,
    ThreeHeadModel,
)


def test_isotonic_calibrator_roundtrip():
    x = np.linspace(0, 1, 20)
    y = x**2
    calibrator = IsotonicCalibrator.fit(x, y)
    restored = IsotonicCalibrator.from_dict(calibrator.to_dict())
    np.testing.assert_allclose(calibrator.predict(x), restored.predict(x))


def test_cleavage_head_checkpoint_roundtrip(dev_embeddings, dev_df, tmp_path):
    head = CleavageHead.new(hidden_dim=16)
    head.fit(
        dev_embeddings["n"],
        dev_embeddings["c"],
        dev_df["cleave_n_prob"].values,
        dev_df["cleave_c_prob"].values,
        epochs=20,
    )
    n_pred, c_pred = head.predict(dev_embeddings["n"][:5], dev_embeddings["c"][:5])

    head.save(tmp_path / "cleavage")
    loaded = CleavageHead.load(tmp_path / "cleavage")
    n_pred2, c_pred2 = loaded.predict(dev_embeddings["n"][:5], dev_embeddings["c"][:5])

    np.testing.assert_allclose(n_pred, n_pred2)
    np.testing.assert_allclose(c_pred, c_pred2)
    assert ((0 <= n_pred) & (n_pred <= 1)).all()


def test_tap_head_uncertainty_and_roundtrip(dev_embeddings, dev_df, tmp_path):
    head = TapHead.fit(
        dev_embeddings["mer"], dev_df["tap_log_ic50_relative"].values, n_bootstrap=15
    )
    mean, std = head.predict_with_uncertainty(dev_embeddings["mer"][:5])
    assert mean.shape == (5,)
    assert (std >= 0).all()

    head.save(tmp_path / "tap")
    loaded = TapHead.load(tmp_path / "tap")
    mean2, std2 = loaded.predict_with_uncertainty(dev_embeddings["mer"][:5])
    np.testing.assert_allclose(mean, mean2)
    np.testing.assert_allclose(std, std2)


def test_mhc_head_score_and_roundtrip(dev_embeddings, dev_df, tmp_path):
    head = MhcHead.fit(
        dev_embeddings["mer"],
        dev_df["mhc_percentile"].values,
        dev_df["hla_allele"].values,
        output_dim=16,
    )
    allele = dev_df["hla_allele"].values[0]
    result = head.score(dev_embeddings["mer"][0], allele)
    assert -1.0001 <= result["cosine_similarity"] <= 1.0001
    assert 0 <= result["presentation_propensity"] <= 1

    head.save(tmp_path / "mhc")
    loaded = MhcHead.load(tmp_path / "mhc")
    result2 = loaded.score(dev_embeddings["mer"][0], allele)
    assert result["cosine_similarity"] == pytest.approx(result2["cosine_similarity"])
    assert result["presentation_propensity"] == pytest.approx(result2["presentation_propensity"])


def test_mhc_head_unknown_allele_raises(trained_heads, dev_embeddings):
    with pytest.raises(KeyError):
        trained_heads.mhc.score(dev_embeddings["mer"][0], "HLA-B*07:02")


def test_three_head_model_checkpoint_roundtrip(trained_heads, dev_embeddings, dev_df, tmp_path):
    trained_heads.save(tmp_path / "ckpt")
    loaded = ThreeHeadModel.load(tmp_path / "ckpt")

    assert loaded.model_version == trained_heads.model_version
    assert loaded.dataset_version_hash == trained_heads.dataset_version_hash

    allele = dev_df["hla_allele"].values[0]
    orig = trained_heads.mhc.score(dev_embeddings["mer"][0], allele)
    restored = loaded.mhc.score(dev_embeddings["mer"][0], allele)
    assert orig["presentation_propensity"] == pytest.approx(restored["presentation_propensity"])
