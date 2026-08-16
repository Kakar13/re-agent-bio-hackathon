import numpy as np
import pytest

from re_agent.e2e_pls.esm3_modal import EmbeddingCache, ESM3Client
from re_agent.e2e_pls.schema import EMBEDDING_DIM


def test_embed_residues_shape_and_determinism(mock_client):
    seq = "MKTAYIAKQ"
    r1 = mock_client.embed_residues(seq)
    r2 = mock_client.embed_residues(seq)
    assert r1.shape == (len(seq), EMBEDDING_DIM)
    assert np.allclose(r1, r2)


def test_embed_residues_empty_sequence(mock_client):
    r = mock_client.embed_residues("")
    assert r.shape == (0, EMBEDDING_DIM)


@pytest.mark.parametrize("recipe", ["cleave_n", "cleave_c", "mean_9mer"])
def test_embed_pooled_shape(mock_client, recipe):
    vec = mock_client.embed_pooled("AYIAKQRQI", "MKT", "SFVK", recipe)
    assert vec.shape == (EMBEDDING_DIM,)


def test_embed_pooled_differs_by_recipe(mock_client):
    n_vec = mock_client.embed_pooled("AYIAKQRQI", "MKT", "SFVK", "cleave_n")
    c_vec = mock_client.embed_pooled("AYIAKQRQI", "MKT", "SFVK", "cleave_c")
    assert not np.allclose(n_vec, c_vec)


def test_generate_masked_only_changes_masked_positions(mock_client):
    seq = "MKTAYIAKQ"
    candidates = mock_client.generate_masked(seq, [2, 5], num_candidates=5)
    assert len(candidates) == 5
    for cand in candidates:
        assert len(cand) == len(seq)
        for i, (a, b) in enumerate(zip(seq, cand, strict=True)):
            if i not in (2, 5):
                assert a == b, f"position {i} should be untouched"
        assert cand[2] != seq[2]
        assert cand[5] != seq[5]


def test_sequence_log_likelihood_deterministic(mock_client):
    seq = "MKTAYIAKQ"
    assert mock_client.sequence_log_likelihood(seq) == mock_client.sequence_log_likelihood(seq)


def test_sequence_log_likelihood_varies_by_sequence(mock_client):
    a = mock_client.sequence_log_likelihood("MKTAYIAKQ")
    b = mock_client.sequence_log_likelihood("GGGGGGGGG")
    assert a != b


def test_embedding_cache_roundtrip_and_growth(tmp_path):
    cache = EmbeddingCache(tmp_path / "cache.dat", dim=8, initial_capacity=2)
    v1 = np.arange(8, dtype="float32")
    v2 = np.arange(8, 16, dtype="float32")
    v3 = np.arange(16, 24, dtype="float32")

    assert cache.get("a") is None
    cache.put("a", v1)
    cache.put("b", v2)
    cache.put("c", v3)  # forces growth beyond initial_capacity=2

    assert np.allclose(cache.get("a"), v1)
    assert np.allclose(cache.get("b"), v2)
    assert np.allclose(cache.get("c"), v3)
    assert "a" in cache
    assert "z" not in cache


def test_embedding_cache_persists_across_instances(tmp_path):
    path = tmp_path / "cache.dat"
    v = np.arange(8, dtype="float32")
    EmbeddingCache(path, dim=8, initial_capacity=4).put("k", v)
    reopened = EmbeddingCache(path, dim=8, initial_capacity=4)
    assert np.allclose(reopened.get("k"), v)


def test_get_or_compute_only_calls_compute_once(tmp_path):
    cache = EmbeddingCache(tmp_path / "cache.dat", dim=8, initial_capacity=4)
    calls = []

    def compute():
        calls.append(1)
        return np.ones(8, dtype="float32")

    v1 = cache.get_or_compute("k", compute)
    v2 = cache.get_or_compute("k", compute)
    assert np.allclose(v1, v2)
    assert len(calls) == 1


def test_client_with_cache_matches_uncached(tmp_path):
    cache = EmbeddingCache(tmp_path / "cache.dat")
    cached_client = ESM3Client(mode="mock", cache=cache)
    uncached_client = ESM3Client(mode="mock")
    v_cached = cached_client.embed_pooled("AYIAKQRQI", "MKT", "SFVK", "mean_9mer")
    v_uncached = uncached_client.embed_pooled("AYIAKQRQI", "MKT", "SFVK", "mean_9mer")
    assert np.allclose(v_cached, v_uncached)
