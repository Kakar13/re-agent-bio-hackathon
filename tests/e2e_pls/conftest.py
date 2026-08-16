import numpy as np
import pytest

from re_agent.e2e_pls import fixtures
from re_agent.e2e_pls.encoder import ProteinEncoder
from re_agent.e2e_pls.model import CleavageHead, MhcHead, TapHead, ThreeHeadModel


@pytest.fixture(scope="session")
def dev_df():
    return fixtures.load_dev_fixture()


@pytest.fixture(scope="session")
def mock_client():
    return ProteinEncoder(mode="mock")


def _embed_all(df, client, recipe):
    return np.stack(
        [
            client.embed_pooled(r.peptide, r.n_flank or "", r.c_flank or "", recipe)
            for r in df.itertuples()
        ]
    )


@pytest.fixture(scope="session")
def dev_embeddings(dev_df, mock_client):
    return {
        "n": _embed_all(dev_df, mock_client, "cleave_n"),
        "c": _embed_all(dev_df, mock_client, "cleave_c"),
        "mer": _embed_all(dev_df, mock_client, "mean_9mer"),
    }


@pytest.fixture(scope="session")
def trained_heads(dev_df, dev_embeddings) -> ThreeHeadModel:
    """Small/fast head fit against the dev fixture -- for wiring tests, not accuracy."""
    cleavage = CleavageHead.new(hidden_dim=32)
    cleavage.fit(
        dev_embeddings["n"],
        dev_embeddings["c"],
        dev_df["cleave_n_prob"].values,
        dev_df["cleave_c_prob"].values,
        epochs=40,
        seed=0,
    )
    tap = TapHead.fit(
        dev_embeddings["mer"], dev_df["tap_log_ic50_relative"].values, n_bootstrap=10, seed=0
    )
    mhc = MhcHead.fit(
        dev_embeddings["mer"],
        dev_df["mhc_percentile"].values,
        dev_df["hla_allele"].values,
        output_dim=16,
    )
    return ThreeHeadModel(cleavage=cleavage, tap=tap, mhc=mhc, dataset_version_hash="test-hash")
