"""CLI: fit/calibrate the three heads and write a checkpoint + metrics.

    uv run python -m re_agent.e2e_pls.train --data path/to/dataset.parquet
    uv run python -m re_agent.e2e_pls.train                      # dev fixture

Reads the Track 1 parquet if `--data` is given (validated against
`schema.py`), otherwise falls back to the checked-in dev fixture so this
runs standalone before Track 1's dataset exists. Trains on the `train`
split, calibrates on `train`'s held-out slice (inside `CleavageHead.fit`),
and reports held-out metrics on the `test` split.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score

from re_agent.e2e_pls import fixtures, schema
from re_agent.e2e_pls.encoder import EmbeddingCache, ProteinEncoder
from re_agent.e2e_pls.model import CleavageHead, MhcHead, TapHead, ThreeHeadModel

DEFAULT_OUTPUT_DIR = Path("results/e2e_pls/checkpoints/latest")
DEFAULT_CACHE_PATH = Path("results/e2e_pls/embedding_cache/embeddings.dat")


def _load_dataset(data_path: str | None) -> pd.DataFrame:
    if data_path is None:
        return fixtures.load_dev_fixture()
    df = pd.read_parquet(data_path)
    result = schema.validate_dataframe(df)
    if not result.ok:
        raise SystemExit(f"dataset at {data_path} failed schema validation: {result.errors}")
    for warning in result.warnings:
        print(f"warning: {warning}")
    return df


def _embed_dataset(df: pd.DataFrame, client: ProteinEncoder) -> dict[str, np.ndarray]:
    n_vecs, c_vecs, mer_vecs = [], [], []
    for row in df.itertuples():
        n_flank, c_flank = row.n_flank or "", row.c_flank or ""
        n_vecs.append(client.embed_pooled(row.peptide, n_flank, c_flank, "cleave_n"))
        c_vecs.append(client.embed_pooled(row.peptide, n_flank, c_flank, "cleave_c"))
        mer_vecs.append(client.embed_pooled(row.peptide, n_flank, c_flank, "mean_9mer"))
    return {"n": np.stack(n_vecs), "c": np.stack(c_vecs), "mer": np.stack(mer_vecs)}


def _labeled_mask(df: pd.DataFrame, base_mask: np.ndarray, *columns: str) -> np.ndarray:
    """`base_mask` narrowed to rows where every one of `columns` is non-null.

    Real data legitimately has partial label coverage -- e.g. DS613 rows
    have no parent-protein context so no cleavage label, and non-DS613
    rows have no measured TAP label -- so each head trains/evaluates only
    on the rows that actually carry its target(s).
    """
    mask = base_mask.copy()
    for col in columns:
        mask &= df[col].notna().values
    return mask


def _bootstrap_ci(fn, *arrays: np.ndarray, n_boot: int = 200, seed: int = 0) -> list[float]:
    rng = np.random.default_rng(seed)
    n = len(arrays[0])
    stats = [fn(*(a[rng.integers(0, n, n)] for a in arrays)) for _ in range(n_boot)]
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return [float(lo), float(hi)]


def _evaluate(
    df: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    test_mask: np.ndarray,
    cleavage: CleavageHead,
    tap: TapHead,
    mhc: MhcHead,
) -> dict:
    if test_mask.sum() == 0:
        return {"warning": "no test-split rows; skipped held-out evaluation"}

    metrics: dict = {"n_test_rows": int(test_mask.sum())}

    cleavage_mask = _labeled_mask(df, test_mask, "cleave_n_prob", "cleave_c_prob")
    if cleavage_mask.sum() == 0:
        metrics["cleavage"] = {"warning": "no test rows with cleavage labels"}
    else:
        n_pred, c_pred = cleavage.predict(
            embeddings["n"][cleavage_mask], embeddings["c"][cleavage_mask]
        )
        cleave_n_true = df["cleave_n_prob"].values[cleavage_mask]
        cleave_c_true = df["cleave_c_prob"].values[cleavage_mask]
        metrics["cleavage"] = {
            "n_terminus_auprc": float(average_precision_score(cleave_n_true >= 0.5, n_pred)),
            "c_terminus_auprc": float(average_precision_score(cleave_c_true >= 0.5, c_pred)),
            "n_terminus_mse": float(np.mean((n_pred - cleave_n_true) ** 2)),
            "c_terminus_mse": float(np.mean((c_pred - cleave_c_true) ** 2)),
            "n_rows": int(cleavage_mask.sum()),
        }

    from re_agent.e2e_pls.label import has_measured_tap

    # TAP evaluation stays measured-only for the same reason training does:
    # scoring against imputed rows would just verify that we agree with our own imputer.
    tap_mask = _labeled_mask(df, test_mask, "tap_log_ic50_relative") & has_measured_tap(df).values
    if tap_mask.sum() == 0:
        metrics["tap"] = {"warning": "no test rows with measured TAP labels"}
    else:
        tap_pred, tap_std = tap.predict_with_uncertainty(embeddings["mer"][tap_mask])
        tap_true = df["tap_log_ic50_relative"].values[tap_mask]
        metrics["tap"] = {
            "spearman": float(spearmanr(tap_pred, tap_true).statistic),
            "spearman_95ci": _bootstrap_ci(
                lambda p, t: spearmanr(p, t).statistic, tap_pred, tap_true
            ),
            "rmse": float(np.sqrt(np.mean((tap_pred - tap_true) ** 2))),
            "mean_bootstrap_uncertainty": float(tap_std.mean()),
            "n_rows": int(tap_mask.sum()),
        }

    mhc_mask = _labeled_mask(df, test_mask, "mhc_percentile")
    if mhc_mask.sum() == 0:
        metrics["mhc"] = {"warning": "no test rows with MHC labels"}
    else:
        hla_alleles = df["hla_allele"].values[mhc_mask]
        mer_vecs = embeddings["mer"][mhc_mask]
        percentile_true = df["mhc_percentile"].values[mhc_mask]
        propensity_pred = np.array(
            [
                mhc.score(mer_vecs[i], hla_alleles[i])["presentation_propensity"]
                for i in range(mhc_mask.sum())
            ]
        )
        target_true = np.clip(1 - percentile_true / 100, 0, 1)
        top_true = percentile_true <= np.percentile(percentile_true, 5)
        top_pred = propensity_pred >= np.percentile(propensity_pred, 95)
        metrics["mhc"] = {
            "spearman": float(spearmanr(propensity_pred, target_true).statistic),
            "top5pct_agreement": float((top_true & top_pred).sum() / max(1, top_true.sum())),
            "n_rows": int(mhc_mask.sum()),
        }

    return metrics


def train(
    data_path: str | None = None,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    encoder_mode: str = "mock",
    cache_path: Path | None = DEFAULT_CACHE_PATH,
    hidden_dim: int = 128,
    cleavage_epochs: int = 300,
    tap_alpha: float = 10.0,
    tap_bootstrap: int = 50,
    mhc_output_dim: int = 64,
    seed: int = 0,
) -> dict:
    t0 = time.monotonic()
    df = _load_dataset(data_path)

    cache = EmbeddingCache(cache_path) if cache_path else None
    client = ProteinEncoder(mode=encoder_mode, cache=cache)
    embeddings = _embed_dataset(df, client)

    train_mask = (df["split"] == "train").values
    test_mask = (df["split"] == "test").values

    from re_agent.e2e_pls.label import has_measured_tap

    cleavage_train_mask = _labeled_mask(df, train_mask, "cleave_n_prob", "cleave_c_prob")
    # TAP: train only on measured rows (DS613). Rows with imputed TAP were filled
    # by our own ridge in label.impute_tap_labels(); training on them would just
    # recover the imputer and teach the head nothing.
    tap_train_mask = (
        _labeled_mask(df, train_mask, "tap_log_ic50_relative") & has_measured_tap(df).values
    )
    mhc_train_mask = _labeled_mask(df, train_mask, "mhc_percentile")

    cleavage = CleavageHead.new(hidden_dim=hidden_dim)
    fit_info = cleavage.fit(
        embeddings["n"][cleavage_train_mask],
        embeddings["c"][cleavage_train_mask],
        df["cleave_n_prob"].values[cleavage_train_mask],
        df["cleave_c_prob"].values[cleavage_train_mask],
        epochs=cleavage_epochs,
        seed=seed,
    )
    tap = TapHead.fit(
        embeddings["mer"][tap_train_mask],
        df["tap_log_ic50_relative"].values[tap_train_mask],
        alpha=tap_alpha,
        n_bootstrap=tap_bootstrap,
        seed=seed,
    )
    mhc = MhcHead.fit(
        embeddings["mer"][mhc_train_mask],
        df["mhc_percentile"].values[mhc_train_mask],
        df["hla_allele"].values[mhc_train_mask],
        output_dim=mhc_output_dim,
    )

    metrics = _evaluate(df, embeddings, test_mask, cleavage, tap, mhc)
    metrics["cleavage_fit"] = fit_info
    metrics["n_train_rows"] = int(train_mask.sum())
    metrics["duration_s"] = time.monotonic() - t0

    heads = ThreeHeadModel(
        cleavage=cleavage, tap=tap, mhc=mhc, dataset_version_hash=schema.dataset_version_hash(df)
    )
    output_dir = Path(output_dir)
    heads.save(output_dir)
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    return metrics


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data", default=None, help="path to Track 1 parquet; omit to use the dev fixture"
    )
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--encoder-mode", choices=["esm2", "mock", "modal"], default="esm2")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--hidden-dim", type=int, default=128)
    parser.add_argument("--cleavage-epochs", type=int, default=300)
    parser.add_argument("--tap-alpha", type=float, default=10.0)
    parser.add_argument("--tap-bootstrap", type=int, default=50)
    parser.add_argument("--mhc-output-dim", type=int, default=64)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    metrics = train(
        data_path=args.data,
        output_dir=Path(args.output_dir),
        encoder_mode=args.encoder_mode,
        cache_path=None if args.no_cache else DEFAULT_CACHE_PATH,
        hidden_dim=args.hidden_dim,
        cleavage_epochs=args.cleavage_epochs,
        tap_alpha=args.tap_alpha,
        tap_bootstrap=args.tap_bootstrap,
        mhc_output_dim=args.mhc_output_dim,
        seed=args.seed,
    )
    print(f"saved checkpoint to {args.output_dir}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
