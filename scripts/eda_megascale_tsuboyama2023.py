"""EDA for LiteFold/MegaScale-Tsuboyama2023. Writes a JSON summary for the canvas."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("data/raw/megascale_tsuboyama2023")
OUT = Path("results/megascale_tsuboyama2023/eda_summary.json")

COLS = [
    "record_id",
    "record_type",
    "figure",
    "table_name",
    "protein_id",
    "mutation_name",
    "sequence",
    "sequence_length",
    "position",
    "position_1",
    "position_2",
    "wt_aa",
    "mut_aa",
    "aa1",
    "aa2",
    "delta_g",
    "reconstructed_delta_g",
    "thermodynamic_delta_g",
    "previous_delta_g",
    "target",
    "score_value",
    "split_bucket",
    "dssp",
    "burial",
    "sasa_sc",
    "pdb_id",
    "wt_name",
    "aa_seq",
    "dg",
    "deltag",
]


def hist_counts(series: pd.Series, bins: list[float]) -> list[dict]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return []
    counts, edges = np.histogram(values, bins=bins)
    out = []
    for i, count in enumerate(counts):
        lo, hi = edges[i], edges[i + 1]
        label = f"{lo:g}–{hi:g}"
        out.append({"label": label, "count": int(count)})
    return out


def describe(series: pd.Series) -> dict:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return {"n": 0}
    return {
        "n": int(values.size),
        "mean": float(values.mean()),
        "std": float(values.std()),
        "min": float(values.min()),
        "p05": float(values.quantile(0.05)),
        "p25": float(values.quantile(0.25)),
        "median": float(values.median()),
        "p75": float(values.quantile(0.75)),
        "p95": float(values.quantile(0.95)),
        "max": float(values.max()),
        "zeros": int((values == 0).sum()),
    }


def main() -> None:
    shards = sorted((ROOT / "data").glob("*.parquet"))
    frames = []
    for path in shards:
        split = "test" if "test" in path.name else "train"
        frame = pd.read_parquet(path, columns=COLS)
        frame["split"] = split
        frames.append(frame)
    df = pd.concat(frames, ignore_index=True)

    numeric_cols = [
        "sequence_length",
        "position",
        "position_1",
        "position_2",
        "delta_g",
        "reconstructed_delta_g",
        "thermodynamic_delta_g",
        "previous_delta_g",
        "target",
        "score_value",
        "burial",
        "sasa_sc",
        "dg",
        "deltag",
        "split_bucket",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    n = len(df)
    split_counts = df["split"].value_counts().to_dict()
    type_counts = df["record_type"].value_counts().to_dict()
    figure_counts = df["figure"].value_counts().to_dict()

    type_by_split = (
        df.groupby(["record_type", "split"]).size().unstack(fill_value=0).to_dict()
    )

    # Missingness of modeling-relevant columns, overall and by type
    model_cols = [
        "protein_id",
        "sequence",
        "sequence_length",
        "mutation_name",
        "wt_aa",
        "mut_aa",
        "aa1",
        "aa2",
        "delta_g",
        "score_value",
        "target",
        "reconstructed_delta_g",
        "thermodynamic_delta_g",
        "dssp",
        "burial",
    ]
    missing_overall = {
        col: {"nulls": int(df[col].isna().sum()), "pct": float(df[col].isna().mean())}
        for col in model_cols
    }
    missing_by_type = {}
    for rec_type, group in df.groupby("record_type"):
        missing_by_type[rec_type] = {
            col: round(float(group[col].isna().mean()), 4) for col in model_cols
        }

    # Protein overlap / leakage
    train_proteins = set(df.loc[df["split"] == "train", "protein_id"].dropna())
    test_proteins = set(df.loc[df["split"] == "test", "protein_id"].dropna())
    protein_overlap = train_proteins & test_proteins

    train_seq = set(df.loc[df["split"] == "train", "sequence"].dropna())
    test_seq = set(df.loc[df["split"] == "test", "sequence"].dropna())

    # Double-mutant protein overlap specifically
    dm = df[df["record_type"] == "double_mutants"]
    dm_train_p = set(dm.loc[dm["split"] == "train", "protein_id"].dropna())
    dm_test_p = set(dm.loc[dm["split"] == "test", "protein_id"].dropna())

    # Per-protein row counts (double mutants)
    dm_per_protein = dm.groupby("protein_id").size()
    protein_row_hist = hist_counts(dm_per_protein, bins=[0, 100, 500, 1000, 2000, 4000, 8000, 20000])

    # ΔG histograms (double mutants + external)
    dg_bins = list(np.arange(-2, 8.5, 0.5))
    dm_dg_hist = hist_counts(dm["delta_g"], bins=dg_bins)
    ext = df[df["record_type"] == "external_delta_g"]
    ext_dg_hist = hist_counts(ext["delta_g"], bins=dg_bins)

    # Sequence length (only external_delta_g has sequences)
    seq_len_hist = hist_counts(ext["sequence_length"], bins=list(range(20, 85, 5)))

    # DSSP on site_features + gemme
    dssp = df["dssp"].dropna().astype(str).str.strip()
    dssp_counts = dssp.value_counts().to_dict()

    # AA pair frequencies for double mutants
    aa1 = dm["aa1"].dropna().astype(str).str.upper()
    aa2 = dm["aa2"].dropna().astype(str).str.upper()
    aa_order = list("ACDEFGHIKLMNPQRSTVWY")
    aa1_counts = aa1.value_counts().reindex(aa_order, fill_value=0).to_dict()
    aa2_counts = aa2.value_counts().reindex(aa_order, fill_value=0).to_dict()

    # wt_aa on site-level tables
    sites = df[df["record_type"].isin(["site_features", "natural_sites", "gemme_natural_sites"])]
    wt_counts = (
        sites["wt_aa"].dropna().astype(str).str.upper().value_counts().reindex(aa_order, fill_value=0).to_dict()
    )

    # score_value vs delta_g identity
    both = df.dropna(subset=["delta_g", "score_value"])
    score_eq_delta = float((both["delta_g"] == both["score_value"]).mean()) if len(both) else None
    target_eq_delta = float((df.dropna(subset=["delta_g", "target"])["delta_g"] == df.dropna(subset=["delta_g", "target"])["target"]).mean()) if df["target"].notna().any() else None

    # reconstructed vs measured ΔG on double mutants
    recon = dm.dropna(subset=["delta_g", "reconstructed_delta_g"])
    if len(recon):
        corr = float(np.corrcoef(recon["delta_g"], recon["reconstructed_delta_g"])[0, 1])
        mae = float((recon["delta_g"] - recon["reconstructed_delta_g"]).abs().mean())
    else:
        corr = mae = None

    # Split bucket distribution (should be ~10% per bucket; 0 = test)
    bucket_counts = df["split_bucket"].value_counts().sort_index().to_dict()
    bucket_counts = {str(int(k)): int(v) for k, v in bucket_counts.items()}

    # Unique proteins per type
    proteins_by_type = {
        rec_type: int(group["protein_id"].nunique())
        for rec_type, group in df.groupby("record_type")
    }

    # Top proteins by row count
    top_proteins = (
        dm.groupby("protein_id").size().sort_values(ascending=False).head(12)
    )
    top_proteins = [{"protein": str(k), "rows": int(v)} for k, v in top_proteins.items()]

    # Duplicate record_ids
    dup_ids = int(df["record_id"].duplicated().sum())

    # Coverage of a usable ML row: double mutant with protein + ΔG
    usable_dm = dm.dropna(subset=["protein_id", "delta_g"])
    usable_ext = ext.dropna(subset=["sequence", "delta_g"])

    # ΔG sign / class balance for a typical stability classifier
    def class_balance(series: pd.Series) -> dict:
        values = pd.to_numeric(series, errors="coerce").dropna()
        return {
            "n": int(values.size),
            "lt_0": int((values < 0).sum()),
            "eq_0": int((values == 0).sum()),
            "gt_0": int((values > 0).sum()),
            "gt_1": int((values > 1).sum()),
            "lt_neg1": int((values < -1).sum()),
        }

    # Per-type numeric summaries
    type_summaries = {}
    for rec_type, group in df.groupby("record_type"):
        type_summaries[rec_type] = {
            "n": int(len(group)),
            "n_train": int((group["split"] == "train").sum()),
            "n_test": int((group["split"] == "test").sum()),
            "n_proteins": int(group["protein_id"].nunique()),
            "delta_g": describe(group["delta_g"]),
            "score_value": describe(group["score_value"]),
            "sequence_length": describe(group["sequence_length"]),
            "has_sequence_pct": float(group["sequence"].notna().mean()),
            "has_delta_g_pct": float(group["delta_g"].notna().mean()),
        }

    summary = {
        "source": "LiteFold/MegaScale-Tsuboyama2023",
        "citation": "Tsuboyama et al., Nature 620:434–444 (2023)",
        "license": "CC-BY-4.0",
        "n_rows": n,
        "n_columns_full": 204,
        "splits": {k: int(v) for k, v in split_counts.items()},
        "split_strategy": "sha256(record_id) % 10; bucket 0 = test, 1-9 = train",
        "record_types": {k: int(v) for k, v in type_counts.items()},
        "figures": {k: int(v) for k, v in figure_counts.items()},
        "type_by_split": {
            split: {k: int(v) for k, v in counts.items()}
            for split, counts in type_by_split.items()
        },
        "n_proteins": int(df["protein_id"].nunique()),
        "proteins_by_type": proteins_by_type,
        "duplicate_record_ids": dup_ids,
        "missing_overall": missing_overall,
        "missing_by_type": missing_by_type,
        "leakage": {
            "train_proteins": len(train_proteins),
            "test_proteins": len(test_proteins),
            "protein_overlap": len(protein_overlap),
            "protein_overlap_pct_of_test": round(len(protein_overlap) / max(len(test_proteins), 1), 4),
            "train_sequences": len(train_seq),
            "test_sequences": len(test_seq),
            "sequence_overlap": len(train_seq & test_seq),
            "double_mutant_protein_overlap": len(dm_train_p & dm_test_p),
            "double_mutant_train_proteins": len(dm_train_p),
            "double_mutant_test_proteins": len(dm_test_p),
        },
        "usable": {
            "double_mutants_with_dg": int(len(usable_dm)),
            "external_with_sequence_and_dg": int(len(usable_ext)),
        },
        "delta_g_double_mutants": describe(dm["delta_g"]),
        "delta_g_external": describe(ext["delta_g"]),
        "reconstructed_vs_measured": {"pearson_r": corr, "mae": mae, "n": int(len(recon))},
        "score_equals_delta_g": score_eq_delta,
        "target_equals_delta_g": target_eq_delta,
        "class_balance_double_mutants": class_balance(dm["delta_g"]),
        "class_balance_external": class_balance(ext["delta_g"]),
        "hist_delta_g_double": dm_dg_hist,
        "hist_delta_g_external": ext_dg_hist,
        "hist_seq_len_external": seq_len_hist,
        "hist_rows_per_protein_dm": protein_row_hist,
        "dssp_counts": {k: int(v) for k, v in dssp_counts.items()},
        "aa1_counts": {k: int(v) for k, v in aa1_counts.items()},
        "aa2_counts": {k: int(v) for k, v in aa2_counts.items()},
        "wt_aa_site_counts": {k: int(v) for k, v in wt_counts.items()},
        "split_buckets": bucket_counts,
        "top_proteins_double_mutants": top_proteins,
        "type_summaries": type_summaries,
        "dm_per_protein": {
            "n_proteins": int(dm_per_protein.size),
            "median_rows": float(dm_per_protein.median()),
            "mean_rows": float(dm_per_protein.mean()),
            "max_rows": int(dm_per_protein.max()),
            "min_rows": int(dm_per_protein.min()),
        },
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in [
        "n_rows", "splits", "record_types", "n_proteins", "leakage",
        "usable", "delta_g_double_mutants", "reconstructed_vs_measured",
        "score_equals_delta_g", "class_balance_double_mutants",
        "dm_per_protein",
    ]}, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
