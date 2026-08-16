"""EDA for aalphabio/open-alphaseq. Writes JSON for the canvas."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("data/raw/open_alphaseq")
OUT = Path("results/open_alphaseq/eda_summary.json")

EXPERIMENTS = {
    "YM_0005": "Anti-CoV Ab × CoV2-RBD mutagenesis",
    "YM_0549": "VHH72 Iter0 × CoV2-RBD",
    "YM_1068": "VHH72 Iter1 × CoV2-RBD",
    "YM_0693": "PP489 Iter0 × TIGIT",
    "YM_0988": "PP489 Iter1 × TIGIT",
    "YM_0852": "Pembro-scFv Iter0 × PD-1",
    "YM_0985": "Pembro-scFv Iter1 × PD-1",
    "YM_0989": "Trastuzumab-scFv CDR3 Iter0 × HER-2",
    "YM_0990": "Trastuzumab-scFv CDR+FW Iter0 × HER-2",
}

CAMPAIGNS = {
    "VHH72 / CoV2-RBD": ["YM_0549", "YM_1068"],
    "PP489 / TIGIT": ["YM_0693", "YM_0988"],
    "Pembro / PD-1": ["YM_0852", "YM_0985"],
    "Trastuzumab / HER-2": ["YM_0989", "YM_0990"],
    "CoV Ab panel": ["YM_0005"],
}


def hist_counts(series: pd.Series, bins: list[float]) -> list[dict]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return []
    counts, edges = np.histogram(values, bins=bins)
    return [
        {"label": f"{edges[i]:g}–{edges[i + 1]:g}", "count": int(counts[i])}
        for i in range(len(counts))
    ]


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
    }


def seq_len(series: pd.Series) -> pd.Series:
    return series.dropna().astype(str).str.len()


def main() -> None:
    frames = []
    for ym in EXPERIMENTS:
        path = ROOT / "data" / ym / "data.parquet"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_parquet(path)
        frame["experiment"] = ym
        frames.append(frame)
    df = pd.concat(frames, ignore_index=True)

    for col in [
        "alphaseq_affinity",
        "affinity_lower_bound",
        "affinity_upper_bound",
        "normalized_affinity",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    n = len(df)
    per_exp = []
    for ym, group in df.groupby("experiment"):
        aff = group["alphaseq_affinity"]
        mata_seq = group["mata_sequence"]
        alpha_seq = group["matalpha_sequence"]
        per_exp.append(
            {
                "experiment": ym,
                "description": EXPERIMENTS[ym],
                "n": int(len(group)),
                "n_affinity": int(aff.notna().sum()),
                "affinity_missing_pct": float(aff.isna().mean()),
                "n_mata_seq": int(mata_seq.nunique(dropna=True)),
                "n_matalpha_seq": int(alpha_seq.nunique(dropna=True)),
                "n_mata_desc": int(group["mata_description"].nunique(dropna=True)),
                "n_matalpha_desc": int(group["matalpha_description"].nunique(dropna=True)),
                "both_sequences_pct": float(
                    (mata_seq.notna() & (mata_seq.astype(str).str.len() > 0)
                     & alpha_seq.notna() & (alpha_seq.astype(str).str.len() > 0)).mean()
                ),
                "above_background_pct": float(group["above_background"].mean())
                if "above_background" in group
                else None,
                "sufficient_replicate_pct": float(group["sufficient_replicate_observations"].mean())
                if "sufficient_replicate_observations" in group
                else None,
                "affinity": describe(aff),
                "mata_len": describe(seq_len(mata_seq)),
                "matalpha_len": describe(seq_len(alpha_seq)),
            }
        )

    usable = df.dropna(subset=["alphaseq_affinity"])
    usable = usable[
        usable["mata_sequence"].notna()
        & (usable["mata_sequence"].astype(str).str.len() > 0)
        & usable["matalpha_sequence"].notna()
        & (usable["matalpha_sequence"].astype(str).str.len() > 0)
    ]

    # Campaign sequence overlap (Iter0 vs Iter1 leakage)
    campaign_overlap = []
    for name, yms in CAMPAIGNS.items():
        if len(yms) < 2:
            continue
        a = df[df["experiment"] == yms[0]]
        b = df[df["experiment"] == yms[1]]
        mata_overlap = set(a["mata_sequence"].dropna()) & set(b["mata_sequence"].dropna())
        alpha_overlap = set(a["matalpha_sequence"].dropna()) & set(b["matalpha_sequence"].dropna())
        campaign_overlap.append(
            {
                "campaign": name,
                "iter0": yms[0],
                "iter1": yms[1],
                "mata_overlap": len(mata_overlap),
                "matalpha_overlap": len(alpha_overlap),
                "iter0_mata": int(a["mata_sequence"].nunique(dropna=True)),
                "iter1_mata": int(b["mata_sequence"].nunique(dropna=True)),
            }
        )

    # Duplicate pair rate
    pair_cols = ["experiment", "mata_sequence", "matalpha_sequence"]
    dup_pairs = int(df.duplicated(subset=pair_cols).sum())

    # Negative-control-ish descriptions
    neg_mask = (
        df["mata_description"].astype(str).str.contains("Neg", case=False, na=False)
        | df["matalpha_description"].astype(str).str.contains("Neg", case=False, na=False)
    )

    summary = {
        "source": "aalphabio/open-alphaseq",
        "license": "open (academic + commercial; see repo LICENSE)",
        "n_rows": n,
        "n_experiments": len(EXPERIMENTS),
        "n_with_affinity": int(df["alphaseq_affinity"].notna().sum()),
        "affinity_missing_pct": float(df["alphaseq_affinity"].isna().mean()),
        "usable_both_seq_and_affinity": int(len(usable)),
        "n_unique_mata_seq": int(df["mata_sequence"].nunique(dropna=True)),
        "n_unique_matalpha_seq": int(df["matalpha_sequence"].nunique(dropna=True)),
        "duplicate_pairs_within_experiment": dup_pairs,
        "negative_control_rows": int(neg_mask.sum()),
        "affinity_all": describe(df["alphaseq_affinity"]),
        "hist_affinity": hist_counts(df["alphaseq_affinity"], bins=list(np.arange(-2, 8.5, 0.5))),
        "hist_mata_len": hist_counts(seq_len(df["mata_sequence"]), bins=list(range(0, 280, 20))),
        "hist_matalpha_len": hist_counts(
            seq_len(df["matalpha_sequence"]), bins=list(range(0, 280, 20))
        ),
        "above_background": {
            "true": int(df["above_background"].sum()) if "above_background" in df else None,
            "pct": float(df["above_background"].mean()) if "above_background" in df else None,
        },
        "sufficient_replicate": {
            "true": int(df["sufficient_replicate_observations"].sum())
            if "sufficient_replicate_observations" in df
            else None,
            "pct": float(df["sufficient_replicate_observations"].mean())
            if "sufficient_replicate_observations" in df
            else None,
        },
        "experiments": per_exp,
        "campaign_overlap": campaign_overlap,
        "label": "log10(estimated Kd in nM); lower = tighter binding; 0 = 1 nM, 3 = 1 µM",
        "split_advice": "Hold out a whole experiment (or campaign). Do not mix Iter0 train with Iter1 test of the same parent.",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2))
    print(json.dumps({
        "n_rows": summary["n_rows"],
        "n_with_affinity": summary["n_with_affinity"],
        "usable_both_seq_and_affinity": summary["usable_both_seq_and_affinity"],
        "affinity_missing_pct": summary["affinity_missing_pct"],
        "n_unique_mata_seq": summary["n_unique_mata_seq"],
        "n_unique_matalpha_seq": summary["n_unique_matalpha_seq"],
        "affinity_all": summary["affinity_all"],
        "campaign_overlap": summary["campaign_overlap"],
        "per_exp_n": {e["experiment"]: e["n"] for e in per_exp},
    }, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
