#!/usr/bin/env python3
"""EDA for the processed AlphaSeq affinity CSV. Writes JSON for the canvas."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

CSV = Path("data/processed/affinity/alphaseq_affinity.csv")
OUT = Path("results/affinity/alphaseq_csv_eda.json")

CAMPAIGNS = {
    "YM_0005": "CoV Ab panel × RBD mutants",
    "YM_0549": "VHH72 Iter0 × CoV2-RBD",
    "YM_1068": "VHH72 Iter1 × CoV2-RBD",
    "YM_0693": "PP489 Iter0 × TIGIT",
    "YM_0988": "PP489 Iter1 × TIGIT",
    "YM_0852": "Pembro-scFv Iter0 × PD-1",
    "YM_0985": "Pembro-scFv Iter1 × PD-1",
    "YM_0989": "Trastuzumab-scFv CDR3 × HER-2",
    "YM_0990": "Trastuzumab-scFv CDR+FW × HER-2",
}

PAIRS = [
    ("VHH72 / CoV2-RBD", "YM_0549", "YM_1068"),
    ("PP489 / TIGIT", "YM_0693", "YM_0988"),
    ("Pembro / PD-1", "YM_0852", "YM_0985"),
    ("Trastuzumab / HER-2", "YM_0989", "YM_0990"),
]


def describe(s: pd.Series) -> dict:
    v = pd.to_numeric(s, errors="coerce").dropna()
    if v.empty:
        return {"n": 0}
    return {
        "n": int(v.size),
        "mean": float(v.mean()),
        "std": float(v.std()),
        "min": float(v.min()),
        "p05": float(v.quantile(0.05)),
        "p25": float(v.quantile(0.25)),
        "median": float(v.median()),
        "p75": float(v.quantile(0.75)),
        "p95": float(v.quantile(0.95)),
        "max": float(v.max()),
    }


def hist(s: pd.Series, bins: list[float]) -> list[dict]:
    v = pd.to_numeric(s, errors="coerce").dropna()
    counts, edges = np.histogram(v, bins=bins)
    return [
        {"label": f"{edges[i]:g}–{edges[i + 1]:g}", "count": int(counts[i])}
        for i in range(len(counts))
    ]


def nonempty(s: pd.Series) -> pd.Series:
    return s.notna() & (s.astype(str).str.len() > 0) & (s.astype(str) != "nan")


def main() -> None:
    print(f"reading {CSV} …")
    df = pd.read_csv(CSV)
    print(f"loaded {len(df):,} rows, {list(df.columns)}")

    binder_ok = nonempty(df["binder_sequence"])
    target_ok = nonempty(df["target_sequence"])
    kd = pd.to_numeric(df["log10_kd"], errors="coerce")
    observed = kd.notna()
    censored = df["censored"].astype(str).str.lower().isin(("true", "1")) | kd.isna()

    kd_nm = pd.to_numeric(df["kd_nm"], errors="coerce")
    both_seq = binder_ok & target_ok
    train_slice = both_seq & observed

    tightness = {
        "lt_1_nM": int((kd_nm < 1).sum()),
        "lt_10_nM": int((kd_nm < 10).sum()),
        "lt_100_nM": int((kd_nm < 100).sum()),
        "lt_1_uM": int((kd_nm < 1_000).sum()),
        "lt_10_uM": int((kd_nm < 10_000).sum()),
        "ge_10_uM": int((kd_nm >= 10_000).sum()),
    }

    pair_key = df["binder_sequence"].fillna("").astype(str) + "||" + df[
        "target_sequence"
    ].fillna("").astype(str)
    pair_counts = pair_key.value_counts()
    n_dup_extra = int((pair_counts - 1).clip(lower=0).sum())

    desc_neg = (
        df["binder_desc"].astype(str).str.contains("Neg", case=False, na=False)
        | df["target_desc"].astype(str).str.contains("Neg", case=False, na=False)
    )

    sources = []
    binder_sets: dict[str, set[str]] = {}
    target_sets: dict[str, set[str]] = {}
    for src, label in CAMPAIGNS.items():
        g = df.loc[df["source"] == src]
        b_ok = nonempty(g["binder_sequence"])
        t_ok = nonempty(g["target_sequence"])
        g_kd = pd.to_numeric(g["log10_kd"], errors="coerce")
        g_nm = pd.to_numeric(g["kd_nm"], errors="coerce")
        binders = set(g.loc[b_ok, "binder_sequence"].astype(str))
        targets = set(g.loc[t_ok, "target_sequence"].astype(str))
        binder_sets[src] = binders
        target_sets[src] = targets
        sources.append(
            {
                "source": src,
                "campaign": label,
                "n": int(len(g)),
                "n_observed": int(g_kd.notna().sum()),
                "n_censored": int(g_kd.isna().sum()),
                "censored_pct": float(g_kd.isna().mean()),
                "n_binder_seq": int(len(binders)),
                "n_target_seq": int(len(targets)),
                "both_seq_and_kd": int((b_ok & t_ok & g_kd.notna()).sum()),
                "both_seq_and_kd_pct": float((b_ok & t_ok & g_kd.notna()).mean()),
                "neg_control_rows": int(
                    (
                        g["binder_desc"].astype(str).str.contains("Neg", case=False, na=False)
                        | g["target_desc"].astype(str).str.contains("Neg", case=False, na=False)
                    ).sum()
                ),
                "affinity": describe(g_kd),
                "kd_nm": describe(g_nm),
                "binder_len": describe(g.loc[b_ok, "binder_sequence"].astype(str).str.len()),
                "target_len": describe(g.loc[t_ok, "target_sequence"].astype(str).str.len()),
                "lt_100_nM": int((g_nm < 100).sum()),
                "lt_1_uM": int((g_nm < 1_000).sum()),
            }
        )

    overlap = []
    for name, a, b in PAIRS:
        overlap.append(
            {
                "campaign": name,
                "iter0": a,
                "iter1": b,
                "binder_overlap": len(binder_sets[a] & binder_sets[b]),
                "target_overlap": len(target_sets[a] & target_sets[b]),
                "iter0_binders": len(binder_sets[a]),
                "iter1_binders": len(binder_sets[b]),
            }
        )

    # Sample rows: tightest observed with both seqs, median-ish, censored, per family
    samples = []

    def add_sample(mask: pd.Series, tag: str) -> None:
        hit = df.loc[mask]
        if hit.empty:
            return
        if "log10_kd" in hit and hit["log10_kd"].notna().any() and tag != "censored":
            row = hit.loc[hit["log10_kd"].idxmin()]
        else:
            row = hit.iloc[0]
        samples.append(
            {
                "tag": tag,
                "source": str(row["source"]),
                "binder_desc": str(row["binder_desc"])[:48],
                "target_desc": str(row["target_desc"])[:48],
                "binder_len": int(len(str(row["binder_sequence"]))) if pd.notna(row["binder_sequence"]) else 0,
                "target_len": int(len(str(row["target_sequence"]))) if pd.notna(row["target_sequence"]) else 0,
                "log10_kd": None if pd.isna(row["log10_kd"]) else float(row["log10_kd"]),
                "kd_nm": None if pd.isna(row["kd_nm"]) else float(row["kd_nm"]),
                "censored": bool(pd.isna(row["log10_kd"])),
                "binder_head": str(row["binder_sequence"])[:18] + "…"
                if pd.notna(row["binder_sequence"]) and str(row["binder_sequence"]) not in ("", "nan")
                else "",
            }
        )

    add_sample(train_slice & (kd_nm < 1), "tightest <1 nM")
    add_sample(train_slice & (kd_nm >= 8_000) & (kd_nm <= 9_000), "median-ish ~8 µM")
    add_sample(censored & both_seq, "censored, both seqs")
    add_sample((df["source"] == "YM_1068") & train_slice, "VHH72 Iter1")
    add_sample((df["source"] == "YM_0988") & train_slice, "PP489 Iter1")
    add_sample((df["source"] == "YM_0985") & train_slice, "Pembro Iter1")
    add_sample((df["source"] == "YM_0005") & train_slice, "CoV panel")

    # Unique binder lengths (deduped)
    uniq_b = df.loc[binder_ok, ["binder_sequence"]].drop_duplicates()
    uniq_t = df.loc[target_ok, ["target_sequence"]].drop_duplicates()

    # Bound coverage
    lo = pd.to_numeric(df["log10_kd_lower"], errors="coerce")
    hi = pd.to_numeric(df["log10_kd_upper"], errors="coerce")

    summary = {
        "file": str(CSV),
        "n_rows": int(len(df)),
        "columns": list(df.columns),
        "null_pct": {c: float(df[c].isna().mean()) for c in df.columns},
        "assay_type": sorted(df["assay_type"].dropna().unique().tolist()),
        "n_observed": int(observed.sum()),
        "n_censored": int(censored.sum()),
        "n_binder_seq": int(df.loc[binder_ok, "binder_sequence"].nunique()),
        "n_target_seq": int(df.loc[target_ok, "target_sequence"].nunique()),
        "n_binder_missing": int((~binder_ok).sum()),
        "n_target_missing": int((~target_ok).sum()),
        "n_both_seq": int(both_seq.sum()),
        "n_train_slice": int(train_slice.sum()),
        "n_neg_control_desc": int(desc_neg.sum()),
        "n_unique_pairs": int(len(pair_counts)),
        "n_replicate_extra_rows": n_dup_extra,
        "n_pairs_with_replicates": int((pair_counts > 1).sum()),
        "affinity": describe(kd),
        "kd_nm": describe(kd_nm),
        "hist_log10_kd": hist(kd, list(np.arange(-1.0, 8.0, 0.5))),
        "hist_binder_len_unique": hist(uniq_b["binder_sequence"].astype(str).str.len(), [100, 120, 130, 200, 240, 250, 260, 600]),
        "hist_target_len_unique": hist(uniq_t["target_sequence"].astype(str).str.len(), [100, 120, 130, 190, 210, 230, 250]),
        "tightness": tightness,
        "bounds": {
            "lower_present": int(lo.notna().sum()),
            "upper_present": int(hi.notna().sum()),
            "both_present": int((lo.notna() & hi.notna()).sum()),
        },
        "sources": sources,
        "campaign_overlap": overlap,
        "samples": samples,
        "units": "log10_kd = log10(Kd in nM); kd_nm = 10**log10_kd; kd_molar = kd_nm * 1e-9",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2))
    print(f"wrote {OUT}")
    print(
        f"train_slice={summary['n_train_slice']:,}  "
        f"censored={summary['n_censored']:,}  "
        f"lt_100nM={tightness['lt_100_nM']:,}  "
        f"lt_1uM={tightness['lt_1_uM']:,}"
    )


if __name__ == "__main__":
    main()
