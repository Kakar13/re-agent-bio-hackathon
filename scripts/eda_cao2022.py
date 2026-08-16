#!/usr/bin/env python3
"""EDA for Cao et al. 2022 IPD experimental affinities. Writes JSON for the canvas."""

from __future__ import annotations

import json
import math
import tarfile
from io import TextIOWrapper
from pathlib import Path

import numpy as np
import pandas as pd

ARCHIVE = Path("data/raw/cao2022/experimental_data_and_analysis.tar.gz")
OUT = Path("results/cao2022/eda_summary.json")

AFFINITY_PREFIX = "supplemental_files/ngs_analysis/affinities/"
SKIP_SUBDIR = "affinities_no_doubly_transformed_thresold"
MAINLINE = (
    "supplemental_files/ngs_analysis/ssm_heatmaps/"
    "mainline_binders_for_manuscript/mainline_affinity.sc"
)

BOOL_COLS = {
    "low_conf",
    "avid_doesnt_agree",
    "binder_4000_nm",
    "binder_400_nm",
    "binder_800_nm",
}


def _to_float(tok: str) -> float:
    if tok in {"inf", "Inf", "INF", "+inf", "+Inf"}:
        return math.inf
    if tok in {"-inf", "-Inf", "-INF"}:
        return -math.inf
    if tok in {"nan", "NaN", "NA", "None"}:
        return math.nan
    return float(tok)


def _to_bool(tok: str) -> bool:
    return tok.lower() in {"true", "1", "t", "yes"}


def parse_sc(text: str) -> pd.DataFrame:
    """Header-driven parse. Last two header tokens are target + description."""
    rows: list[dict] = []
    header: list[str] | None = None
    n_core = 0
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if header is None:
            if "kd_lb" not in parts or "description" not in parts:
                continue
            header = parts
            n_core = len(header) - 2
            continue
        if len(parts) < n_core + 2:
            continue
        rec: dict = {}
        ok = True
        for i, name in enumerate(header[:n_core]):
            tok = parts[i]
            try:
                rec[name] = (
            _to_bool(tok)
            if name in BOOL_COLS or name.startswith("binder_")
            else _to_float(tok)
        )
            except ValueError:
                ok = False
                break
        if not ok:
            continue
        rec["target"] = parts[n_core]
        rec["description"] = " ".join(parts[n_core + 1 :])
        rows.append(rec)
    return pd.DataFrame(rows)


def midpoint_kd(lb: pd.Series, ub: pd.Series) -> pd.Series:
    """Geometric midpoint when both bounds finite; else the finite side."""
    lo = lb.to_numpy(dtype=float)
    hi = ub.to_numpy(dtype=float)
    out = np.full(len(lo), np.nan)
    both = np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi > 0)
    out[both] = np.sqrt(lo[both] * hi[both])
    only_lo = np.isfinite(lo) & ~np.isfinite(hi) & (lo > 0)
    out[only_lo] = lo[only_lo]
    only_hi = ~np.isfinite(lo) & np.isfinite(hi) & (hi > 0)
    out[only_hi] = hi[only_hi]
    return pd.Series(out)


def describe(s: pd.Series) -> dict:
    v = pd.to_numeric(s, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if v.empty:
        return {"n": 0}
    return {
        "n": int(v.size),
        "mean": float(v.mean()),
        "std": float(v.std()) if v.size > 1 else 0.0,
        "min": float(v.min()),
        "p05": float(v.quantile(0.05)),
        "p25": float(v.quantile(0.25)),
        "median": float(v.median()),
        "p75": float(v.quantile(0.75)),
        "p95": float(v.quantile(0.95)),
        "max": float(v.max()),
    }


def hist_log10(s: pd.Series, bins: list[float]) -> list[dict]:
    v = pd.to_numeric(s, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    v = v[v > 0]
    if v.empty:
        return []
    counts, edges = np.histogram(np.log10(v), bins=bins)
    return [
        {"label": f"{edges[i]:g}–{edges[i + 1]:g}", "count": int(counts[i])}
        for i in range(len(counts))
    ]


def is_ssm(name: str, desc: pd.Series) -> bool:
    if name.endswith("_ssm") or name.endswith("_ssm_2"):
        return True
    return bool(desc.astype(str).str.contains(r"__\d+__[A-Z]", regex=True).mean() > 0.5)


def main() -> None:
    if not ARCHIVE.exists():
        raise FileNotFoundError(ARCHIVE)

    frames: list[pd.DataFrame] = []
    mainline = pd.DataFrame()
    archive_files = 0
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        for m in tar.getmembers():
            if not m.isfile():
                continue
            archive_files += 1
            name = m.name
            is_main = name == MAINLINE
            in_aff = name.startswith(AFFINITY_PREFIX) and name.endswith(".sc")
            if in_aff and SKIP_SUBDIR in name:
                continue
            if not (in_aff or is_main):
                continue
            extracted = tar.extractfile(m)
            if extracted is None:
                continue
            text = TextIOWrapper(extracted, encoding="utf-8", errors="replace").read()
            frame = parse_sc(text)
            frame["file"] = Path(name).name
            frame["kind"] = "mainline" if is_main else "affinity"
            if is_main:
                mainline = frame
            else:
                frames.append(frame)
            print(f"  {Path(name).name:28} {len(frame):7,} rows")

    df = pd.concat(frames, ignore_index=True)
    df["kd_mid"] = midpoint_kd(df["kd_lb"], df["kd_ub"])
    df["lb_finite"] = np.isfinite(df["kd_lb"].to_numpy(dtype=float))
    df["ub_finite"] = np.isfinite(df["kd_ub"].to_numpy(dtype=float))
    df["both_finite"] = df["lb_finite"] & df["ub_finite"]
    df["unbound"] = ~df["lb_finite"] & ~df["ub_finite"]
    df["lower_only"] = df["lb_finite"] & ~df["ub_finite"]
    df["is_ssm"] = df.apply(
        lambda r: is_ssm(str(r["target"]), pd.Series([r["description"]])), axis=1
    )
    # faster ssm flag
    df["is_ssm"] = df["target"].astype(str).str.contains(r"_ssm", regex=True) | df[
        "file"
    ].str.contains(r"_ssm")

    targets = []
    for tgt, g in df.groupby("target", sort=False):
        mid = g["kd_mid"]
        targets.append(
            {
                "target": str(tgt),
                "file": str(g["file"].iloc[0]),
                "n": int(len(g)),
                "n_unique_desc": int(g["description"].nunique()),
                "is_ssm": bool(g["is_ssm"].mean() > 0.5),
                "both_finite": int(g["both_finite"].sum()),
                "lower_only": int(g["lower_only"].sum()),
                "unbound": int(g["unbound"].sum()),
                "low_conf": int(g["low_conf"].sum()) if "low_conf" in g else 0,
                "kd_mid": describe(mid),
                "lt_10_nM": int((mid < 10).sum()),
                "lt_100_nM": int((mid < 100).sum()),
                "lt_1_uM": int((mid < 1_000).sum()),
                "lt_10_uM": int((mid < 10_000).sum()),
            }
        )
    targets.sort(key=lambda r: r["n"], reverse=True)

    mid = df["kd_mid"]
    tightness = {
        "lt_1_nM": int((mid < 1).sum()),
        "lt_10_nM": int((mid < 10).sum()),
        "lt_100_nM": int((mid < 100).sum()),
        "lt_1_uM": int((mid < 1_000).sum()),
        "lt_10_uM": int((mid < 10_000).sum()),
        "ge_10_uM": int((mid >= 10_000).sum()),
        "no_midpoint": int(mid.isna().sum()),
    }

    samples = []
    finite = df.loc[df["both_finite"]].copy()
    if not finite.empty:
        tight = finite.loc[finite["kd_mid"].idxmin()]
        samples.append(
            {
                "tag": "tightest both-bounds",
                "target": str(tight["target"]),
                "description": str(tight["description"])[:72],
                "kd_lb": float(tight["kd_lb"]),
                "kd_ub": float(tight["kd_ub"]),
                "kd_mid": float(tight["kd_mid"]),
                "low_conf": bool(tight["low_conf"]),
            }
        )
    lo_only = df.loc[df["lower_only"] & np.isfinite(df["kd_lb"].to_numpy(dtype=float))]
    if not lo_only.empty:
        row = lo_only.loc[lo_only["kd_lb"].idxmin()]
        samples.append(
            {
                "tag": "tightest lower-only (ub=inf)",
                "target": str(row["target"]),
                "description": str(row["description"])[:72],
                "kd_lb": float(row["kd_lb"]),
                "kd_ub": None,
                "kd_mid": float(row["kd_mid"]) if pd.notna(row["kd_mid"]) else None,
                "low_conf": bool(row["low_conf"]),
            }
        )
    unbound = df.loc[df["unbound"]]
    if not unbound.empty:
        row = unbound.iloc[0]
        samples.append(
            {
                "tag": "fully censored (inf/inf)",
                "target": str(row["target"]),
                "description": str(row["description"])[:72],
                "kd_lb": None,
                "kd_ub": None,
                "kd_mid": None,
                "low_conf": bool(row["low_conf"]),
            }
        )
    ssm = df.loc[df["is_ssm"] & df["both_finite"]]
    if not ssm.empty:
        row = ssm.iloc[0]
        samples.append(
            {
                "tag": "SSM point mutant",
                "target": str(row["target"]),
                "description": str(row["description"])[:72],
                "kd_lb": float(row["kd_lb"]),
                "kd_ub": float(row["kd_ub"]),
                "kd_mid": float(row["kd_mid"]),
                "low_conf": bool(row["low_conf"]),
            }
        )

    mainline_stats = {}
    if not mainline.empty:
        mainline["kd_mid"] = midpoint_kd(mainline["kd_lb"], mainline["kd_ub"])
        mainline_stats = {
            "n": int(len(mainline)),
            "n_unique_desc": int(mainline["description"].nunique()),
            "targets": sorted(mainline["target"].dropna().unique().tolist()),
            "kd_mid": describe(mainline["kd_mid"]),
            "both_finite": int(
                (
                    np.isfinite(mainline["kd_lb"].to_numpy(dtype=float))
                    & np.isfinite(mainline["kd_ub"].to_numpy(dtype=float))
                ).sum()
            ),
        }

    designed = df.loc[~df["is_ssm"]]
    ssm_df = df.loc[df["is_ssm"]]

    summary = {
        "paper": "Cao et al., Nature 2022 — Robust de novo design of protein binding proteins",
        "archive": str(ARCHIVE),
        "archive_bytes": ARCHIVE.stat().st_size,
        "archive_file_count": archive_files,
        "n_affinity_rows": int(len(df)),
        "n_unique_designs": int(df["description"].nunique()),
        "n_targets": int(df["target"].nunique()),
        "n_designed_rows": int(len(designed)),
        "n_ssm_rows": int(len(ssm_df)),
        "n_unique_designed": int(designed["description"].nunique()),
        "n_unique_ssm": int(ssm_df["description"].nunique()),
        "censoring": {
            "both_finite": int(df["both_finite"].sum()),
            "lower_only": int(df["lower_only"].sum()),
            "unbound": int(df["unbound"].sum()),
            "low_conf": int(df["low_conf"].sum()),
        },
        "kd_mid": describe(mid),
        "hist_log10_kd_mid": hist_log10(mid, list(np.arange(-1.0, 7.0, 0.5))),
        "tightness": tightness,
        "has_sequence_column": False,
        "targets": targets,
        "mainline": mainline_stats,
        "samples": samples,
        "units": (
            "kd_lb / kd_ub are FACS/NGS enrichment-derived nM estimates, "
            "censored at the sort dynamic range (inf). kd_mid is the geometric "
            "mean when both bounds are finite, else the finite side. "
            "Not SPR/BLI. description is a design name, not a sequence."
        ),
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {OUT}")
    print(
        f"rows={len(df):,}  designs={summary['n_unique_designs']:,}  "
        f"targets={summary['n_targets']}  both_finite={summary['censoring']['both_finite']:,}  "
        f"unbound={summary['censoring']['unbound']:,}"
    )


if __name__ == "__main__":
    main()
