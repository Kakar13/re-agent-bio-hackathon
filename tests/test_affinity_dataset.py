"""Smoke test: build a 1000-row AlphaSeq table and assert the export schema."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "build_affinity_dataset.py"
OUT = REPO / "data" / "processed" / "affinity" / "alphaseq_smoke.csv"

EXPECTED = [
    "binder_sequence",
    "target_sequence",
    "binder_desc",
    "target_desc",
    "log10_kd",
    "log10_kd_lower",
    "log10_kd_upper",
    "kd_nm",
    "kd_molar",
    "censored",
    "source",
    "assay_type",
]


def test_alphaseq_smoke_schema() -> None:
    cmd = [
        sys.executable,
        str(SCRIPT),
        "alphaseq",
        "--max-rows",
        "1000",
        "--out",
        str(OUT),
        "--cache",
        str(REPO / "data" / "raw"),
    ]
    subprocess.run(cmd, check=True, cwd=REPO)
    df = pd.read_csv(OUT)
    assert list(df.columns)[: len(EXPECTED)] == EXPECTED, list(df.columns)
    assert len(df) == 1000
    assert df["binder_sequence"].notna().any()
    assert df["target_sequence"].notna().any()
    assert df["source"].nunique() >= 1
    assert set(df["assay_type"].dropna()) == {"alphaseq_yeast_mating"}
    # Default keeps censored rows; they must be marked, not dropped.
    assert "censored" in df.columns
    # Units: log10_kd is log10(nM). 1 nM → kd_nm=1, kd_molar=1e-9.
    observed = df.loc[~df["censored"].astype(bool) & df["log10_kd"].notna()]
    if not observed.empty:
        row = observed.iloc[0]
        expected_nm = 10.0 ** float(row["log10_kd"])
        assert abs(float(row["kd_nm"]) - expected_nm) / expected_nm < 1e-6
        assert abs(float(row["kd_molar"]) - expected_nm * 1e-9) / (expected_nm * 1e-9) < 1e-6


def test_binders_only_reports_drop() -> None:
    out = REPO / "data" / "processed" / "affinity" / "alphaseq_smoke_binders.csv"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "alphaseq",
            "--max-rows",
            "1000",
            "--binders-only",
            "--out",
            str(out),
            "--cache",
            str(REPO / "data" / "raw"),
        ],
        check=True,
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    df = pd.read_csv(out)
    assert not bool(df["censored"].astype(bool).any())
    assert "--binders-only" in proc.stdout or "--binders-only" in proc.stderr


if __name__ == "__main__":
    test_alphaseq_smoke_schema()
    test_binders_only_reports_drop()
    print("smoke ok")
