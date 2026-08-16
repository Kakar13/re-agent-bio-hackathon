#!/usr/bin/env python3
"""Build a flat sequence → binding-affinity table from public sources.

Subcommands:
  alphaseq   aalphabio/open-alphaseq (primary)
  cao        Cao et al. 2022 IPD dump (inspect first; parser after headers)

log10_kd in AlphaSeq is log10(estimated Kd in nM). kd_nm = 10**log10_kd;
kd_molar = kd_nm * 1e-9. NaN affinities are censored hard negatives — kept
by default.
"""

from __future__ import annotations

import argparse
import sys
import tarfile
import warnings
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = REPO / "data" / "raw"
DEFAULT_OUT_DIR = REPO / "data" / "processed" / "affinity"

ALPHASEQ_REPO = "aalphabio/open-alphaseq"
CAO_BASE = (
    "https://files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/supplemental_files/"
)
CAO_ARCHIVE = "experimental_data_and_analysis.tar.gz"

# Published names plus the aliases actually shipped on Hugging Face.
BINDER_SEQ = ("mata_sequence",)
TARGET_SEQ = ("matalpha_sequence",)
BINDER_DESC = ("mata_description",)
TARGET_DESC = ("matalpha_description",)
LOG10_KD = ("alphaseq_affinity",)
LOG10_KD_LO = ("alphaseq_affinity_lower_bound", "affinity_lower_bound")
LOG10_KD_HI = ("alphaseq_affinity_upper_bound", "affinity_upper_bound")

EXPORT_COLUMNS = [
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


def _pick(columns: list[str], aliases: tuple[str, ...]) -> str | None:
    present = {c.lower(): c for c in columns}
    for alias in aliases:
        if alias.lower() in present:
            return present[alias.lower()]
    return None


def _warn_missing(config: str, columns: list[str], needed: dict[str, str | None]) -> bool:
    missing = [name for name, col in needed.items() if col is None]
    if not missing:
        return False
    print(
        f"WARNING: config {config!r} is missing expected columns {missing}.\n"
        f"  Actual columns: {columns}\n"
        f"  Skipping this config rather than writing a wrong file.",
        file=sys.stderr,
    )
    return True


def _distribution_summary(df: pd.DataFrame, title: str) -> None:
    kd = pd.to_numeric(df["log10_kd"], errors="coerce")
    observed = kd.dropna()
    print(f"\n=== {title} ===")
    print(f"rows                {len(df):,}")
    print(f"censored (NaN Kd)   {int(df['censored'].sum()):,}  ({df['censored'].mean():.1%})")
    print(f"unique binders      {df['binder_sequence'].nunique(dropna=True):,}")
    print(f"unique targets      {df['target_sequence'].nunique(dropna=True):,}")
    if observed.empty:
        print("log10_kd            no observed values")
    else:
        print(
            "log10_kd (nM)       "
            f"n={len(observed):,}  min={observed.min():.3f}  "
            f"median={observed.median():.3f}  max={observed.max():.3f}"
        )
    print(f"sources             {sorted(df['source'].dropna().unique().tolist())}")
    print(f"assay_type          {sorted(df['assay_type'].dropna().unique().tolist())}")


def _normalize_alphaseq(frame: pd.DataFrame, source: str) -> pd.DataFrame | None:
    cols = list(frame.columns)
    mapping = {
        "binder_sequence": _pick(cols, BINDER_SEQ),
        "target_sequence": _pick(cols, TARGET_SEQ),
        "binder_desc": _pick(cols, BINDER_DESC),
        "target_desc": _pick(cols, TARGET_DESC),
        "log10_kd": _pick(cols, LOG10_KD),
        "log10_kd_lower": _pick(cols, LOG10_KD_LO),
        "log10_kd_upper": _pick(cols, LOG10_KD_HI),
    }
    required = {k: mapping[k] for k in ("binder_sequence", "target_sequence", "log10_kd")}
    if _warn_missing(source, cols, required):
        return None

    out = pd.DataFrame(
        {
            "binder_sequence": frame[mapping["binder_sequence"]],
            "target_sequence": frame[mapping["target_sequence"]],
            "binder_desc": frame[mapping["binder_desc"]] if mapping["binder_desc"] else "",
            "target_desc": frame[mapping["target_desc"]] if mapping["target_desc"] else "",
            "log10_kd": pd.to_numeric(frame[mapping["log10_kd"]], errors="coerce"),
            "log10_kd_lower": (
                pd.to_numeric(frame[mapping["log10_kd_lower"]], errors="coerce")
                if mapping["log10_kd_lower"]
                else pd.NA
            ),
            "log10_kd_upper": (
                pd.to_numeric(frame[mapping["log10_kd_upper"]], errors="coerce")
                if mapping["log10_kd_upper"]
                else pd.NA
            ),
        }
    )
    out["kd_nm"] = 10.0 ** out["log10_kd"]
    out["kd_molar"] = out["kd_nm"] * 1e-9
    out["censored"] = out["log10_kd"].isna()
    out["source"] = source
    out["assay_type"] = "alphaseq_yeast_mating"
    return out[EXPORT_COLUMNS]


def _list_alphaseq_configs(cache: Path) -> list[str]:
    try:
        from datasets import get_dataset_config_names

        names = get_dataset_config_names(ALPHASEQ_REPO)
        if names:
            return list(names)
    except Exception as exc:
        print(f"WARNING: could not list Hugging Face configs ({exc}).", file=sys.stderr)

    local = cache / "open_alphaseq" / "data"
    if local.is_dir():
        found = sorted(p.name for p in local.iterdir() if (p / "data.parquet").exists())
        if found:
            print(f"Using local AlphaSeq configs from {local}: {found}", file=sys.stderr)
            return found
    raise RuntimeError(
        "Could not enumerate AlphaSeq configs from Hugging Face or "
        f"{cache / 'open_alphaseq' / 'data'}."
    )


def _load_alphaseq_config(name: str, cache: Path) -> pd.DataFrame:
    local = cache / "open_alphaseq" / "data" / name / "data.parquet"
    if local.exists():
        return pd.read_parquet(local)

    from datasets import load_dataset

    ds = load_dataset(ALPHASEQ_REPO, name, split="all", cache_dir=str(cache / "hf"))
    return ds.to_pandas()


def _collapse_replicates(df: pd.DataFrame) -> pd.DataFrame:
    keys = ["binder_sequence", "target_sequence"]
    grouped = df.groupby(keys, dropna=False, sort=False)
    n_in = len(df)
    n_groups = grouped.ngroups

    def _agg(g: pd.DataFrame) -> pd.Series:
        kd = pd.to_numeric(g["log10_kd"], errors="coerce")
        return pd.Series(
            {
                "binder_desc": g["binder_desc"].dropna().astype(str).iloc[0] if g["binder_desc"].notna().any() else "",
                "target_desc": g["target_desc"].dropna().astype(str).iloc[0] if g["target_desc"].notna().any() else "",
                "log10_kd": kd.median(),
                "log10_kd_lower": pd.to_numeric(g["log10_kd_lower"], errors="coerce").median(),
                "log10_kd_upper": pd.to_numeric(g["log10_kd_upper"], errors="coerce").median(),
                "log10_kd_spread": float(kd.max() - kd.min()) if kd.notna().any() else pd.NA,
                "n_replicates": len(g),
                "source": ",".join(sorted(g["source"].dropna().unique())),
                "assay_type": g["assay_type"].iloc[0],
            }
        )

    out = grouped.apply(_agg, include_groups=False).reset_index()
    out["kd_nm"] = 10.0 ** pd.to_numeric(out["log10_kd"], errors="coerce")
    out["kd_molar"] = out["kd_nm"] * 1e-9
    out["censored"] = pd.to_numeric(out["log10_kd"], errors="coerce").isna()
    print(
        f"Collapsed technical replicates: {n_in:,} rows → {n_groups:,} unique pairs "
        f"(median log10_kd; log10_kd_spread is max−min)."
    )
    extra = [c for c in ("log10_kd_spread", "n_replicates") if c in out.columns]
    return out[EXPORT_COLUMNS + extra]


def cmd_alphaseq(args: argparse.Namespace) -> Path:
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    configs = _list_alphaseq_configs(cache)
    print(f"AlphaSeq configs ({len(configs)}): {configs}")

    frames: list[pd.DataFrame] = []
    remaining = args.max_rows
    for name in configs:
        raw = _load_alphaseq_config(name, cache)
        if remaining is not None:
            raw = raw.head(remaining)
        norm = _normalize_alphaseq(raw, source=name)
        if norm is None:
            continue
        frames.append(norm)
        if remaining is not None:
            remaining -= len(norm)
            if remaining <= 0:
                break

    if not frames:
        raise SystemExit("No AlphaSeq configs produced a valid table.")

    df = pd.concat(frames, ignore_index=True)
    if args.max_rows is not None:
        df = df.head(args.max_rows)

    n_before_binders = len(df)
    n_censored = int(df["censored"].sum())
    if args.binders_only:
        df = df.loc[~df["censored"]].copy()
        print(
            f"--binders-only: dropped {n_before_binders - len(df):,} censored rows "
            f"({n_censored:,} NaN Kd). Remaining {len(df):,}."
        )
    else:
        print(
            f"Keeping {n_censored:,} censored (NaN Kd) rows as hard negatives. "
            "Pass --binders-only to drop them."
        )

    if args.collapse_replicates:
        df = _collapse_replicates(df)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    _distribution_summary(df, f"wrote {out}")
    return out


def _strip_url(url: str) -> str:
    return url.replace("%20", "").replace(" ", "")


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"cache hit {dest} ({dest.stat().st_size:,} bytes)")
        return dest
    clean = _strip_url(url)
    print(f"GET {clean}")
    with requests.get(clean, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with tmp.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                if chunk:
                    fh.write(chunk)
        tmp.replace(dest)
    print(f"  wrote {dest.stat().st_size:,} bytes")
    return dest


def _preview_text(member: tarfile.TarInfo, extracted: bytes, n: int = 8) -> None:
    text = extracted.decode("utf-8", errors="replace")
    lines = text.splitlines()
    print(f"  --- {member.name} ({member.size:,} bytes, {len(lines)} lines) ---")
    for line in lines[:n]:
        print(f"  {line[:200]}")
    if len(lines) > n:
        print(f"  … {len(lines) - n} more lines")


def cmd_cao(args: argparse.Namespace) -> None:
    cache = Path(args.cache) / "cao2022"
    archive = _download(CAO_BASE + CAO_ARCHIVE, cache / CAO_ARCHIVE)

    if not args.inspect and not args.export:
        print("Cao parser is inspect-first. Re-run with --inspect to see archive headers.")
        print(f"Cached archive: {archive}")
        return

    print(f"\n=== Cao archive listing: {archive.name} ===")
    text_suffixes = {".csv", ".tsv", ".txt", ".sc", ".dat"}
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        for m in members:
            kind = "dir" if m.isdir() else "file"
            print(f"  {kind:4} {m.size:12,}  {m.name}")

        if args.inspect:
            print("\n=== Previews of tabular files ===")
            for m in members:
                if not m.isfile():
                    continue
                suffix = Path(m.name).suffix.lower()
                if suffix not in text_suffixes:
                    continue
                extracted = tar.extractfile(m)
                if extracted is None:
                    continue
                raw = extracted.read(64_000)
                _preview_text(m, raw)

    if args.export:
        print(
            "\nWARNING: Cao export is not enabled until --inspect headers are reviewed.\n"
            "Do not guess a parser. Re-run after confirming column names.",
            file=sys.stderr,
        )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    al = sub.add_parser("alphaseq", help="Build the open-alphaseq affinity CSV")
    al.add_argument(
        "--out",
        default=str(DEFAULT_OUT_DIR / "alphaseq_affinity.csv"),
        help="Output CSV path",
    )
    al.add_argument("--cache", default=str(DEFAULT_CACHE), help="Download / parquet cache")
    al.add_argument("--max-rows", type=int, default=None, help="Cap rows (smoke tests)")
    al.add_argument(
        "--binders-only",
        action="store_true",
        help="Drop censored (NaN Kd) rows. Opt-in; default keeps hard negatives.",
    )
    al.add_argument(
        "--collapse-replicates",
        action="store_true",
        help="Median-aggregate duplicate (binder, target) pairs; report spread.",
    )
    al.set_defaults(func=cmd_alphaseq)

    cao = sub.add_parser("cao", help="Cao et al. 2022 IPD dump (inspect first)")
    cao.add_argument("--cache", default=str(DEFAULT_CACHE), help="Download cache")
    cao.add_argument(
        "--inspect",
        action="store_true",
        help="Download (cached), list the archive, preview tabular files",
    )
    cao.add_argument(
        "--export",
        action="store_true",
        help="Export CSV (refuses until inspect headers are wired)",
    )
    cao.add_argument("--out", default=str(DEFAULT_OUT_DIR / "cao_affinity.csv"))
    cao.set_defaults(func=cmd_cao)
    return p


def main(argv: list[str] | None = None) -> None:
    warnings.filterwarnings("ignore", category=FutureWarning)
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
