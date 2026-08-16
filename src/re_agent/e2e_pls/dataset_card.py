"""Track 1: dataset manifest, statistics, dataset card, and quality gates.

`build_dataset_card` is the hard gate: it re-runs `schema.validate_dataframe`
and raises on any error (duplicate row IDs, coordinate/flank inconsistency,
missing provenance, split-group overlap, noncanonical residues, invalid
enums) before writing anything. Everything it writes goes under gitignored
`results/dataset/` -- only code and a manifest are meant to be committed,
never the parquet or raw source mirrors.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from re_agent.e2e_pls import schema
from re_agent.e2e_pls.data import SourceProvenance

LABEL_COLUMNS = (
    "cleave_n_prob",
    "cleave_c_prob",
    "tap_log_ic50_relative",
    "mhc_affinity_nm",
    "mhc_percentile",
)

CLAIM_DISCLAIMER = (
    "Experimental cleavage/TAP/MHC-I surrogate labels from teacher models (Pepsickle, MHCflurry) "
    "plus 613 measured TAP-binding peptides (DS613). Not measured immunogenicity, T-cell response, "
    "or presentation ground truth."
)


def _sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _compute_statistics(df: pd.DataFrame) -> dict:
    from re_agent.e2e_pls.label import TAP_RIDGE_TAG

    n_tap_measured = int((df["study_id"] == "ds613").sum())
    tap_imputed_mask = df["label_model_version"].fillna("").str.contains(TAP_RIDGE_TAG)
    n_tap_imputed = int(tap_imputed_mask.sum())

    return {
        "by_split": df["split"].value_counts().to_dict(),
        "by_source_domain": df["source_domain"].value_counts().to_dict(),
        "by_label_origin": df["label_origin"].value_counts().to_dict(),
        "n_unique_parent_sequences": int(df["parent_sequence_id"].nunique()),
        "n_unique_hla_alleles": int(df["hla_allele"].nunique()),
        "tap_provenance": {
            "measured_ds613": n_tap_measured,
            "imputed_by_ridge": n_tap_imputed,
            "missing": int(df["tap_log_ic50_relative"].isna().sum()),
        },
        "label_coverage": {
            col: {"n_present": int(df[col].notna().sum()), "n_missing": int(df[col].isna().sum())}
            for col in LABEL_COLUMNS
        },
        "label_distributions": {
            col: (
                {
                    "mean": float(df[col].mean()),
                    "std": float(df[col].std()),
                    "min": float(df[col].min()),
                    "max": float(df[col].max()),
                }
                if df[col].notna().any()
                else None
            )
            for col in LABEL_COLUMNS
        },
    }


def _render_card(manifest: dict, stats: dict) -> str:
    lines = [
        "# E2E-PLS Dataset Card",
        "",
        f"- Schema version: `{manifest['schema_version']}`",
        f"- Dataset version hash: `{manifest['dataset_version_hash']}`",
        f"- Built at: {manifest['built_at']}",
        f"- Rows: {manifest['n_rows']}",
        f"- Parquet checksum (sha256): `{manifest['parquet_checksum_sha256']}`",
        "",
        "## Sources",
        "",
    ]
    for s in manifest["sources"]:
        lines.append(f"- **{s['name']}** ({s['n_records']} records) -- {s['license']}")
        lines.append(f"  - {s['url']}")
        if s.get("notes"):
            lines.append(f"  - {s['notes']}")

    lines += ["", "## Split sizes", ""]
    for split, n in stats["by_split"].items():
        lines.append(f"- {split}: {n}")

    lines += ["", "## Source domain", ""]
    for domain, n in stats["by_source_domain"].items():
        lines.append(f"- {domain}: {n}")

    lines += [
        "",
        "## Leakage audit",
        "",
        "PASSED -- no protein/peptide cluster spans multiple splits (hard build gate).",
    ]

    lines += ["", "## Label coverage", ""]
    for col, cov in stats["label_coverage"].items():
        lines.append(f"- {col}: {cov['n_present']} present, {cov['n_missing']} missing")

    tap_prov = stats.get("tap_provenance", {})
    if tap_prov:
        lines += [
            "",
            "## TAP provenance",
            "",
            f"- Measured (DS613): {tap_prov['measured_ds613']}",
            f"- Imputed by ridge (one-hot 9-mer, trained on DS613): {tap_prov['imputed_by_ridge']}",
            f"- Missing: {tap_prov['missing']}",
            "",
            "Imputed TAP values spread the 613-peptide measured signal across the "
            "dataset for coverage; they are NOT new measurements. Downstream TapHead "
            "training filters to measured-only rows to avoid learning its own imputer.",
        ]

    lines += ["", "## Label distributions", ""]
    for col, dist in stats["label_distributions"].items():
        if dist is None:
            lines.append(f"- {col}: no values present")
            continue
        line = f"- {col}: mean={dist['mean']:.4f} std={dist['std']:.4f} "
        line += f"min={dist['min']:.4f} max={dist['max']:.4f}"
        lines.append(line)

    if manifest.get("validation_warnings"):
        lines += ["", "## Validation warnings (non-fatal)", ""]
        lines += [f"- {w}" for w in manifest["validation_warnings"]]

    lines += ["", "## Claim disclaimer", "", CLAIM_DISCLAIMER]
    return "\n".join(lines) + "\n"


def build_dataset_card(
    df: pd.DataFrame,
    output_dir: str | Path,
    sources: list[SourceProvenance],
    build_config: dict,
) -> dict:
    """Hard-gates `df` against `schema.validate_dataframe`, then writes the
    parquet, manifest, dataset card, and statistics. Returns the manifest.
    """
    result = schema.validate_dataframe(df)
    if not result.ok:
        raise ValueError(f"dataset failed quality gates: {result.errors}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    parquet_path = output_dir / "dataset.parquet"
    df.to_parquet(parquet_path, index=False)
    checksum = _sha256_file(parquet_path)
    stats = _compute_statistics(df)

    manifest = {
        "schema_version": schema.SCHEMA_VERSION,
        "dataset_version_hash": schema.dataset_version_hash(df),
        "built_at": datetime.now(UTC).isoformat(),
        "n_rows": len(df),
        "parquet_checksum_sha256": checksum,
        "parquet_path": str(parquet_path),
        "build_config": build_config,
        "sources": [asdict(s) for s in sources],
        "validation_warnings": result.warnings,
        "statistics": stats,
    }
    (output_dir / "dataset_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    (output_dir / "dataset_card.md").write_text(_render_card(manifest, stats))
    return manifest
