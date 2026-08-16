"""CLI: run the full Track 1 dataset build.

    uv run python -m re_agent.e2e_pls.build_dataset --target-rows 500000

Pipeline: fetch UniProt human/viral/bacterial + RCSB de novo parent proteins
and DS613 -> tile into 9-mer windows -> label with Pepsickle (cleavage) and
MHCflurry (MHC-I binding) -> quantile-stratify sample down to --target-rows
-> assign protein-cluster splits -> merge in DS613 -> hard-gate against
schema.py -> write parquet + manifest + dataset card under --output-dir
(gitignored).

TAP coverage stays at DS613's 613 measured peptides regardless of scale --
no bulk-downloadable, continuous-valued, larger measured TAP source exists
(checked PyPI, IEDB's query API, and MHCBN; see conversation/commit history).
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

from re_agent.e2e_pls import data, dataset_card, label, schema

DEFAULT_OUTPUT_DIR = Path("results/dataset")
DEFAULT_RAW_DIR = Path("data/raw")


def build(
    target_rows: int = 500_000,
    n_human_proteins: int = 1500,
    n_viral_proteins: int = 600,
    n_bacterial_proteins: int = 600,
    n_de_novo_proteins: int = 800,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
    seed: int = 0,
) -> dict:
    t0 = time.monotonic()
    print(
        f"fetching {n_human_proteins} human + {n_viral_proteins} viral + "
        f"{n_bacterial_proteins} bacterial + {n_de_novo_proteins} de novo parent proteins..."
    )
    human_seqs = data.fetch_uniprot_human(n_human_proteins, seed=seed)
    viral_seqs = data.fetch_uniprot_viral(n_viral_proteins, seed=seed)
    bacterial_seqs = data.fetch_uniprot_bacterial(n_bacterial_proteins, seed=seed)
    de_novo_seqs = data.fetch_rcsb_de_novo(n_de_novo_proteins, seed=seed)
    fetch_s = time.monotonic() - t0
    print(
        f"  got {len(human_seqs)} human, {len(viral_seqs)} viral, "
        f"{len(bacterial_seqs)} bacterial, {len(de_novo_seqs)} de novo (fetch {fetch_s:.1f}s)"
    )

    ds613_df, ds613_provenance = data.fetch_ds613(raw_dir / "ds613")

    sources_by_domain = {
        "natural_human": human_seqs,
        "natural_viral": viral_seqs,
        "natural_bacterial": bacterial_seqs,
        "de_novo": de_novo_seqs,
    }
    all_parent_seqs = {pid: seq for seqs in sources_by_domain.values() for pid, seq in seqs.items()}

    pool = data.build_candidate_pool(sources_by_domain)
    pool = data.finalize_candidate_columns(pool)
    print(f"candidate pool: {len(pool)} unique 9-mers")

    t1 = time.monotonic()
    cleavage_table = label.run_pepsickle_on_proteins(all_parent_seqs)
    pool = label.attach_cleavage_labels(pool, cleavage_table)
    print(f"pepsickle labeling done ({time.monotonic() - t1:.1f}s)")

    t1 = time.monotonic()
    pool = label.attach_mhc_labels(pool)
    print(f"mhcflurry labeling done ({time.monotonic() - t1:.1f}s)")

    pool = data.quantile_stratified_sample(pool, target_rows, score_col="mhc_percentile", seed=seed)
    pool = data.assign_clusters_and_splits(pool, seed=seed)
    print(f"sampled to {len(pool)} rows")

    ds613_rows = data.build_ds613_rows(ds613_df, seed=seed)
    combined = pd.concat([pool, ds613_rows], ignore_index=True)[list(schema.REQUIRED_COLUMNS)]
    combined = label.attach_mhc_labels(combined, only_missing=True)  # fills DS613 rows only
    combined = label.stamp_label_model_version(combined)
    combined = label.impute_tap_labels(combined)  # spreads DS613 signal via one-hot ridge
    print(f"tap imputed for {combined['tap_log_ic50_relative'].notna().sum()} rows total")

    def _sample_ids(seqs: dict[str, str]) -> str:
        return ", ".join(sorted(seqs)[:3]) + (", ..." if len(seqs) > 3 else "")

    sources = [
        data.SourceProvenance(
            name="UniProt reviewed human reference",
            url=data.UNIPROT_SEARCH_URL,
            retrieved_at=data.now_iso(),
            sha256="",
            license=data.LICENSE_NATURAL,
            n_records=len(human_seqs),
            notes=f"Sample proteins: {_sample_ids(human_seqs)}",
        ),
        data.SourceProvenance(
            name="UniProt viral reference (human-pathogenic species)",
            url=data.UNIPROT_SEARCH_URL,
            retrieved_at=data.now_iso(),
            sha256="",
            license=data.LICENSE_NATURAL,
            n_records=len(viral_seqs),
            notes=f"organism_ids={data.VIRAL_ORGANISM_IDS}; reviewed+unreviewed. "
            f"Sample: {_sample_ids(viral_seqs)}",
        ),
        data.SourceProvenance(
            name="UniProt bacterial reference (human-pathogenic species)",
            url=data.UNIPROT_SEARCH_URL,
            retrieved_at=data.now_iso(),
            sha256="",
            license=data.LICENSE_NATURAL,
            n_records=len(bacterial_seqs),
            notes=f"organism_ids={data.BACTERIAL_ORGANISM_IDS}; reviewed+unreviewed. "
            f"Sample: {_sample_ids(bacterial_seqs)}",
        ),
        data.SourceProvenance(
            name="RCSB PDB de novo designs",
            url=data.RCSB_SEARCH_URL,
            retrieved_at=data.now_iso(),
            sha256="",
            license=data.LICENSE_DE_NOVO,
            n_records=len(de_novo_seqs),
        ),
        ds613_provenance,
    ]
    build_config = {
        "target_rows": target_rows,
        "n_human_proteins": n_human_proteins,
        "n_viral_proteins": n_viral_proteins,
        "n_bacterial_proteins": n_bacterial_proteins,
        "n_de_novo_proteins": n_de_novo_proteins,
        "seed": seed,
        "peptide_len": data.DEFAULT_PEPTIDE_LEN,
        "flank_len": data.DEFAULT_FLANK_LEN,
        "hla_allele": data.DEFAULT_HLA,
    }

    manifest = dataset_card.build_dataset_card(combined, output_dir, sources, build_config)
    manifest["duration_s"] = time.monotonic() - t0
    print(f"done in {manifest['duration_s']:.1f}s -- {manifest['n_rows']} rows -> {output_dir}")
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-rows", type=int, default=500_000)
    parser.add_argument("--n-human-proteins", type=int, default=1500)
    parser.add_argument("--n-viral-proteins", type=int, default=600)
    parser.add_argument("--n-bacterial-proteins", type=int, default=600)
    parser.add_argument("--n-de-novo-proteins", type=int, default=800)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--raw-dir", default=str(DEFAULT_RAW_DIR))
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    build(
        target_rows=args.target_rows,
        n_human_proteins=args.n_human_proteins,
        n_viral_proteins=args.n_viral_proteins,
        n_bacterial_proteins=args.n_bacterial_proteins,
        n_de_novo_proteins=args.n_de_novo_proteins,
        output_dir=Path(args.output_dir),
        raw_dir=Path(args.raw_dir),
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
