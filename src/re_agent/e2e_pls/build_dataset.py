"""CLI: run the full Track 1 dataset build.

    uv run python -m re_agent.e2e_pls.build_dataset --target-rows 50000

Pipeline: fetch UniProt natural + RCSB de novo parent proteins and DS613 ->
tile into 9-mer windows -> label with Pepsickle (cleavage) and MHCflurry
(MHC-I binding) -> quantile-stratify sample down to --target-rows -> assign
protein-cluster splits -> merge in DS613 -> hard-gate against schema.py ->
write parquet + manifest + dataset card under --output-dir (gitignored).
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
    target_rows: int = 50_000,
    n_natural_proteins: int = 200,
    n_de_novo_proteins: int = 150,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
    seed: int = 0,
) -> dict:
    t0 = time.monotonic()
    print(
        f"fetching {n_natural_proteins} natural + {n_de_novo_proteins} de novo parent proteins..."
    )
    natural_seqs = data.fetch_uniprot_natural(n_natural_proteins, seed=seed)
    de_novo_seqs = data.fetch_rcsb_de_novo(n_de_novo_proteins, seed=seed)
    fetch_s = time.monotonic() - t0
    print(f"  got {len(natural_seqs)} natural, {len(de_novo_seqs)} de novo (fetch {fetch_s:.1f}s)")

    ds613_df, ds613_provenance = data.fetch_ds613(raw_dir / "ds613")

    pool = data.build_candidate_pool(natural_seqs, de_novo_seqs)
    pool = data.finalize_candidate_columns(pool)
    print(f"candidate pool: {len(pool)} unique 9-mers")

    t1 = time.monotonic()
    cleavage_table = label.run_pepsickle_on_proteins({**natural_seqs, **de_novo_seqs})
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

    natural_urls = ", ".join(sorted(natural_seqs)[:3]) + ", ..."
    sources = [
        data.SourceProvenance(
            name="UniProt reviewed human reference",
            url=data.UNIPROT_SEARCH_URL,
            retrieved_at=data.now_iso(),
            sha256="",
            license=data.LICENSE_NATURAL,
            n_records=len(natural_seqs),
            notes=f"Sample proteins: {natural_urls}",
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
        "n_natural_proteins": n_natural_proteins,
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
    parser.add_argument("--target-rows", type=int, default=50_000)
    parser.add_argument("--n-natural-proteins", type=int, default=200)
    parser.add_argument("--n-de-novo-proteins", type=int, default=150)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--raw-dir", default=str(DEFAULT_RAW_DIR))
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    build(
        target_rows=args.target_rows,
        n_natural_proteins=args.n_natural_proteins,
        n_de_novo_proteins=args.n_de_novo_proteins,
        output_dir=Path(args.output_dir),
        raw_dir=Path(args.raw_dir),
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
