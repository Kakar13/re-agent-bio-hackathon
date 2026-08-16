"""CLI: run the full Track 1 dataset build.

    uv run python -m re_agent.e2e_pls.build_dataset --target-rows 3000000
    uv run python -m re_agent.e2e_pls.build_dataset --expand-by 2000000

Pipeline: fetch UniProt human/viral/bacterial + RCSB de novo parent proteins
and DS613 -> tile into 9-mer windows -> label with MHCflurry (MHC-I binding)
-> quantile-stratify sample down to --target-rows -> Pepsickle (cleavage) on
the kept parents -> assign protein-cluster splits -> merge in DS613 ->
hard-gate against schema.py -> write parquet + manifest + dataset card under
--output-dir (gitignored).

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
DEFAULT_TARGET_ROWS = 3_000_000
DEFAULT_N_HUMAN = 16_000
DEFAULT_N_VIRAL = 6_000
DEFAULT_N_BACTERIAL = 8_000
DEFAULT_N_DE_NOVO = 2_000


def _load_parent_sources(
    n_human_proteins: int,
    n_viral_proteins: int,
    n_bacterial_proteins: int,
    n_de_novo_proteins: int,
    uniprot_min_length: int,
    uniprot_max_length: int,
    rcsb_min_length: int,
    rcsb_max_length: int,
    raw_dir: Path,
    seed: int,
) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    seq_cache = Path(raw_dir) / "sequences"
    fetch_kw = {
        "seed": seed,
        "cache_dir": seq_cache,
        "min_length": uniprot_min_length,
        "max_length": uniprot_max_length,
    }
    length_tag = f"len{uniprot_min_length}-{uniprot_max_length}"
    sources_by_domain = {
        "natural_human": data.fetch_uniprot_human(
            n_human_proteins, cache_name=f"uniprot_human_{length_tag}", **fetch_kw
        ),
        "natural_viral": data.fetch_uniprot_viral(
            n_viral_proteins, cache_name=f"uniprot_viral_{length_tag}", **fetch_kw
        ),
        "natural_bacterial": data.fetch_uniprot_bacterial(
            n_bacterial_proteins, cache_name=f"uniprot_bacterial_{length_tag}", **fetch_kw
        ),
        "de_novo": data.fetch_rcsb_de_novo(
            n_de_novo_proteins,
            min_length=rcsb_min_length,
            max_length=rcsb_max_length,
            seed=seed,
            cache_dir=seq_cache,
            cache_name=f"rcsb_de_novo_len{rcsb_min_length}-{rcsb_max_length}",
        ),
    }
    all_parent_seqs = {pid: seq for seqs in sources_by_domain.values() for pid, seq in seqs.items()}
    return sources_by_domain, all_parent_seqs


def _source_records(
    sources_by_domain: dict[str, dict[str, str]], ds613_provenance: data.SourceProvenance
) -> list[data.SourceProvenance]:
    def _sample_ids(seqs: dict[str, str]) -> str:
        return ", ".join(sorted(seqs)[:3]) + (", ..." if len(seqs) > 3 else "")

    human_seqs = sources_by_domain["natural_human"]
    viral_seqs = sources_by_domain["natural_viral"]
    bacterial_seqs = sources_by_domain["natural_bacterial"]
    de_novo_seqs = sources_by_domain["de_novo"]
    return [
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


def build(
    target_rows: int = DEFAULT_TARGET_ROWS,
    n_human_proteins: int = DEFAULT_N_HUMAN,
    n_viral_proteins: int = DEFAULT_N_VIRAL,
    n_bacterial_proteins: int = DEFAULT_N_BACTERIAL,
    n_de_novo_proteins: int = DEFAULT_N_DE_NOVO,
    uniprot_min_length: int = data.DEFAULT_UNIPROT_MIN_LENGTH,
    uniprot_max_length: int = data.DEFAULT_UNIPROT_MAX_LENGTH,
    rcsb_min_length: int = data.DEFAULT_RCSB_MIN_LENGTH,
    rcsb_max_length: int = data.DEFAULT_RCSB_MAX_LENGTH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
    seed: int = 0,
) -> dict:
    t0 = time.monotonic()
    print(
        f"fetching {n_human_proteins} human + {n_viral_proteins} viral + "
        f"{n_bacterial_proteins} bacterial + {n_de_novo_proteins} de novo parent proteins..."
    )
    sources_by_domain, all_parent_seqs = _load_parent_sources(
        n_human_proteins,
        n_viral_proteins,
        n_bacterial_proteins,
        n_de_novo_proteins,
        uniprot_min_length,
        uniprot_max_length,
        rcsb_min_length,
        rcsb_max_length,
        raw_dir,
        seed,
    )
    fetch_s = time.monotonic() - t0
    print(
        f"  got {len(sources_by_domain['natural_human'])} human, "
        f"{len(sources_by_domain['natural_viral'])} viral, "
        f"{len(sources_by_domain['natural_bacterial'])} bacterial, "
        f"{len(sources_by_domain['de_novo'])} de novo (fetch {fetch_s:.1f}s)"
    )

    ds613_df, ds613_provenance = data.fetch_ds613(raw_dir / "ds613")

    pool = data.build_candidate_pool(sources_by_domain)
    pool["hla_allele"] = data.DEFAULT_HLA
    print(f"candidate pool: {len(pool)} unique 9-mers from {len(all_parent_seqs)} parents")

    t1 = time.monotonic()
    pool = label.attach_mhc_labels(pool)
    print(f"mhcflurry labeling done ({time.monotonic() - t1:.1f}s)")

    pool = data.quantile_stratified_sample(pool, target_rows, score_col="mhc_percentile", seed=seed)
    print(f"sampled to {len(pool)} rows")

    pool = data.finalize_candidate_columns(pool)
    pool = data.assign_clusters_and_splits(pool, seed=seed)

    kept_ids = set(pool["parent_sequence_id"])
    kept_seqs = {k: v for k, v in all_parent_seqs.items() if k in kept_ids}
    t1 = time.monotonic()
    cleavage_table = label.run_pepsickle_on_proteins(kept_seqs)
    pool = label.attach_cleavage_labels(pool, cleavage_table)
    print(f"pepsickle labeling done on {len(kept_seqs)} parents ({time.monotonic() - t1:.1f}s)")

    ds613_rows = data.build_ds613_rows(ds613_df, seed=seed)
    combined = pd.concat([pool, ds613_rows], ignore_index=True)[list(schema.REQUIRED_COLUMNS)]
    combined = label.attach_mhc_labels(combined, only_missing=True)  # fills DS613 rows only
    combined = label.stamp_label_model_version(combined)
    combined = label.impute_tap_labels(combined)  # spreads DS613 signal via one-hot ridge
    print(f"tap imputed for {combined['tap_log_ic50_relative'].notna().sum()} rows total")

    sources = _source_records(sources_by_domain, ds613_provenance)
    build_config = {
        "target_rows": target_rows,
        "n_human_proteins": n_human_proteins,
        "n_viral_proteins": n_viral_proteins,
        "n_bacterial_proteins": n_bacterial_proteins,
        "n_de_novo_proteins": n_de_novo_proteins,
        "uniprot_min_length": uniprot_min_length,
        "uniprot_max_length": uniprot_max_length,
        "rcsb_min_length": rcsb_min_length,
        "rcsb_max_length": rcsb_max_length,
        "seed": seed,
        "peptide_len": data.DEFAULT_PEPTIDE_LEN,
        "flank_len": data.DEFAULT_FLANK_LEN,
        "hla_allele": data.DEFAULT_HLA,
    }

    manifest = dataset_card.build_dataset_card(combined, output_dir, sources, build_config)
    manifest["duration_s"] = time.monotonic() - t0
    print(f"done in {manifest['duration_s']:.1f}s -- {manifest['n_rows']} rows -> {output_dir}")
    return manifest


def expand(
    add_rows: int,
    existing_path: Path,
    n_human_proteins: int = DEFAULT_N_HUMAN,
    n_viral_proteins: int = DEFAULT_N_VIRAL,
    n_bacterial_proteins: int = DEFAULT_N_BACTERIAL,
    n_de_novo_proteins: int = DEFAULT_N_DE_NOVO,
    uniprot_min_length: int = data.DEFAULT_UNIPROT_MIN_LENGTH,
    uniprot_max_length: int = data.DEFAULT_UNIPROT_MAX_LENGTH,
    rcsb_min_length: int = data.DEFAULT_RCSB_MIN_LENGTH,
    rcsb_max_length: int = data.DEFAULT_RCSB_MAX_LENGTH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
    seed: int = 0,
) -> dict:
    """Append `add_rows` new unique 9-mers to an existing parquet.

    Reuses the same parent-protein caches, drops peptides already in the
    dataset, MHC-labels the leftovers, stratified-samples `add_rows`, and
    inherits each parent protein's existing train/val/test split.
    """
    t0 = time.monotonic()
    existing = pd.read_parquet(existing_path)
    result = schema.validate_dataframe(existing)
    if not result.ok:
        raise SystemExit(f"existing dataset failed schema validation: {result.errors}")
    existing_peptides = set(existing["peptide"].astype(str))
    print(f"existing dataset: {len(existing)} rows, {len(existing_peptides)} unique peptides")

    sources_by_domain, all_parent_seqs = _load_parent_sources(
        n_human_proteins,
        n_viral_proteins,
        n_bacterial_proteins,
        n_de_novo_proteins,
        uniprot_min_length,
        uniprot_max_length,
        rcsb_min_length,
        rcsb_max_length,
        raw_dir,
        seed,
    )
    print(
        f"  parents: {len(sources_by_domain['natural_human'])} human, "
        f"{len(sources_by_domain['natural_viral'])} viral, "
        f"{len(sources_by_domain['natural_bacterial'])} bacterial, "
        f"{len(sources_by_domain['de_novo'])} de novo"
    )

    pool = data.build_candidate_pool(sources_by_domain)
    leftover = pool.loc[~pool["peptide"].isin(existing_peptides)].copy()
    leftover["hla_allele"] = data.DEFAULT_HLA
    print(f"leftover unique 9-mers: {len(leftover)}")
    if len(leftover) == 0:
        raise SystemExit("no unused 9-mers left in the parent pool; fetch more proteins")

    t1 = time.monotonic()
    leftover = label.attach_mhc_labels(leftover)
    print(f"mhcflurry labeling done ({time.monotonic() - t1:.1f}s)")

    added = data.quantile_stratified_sample(
        leftover, add_rows, score_col="mhc_percentile", seed=seed + 1
    )
    print(f"sampled {len(added)} new rows")

    added = data.finalize_candidate_columns(added)
    added = data.inherit_or_assign_splits(added, existing, seed=seed + 1)

    kept_ids = set(added["parent_sequence_id"])
    kept_seqs = {k: v for k, v in all_parent_seqs.items() if k in kept_ids}
    t1 = time.monotonic()
    cleavage_table = label.run_pepsickle_on_proteins(kept_seqs)
    added = label.attach_cleavage_labels(added, cleavage_table)
    print(f"pepsickle labeling done on {len(kept_seqs)} parents ({time.monotonic() - t1:.1f}s)")

    added = added[list(schema.REQUIRED_COLUMNS)]
    added = label.stamp_label_model_version(added)
    combined = pd.concat([existing, added], ignore_index=True)
    combined = label.impute_tap_labels(combined)
    print(f"combined: {len(combined)} rows")

    _, ds613_provenance = data.fetch_ds613(raw_dir / "ds613")
    sources = _source_records(sources_by_domain, ds613_provenance)
    teacher_n = int((combined["label_origin"] == "teacher").sum())
    build_config = {
        "mode": "expand",
        "expand_by": add_rows,
        "existing_path": str(existing_path),
        "existing_n_rows": len(existing),
        "n_human_proteins": n_human_proteins,
        "n_viral_proteins": n_viral_proteins,
        "n_bacterial_proteins": n_bacterial_proteins,
        "n_de_novo_proteins": n_de_novo_proteins,
        "uniprot_min_length": uniprot_min_length,
        "uniprot_max_length": uniprot_max_length,
        "rcsb_min_length": rcsb_min_length,
        "rcsb_max_length": rcsb_max_length,
        "seed": seed,
        "peptide_len": data.DEFAULT_PEPTIDE_LEN,
        "flank_len": data.DEFAULT_FLANK_LEN,
        "hla_allele": data.DEFAULT_HLA,
        "target_rows": teacher_n,
    }

    manifest = dataset_card.build_dataset_card(combined, output_dir, sources, build_config)
    manifest["duration_s"] = time.monotonic() - t0
    print(f"done in {manifest['duration_s']:.1f}s -- {manifest['n_rows']} rows -> {output_dir}")
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-rows", type=int, default=DEFAULT_TARGET_ROWS)
    parser.add_argument(
        "--expand-by",
        type=int,
        default=None,
        help="append this many new unique 9-mers to the existing parquet instead of rebuilding",
    )
    parser.add_argument(
        "--existing",
        default=None,
        help="parquet to expand (default: <output-dir>/dataset.parquet)",
    )
    parser.add_argument("--n-human-proteins", type=int, default=DEFAULT_N_HUMAN)
    parser.add_argument("--n-viral-proteins", type=int, default=DEFAULT_N_VIRAL)
    parser.add_argument("--n-bacterial-proteins", type=int, default=DEFAULT_N_BACTERIAL)
    parser.add_argument("--n-de-novo-proteins", type=int, default=DEFAULT_N_DE_NOVO)
    parser.add_argument("--uniprot-min-length", type=int, default=data.DEFAULT_UNIPROT_MIN_LENGTH)
    parser.add_argument("--uniprot-max-length", type=int, default=data.DEFAULT_UNIPROT_MAX_LENGTH)
    parser.add_argument("--rcsb-min-length", type=int, default=data.DEFAULT_RCSB_MIN_LENGTH)
    parser.add_argument("--rcsb-max-length", type=int, default=data.DEFAULT_RCSB_MAX_LENGTH)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--raw-dir", default=str(DEFAULT_RAW_DIR))
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    common = dict(
        n_human_proteins=args.n_human_proteins,
        n_viral_proteins=args.n_viral_proteins,
        n_bacterial_proteins=args.n_bacterial_proteins,
        n_de_novo_proteins=args.n_de_novo_proteins,
        uniprot_min_length=args.uniprot_min_length,
        uniprot_max_length=args.uniprot_max_length,
        rcsb_min_length=args.rcsb_min_length,
        rcsb_max_length=args.rcsb_max_length,
        output_dir=Path(args.output_dir),
        raw_dir=Path(args.raw_dir),
        seed=args.seed,
    )
    if args.expand_by is not None:
        existing = Path(args.existing) if args.existing else Path(args.output_dir) / "dataset.parquet"
        expand(add_rows=args.expand_by, existing_path=existing, **common)
    else:
        build(target_rows=args.target_rows, **common)


if __name__ == "__main__":
    main()
