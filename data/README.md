# Affinity training data

Flat `sequence → binding affinity` tables for a nanobody/binder model.
Built by `scripts/build_affinity_dataset.py`. Splits: `scripts/make_splits.py`.

Do **not** random-split these rows. Cluster binders (MMseqs2, identity ≥ 0.9,
coverage ≥ 0.8) or hold out whole targets.

## Primary: open-alphaseq

| | |
|---|---|
| Source | [aalphabio/open-alphaseq](https://huggingface.co/datasets/aalphabio/open-alphaseq) |
| License | MIT (academic + commercial) |
| Assay | AlphaSeq yeast mating (`assay_type=alphaseq_yeast_mating`) |
| Label | `log10_kd` = log₁₀(estimated Kd **in nM**). Lower is tighter. |
| Local parquet | `data/raw/open_alphaseq/` (gitignored) |
| Export CSV | `data/processed/affinity/alphaseq_affinity.csv` (gitignored) |

```bash
uv sync --extra affinity
uv run python scripts/build_affinity_dataset.py alphaseq
uv run python scripts/make_splits.py \
  --in data/processed/affinity/alphaseq_affinity.csv \
  --out data/processed/affinity/alphaseq_affinity.splits.csv \
  --by binder
```

### Units

The Hugging Face column `alphaseq_affinity` is log₁₀(Kd in **nanomolar**), not
molar. The prompt formula `kd_molar = 10 ** log10_kd` is therefore off by 10⁹.

| column | definition |
|---|---|
| `log10_kd` | log₁₀(Kd / nM) |
| `kd_nm` | `10 ** log10_kd` |
| `kd_molar` | `kd_nm * 1e-9` |
| `censored` | `True` when `log10_kd` is NaN (no detectable binding — hard negative) |

Default export **keeps** censored rows. `--binders-only` drops them and prints
how many were lost. `--collapse-replicates` median-aggregates duplicate
`(binder_sequence, target_sequence)` pairs and adds `log10_kd_spread`.

### Distribution (full 9-config concat, from local parquet EDA)

| metric | value |
|---|---|
| rows | 1,705,908 |
| observed Kd (`log10_kd` not NaN) | 1,378,079 (80.8%) |
| censored / hard negatives | 327,829 (19.2%) |
| unique binders (`mata_sequence`) | 251,418 |
| unique targets (`matalpha_sequence`) | 2,441 |
| log10 Kd n / min / median / max | 1,378,079 / −0.989 / 3.910 / 7.459 |
| median Kd | ~8.1 µM (`10**3.910` nM) |
| duplicate pairs (technical replicates) | 549,462 |
| usable both-sequences + Kd | 728,126 |

Configs are enumerated at runtime (`datasets.get_dataset_config_names` or the
local `data/raw/open_alphaseq/data/*/data.parquet` cache). Do not hardcode YM
IDs. Published bound-column names (`alphaseq_affinity_*_bound`) are aliased to
the names actually shipped (`affinity_lower_bound` / `affinity_upper_bound`).

This is **parent optimization** (VHH72, PP489, pembrolizumab, trastuzumab, CoV
scFv panel), not de novo design. Iter1 libraries are designed from Iter0 of the
same parent — a random row split leaks.

## Secondary: Cao et al. 2022 (de novo minibinders)

UW IPD dump, not Nature SI. Parsed 27 `{target}.sc` affinity tables from
`experimental_data_and_analysis.tar.gz` (233 MB, 706 files).

```bash
uv run python scripts/build_affinity_dataset.py cao --inspect
uv run python scripts/eda_cao2022.py
```

| metric | value |
|---|---|
| affinity rows | 997,343 |
| unique design names | 995,735 |
| de novo design rows | 754,716 |
| SSM point-mutant rows | 242,627 |
| targets | 27 (12 designed + 15 SSM) |
| both Kd bounds finite | 130,718 (13.1%) |
| lower-only (`ub=inf`) | 224,290 (22.5%) |
| fully unbound (`inf/inf`) | 642,335 (64.4%) |
| kd_mid n / min / median / max | 346,664 / 0.11 nM / 18.7 µM / 1.4 mM |
| sequences | **1,070,060** in `scripts_and_main_pdbs.tar.gz` (`design_models_sequence/*.seq`) |
| affinity ↔ sequence join | **997,343 / 997,343 exact** |

`kd_lb` / `kd_ub` are FACS/NGS enrichment-derived nM estimates, censored at the
sort dynamic range. `kd_mid` is the geometric mean when both bounds are finite.
This is **not** SPR/BLI. True SPR Kd exists only for the ~13 optimized binders
in Table 1 of the paper. Binders are minibinders (36–72 aa, median 64).
`design_models_pdb.tar.gz` (68 GB) and the silent dump (49 GB) are coordinates
— not required for a sequence→Kd table. Never pool with AlphaSeq without
calibration.

## Column contract

| column | meaning |
|---|---|
| `binder_sequence` | A-library (usually VHH / scFv) |
| `target_sequence` | Alpha-library antigen |
| `binder_desc` / `target_desc` | library descriptions |
| `log10_kd` / `_lower` / `_upper` | assay point + bounds |
| `kd_nm` / `kd_molar` | derived; see Units |
| `censored` | NaN Kd |
| `source` | AlphaSeq config name (e.g. `YM_0549`) |
| `assay_type` | `alphaseq_yeast_mating` (or Cao, later) |

Binder/target column order is fixed. Do not swap at train vs inference.
