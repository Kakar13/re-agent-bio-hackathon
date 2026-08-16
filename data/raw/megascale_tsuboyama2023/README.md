# LiteFold/MegaScale-Tsuboyama2023

Local copy of [LiteFold/MegaScale-Tsuboyama2023](https://huggingface.co/datasets/LiteFold/MegaScale-Tsuboyama2023) (CC-BY-4.0). Parquet shards are gitignored.

```bash
uv run python scripts/download_megascale_tsuboyama2023.py
```

Lands here:

- `data/train-*.parquet`, `data/test-*.parquet`
- `metadata/column_mapping.parquet`, `metadata/source_tables.parquet`
- `dataset_summary.json`, `_MANIFEST.json`

EDA: `uv run python scripts/eda_megascale_tsuboyama2023.py`
