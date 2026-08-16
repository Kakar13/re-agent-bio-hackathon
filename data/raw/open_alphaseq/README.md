# aalphabio/open-alphaseq

Local copy of [aalphabio/open-alphaseq](https://huggingface.co/datasets/aalphabio/open-alphaseq) (MIT). Nine AlphaSeq yeast-mating PPI experiments. Parquet shards are gitignored.

```bash
uv run python scripts/download_open_alphaseq.py
uv run python scripts/eda_open_alphaseq.py
```

Label is `alphaseq_affinity` = log₁₀(estimated Kd in nM). Lower is tighter binding (0 = 1 nM, 3 = 1 µM). This is **binding affinity**, not folding ΔG.
