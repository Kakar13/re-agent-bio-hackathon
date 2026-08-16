"""Download LiteFold/MegaScale-Tsuboyama2023 parquet into data/raw/."""

from __future__ import annotations

import urllib.request
from pathlib import Path

ROOT = Path("data/raw/megascale_tsuboyama2023")
BASE = "https://huggingface.co/datasets/LiteFold/MegaScale-Tsuboyama2023/resolve/main"
FILES = [
    "dataset_summary.json",
    "_MANIFEST.json",
    "metadata/source_tables.parquet",
    "metadata/column_mapping.parquet",
    "data/test-00000-of-00001.parquet",
    "data/train-00000-of-00002.parquet",
    "data/train-00001-of-00002.parquet",
]


def main() -> None:
    for rel in FILES:
        dest = ROOT / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and dest.stat().st_size > 0:
            print(f"exists {rel} ({dest.stat().st_size} bytes)")
            continue
        url = f"{BASE}/{rel}"
        print(f"GET {rel}")
        urllib.request.urlretrieve(url, dest)
        print(f"  wrote {dest.stat().st_size} bytes")


if __name__ == "__main__":
    main()
