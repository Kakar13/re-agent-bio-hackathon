"""Download aalphabio/open-alphaseq parquet into data/raw/."""

from __future__ import annotations

import urllib.request
from pathlib import Path

ROOT = Path("data/raw/open_alphaseq")
BASE = "https://huggingface.co/datasets/aalphabio/open-alphaseq/resolve/main"
EXPERIMENTS = [
    "YM_0005",
    "YM_0549",
    "YM_0693",
    "YM_0852",
    "YM_0985",
    "YM_0988",
    "YM_0989",
    "YM_0990",
    "YM_1068",
]
FILES = ["README.md", "LICENSE"] + [
    item
    for ym in EXPERIMENTS
    for item in (f"data/{ym}/data.parquet", f"data/{ym}/README.md")
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
