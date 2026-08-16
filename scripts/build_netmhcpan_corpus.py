#!/usr/bin/env python3
"""Build an MHCflurry-independent NetMHCpan HLA-A*02:01 corpus locally.

For a cloud-detached build, use ``scripts/build_netmhcpan_corpus_modal.py``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from re_agent.e2e_pls.netmhcpan_corpus import (  # noqa: E402
    build_corpus_artifacts,
    build_pda_corpus_artifacts,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/netmhcpan-a0201"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/netmhcpan-a0201"))
    parser.add_argument("--target-rows", type=int, default=10_000)
    parser.add_argument("--pda-challenge-rows", type=int, default=1_000)
    parser.add_argument(
        "--pda-designs",
        type=Path,
        default=Path("data/processed/pda_designs.parquet"),
    )
    parser.add_argument("--n-human", type=int, default=120)
    parser.add_argument("--n-viral", type=int, default=60)
    parser.add_argument("--n-bacterial", type=int, default=80)
    parser.add_argument("--n-de-novo", type=int, default=40)
    parser.add_argument("--api-batch-size", type=int, default=500)
    parser.add_argument(
        "--pda-only-full",
        action="store_true",
        help="label every PDA parent and write the distinct 5-fold PDA training corpus",
    )
    parser.add_argument("--parent-batch-size", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if args.pda_only_full:
        manifest = build_pda_corpus_artifacts(
            output_dir=args.output_dir,
            cache_dir=args.cache_dir,
            pda_designs_path=args.pda_designs,
            api_batch_size=args.api_batch_size,
            parent_batch_size=args.parent_batch_size,
            seed=args.seed,
        )
    else:
        manifest = build_corpus_artifacts(
            output_dir=args.output_dir,
            cache_dir=args.cache_dir,
            target_rows=args.target_rows,
            pda_challenge_rows=args.pda_challenge_rows,
            n_human_proteins=args.n_human,
            n_viral_proteins=args.n_viral,
            n_bacterial_proteins=args.n_bacterial,
            n_de_novo_proteins=args.n_de_novo,
            pda_designs_path=args.pda_designs,
            api_batch_size=args.api_batch_size,
            seed=args.seed,
        )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
