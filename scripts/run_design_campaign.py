#!/usr/bin/env python3
"""Plan or execute the bounded Proto design campaign."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from re_agent.design.campaign import ProtoLocalRunner, run_campaign  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spec",
        type=Path,
        default=ROOT / "docs" / "design_spec.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results" / "design_campaigns" / "il7ra-rfd3-20260815",
    )
    parser.add_argument(
        "--approve-compute",
        action="store_true",
        help="Actually run GPU tools; otherwise only write the versioned plan",
    )
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    manifest = run_campaign(
        args.spec.resolve(),
        args.output_dir.resolve(),
        ProtoLocalRunner(),
        approved=args.approve_compute,
        device=args.device,
    )
    print(manifest.model_dump_json(indent=2))


if __name__ == "__main__":
    main()
