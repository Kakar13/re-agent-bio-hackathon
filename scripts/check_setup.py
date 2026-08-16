#!/usr/bin/env python3
"""Print a weekend-ready setup report."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from re_agent.checks import immuno_checks, run_checks  # noqa: E402


def report(title: str, checks: list) -> int:
    width = max(len(check.name) for check in checks)
    failed = 0
    print(title)
    print("-" * 56)
    for check in checks:
        mark = "ok" if check.ok else "!!"
        if not check.ok:
            failed += 1
        print(f"[{mark}] {check.name:<{width}}  {check.detail}")
    print("-" * 56)
    return failed


def main() -> int:
    parser = argparse.ArgumentParser(description="re:AGENT setup report")
    parser.add_argument(
        "--immuno", action="store_true", help="also check the immunogenicity pipeline"
    )
    args = parser.parse_args()

    failed = report("re:AGENT setup", run_checks())
    if args.immuno:
        print()
        failed += report("immunogenicity pipeline", immuno_checks())

    if failed:
        print(f"{failed} item(s) still need setup. See README.md.")
        return 1
    print("Ready to build.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
