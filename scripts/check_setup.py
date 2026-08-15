#!/usr/bin/env python3
"""Print a weekend-ready setup report."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from re_agent.checks import run_checks  # noqa: E402


def main() -> int:
    checks = run_checks()
    width = max(len(check.name) for check in checks)
    failed = 0
    print("re:AGENT setup")
    print("-" * 48)
    for check in checks:
        mark = "ok" if check.ok else "!!"
        if not check.ok:
            failed += 1
        print(f"[{mark}] {check.name:<{width}}  {check.detail}")
    print("-" * 48)
    if failed:
        print(f"{failed} item(s) still need setup. See README.md.")
        return 1
    print("Ready to build.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
