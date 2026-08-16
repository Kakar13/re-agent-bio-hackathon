#!/usr/bin/env python3
"""Single entry point for the immunogenicity pipeline.

`pyproject.toml` sets `package = false`, so `src/` is never installed into the
environment and `python -m re_agent.immuno.<stage>` cannot find the package
without an exported PYTHONPATH. This wrapper locates `src` relative to itself and
forwards to each stage's own `main`, so every stage is reachable from a clean
shell. Arguments after the stage name are passed straight through.

    uv run python scripts/immuno.py data
    uv run python scripts/immuno.py train --arm mean_teacher --epochs 30
    uv run python scripts/immuno.py report data/raw/immuno_tests/flu_m1_window.fasta
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

STAGES = {
    "data": ("re_agent.immuno.data", "build the labeled, de novo, and reference window tables"),
    "embed": ("re_agent.immuno.embed", "cache frozen ESM-2 embeddings for all three tables"),
    "train": ("re_agent.immuno.train", "train the baseline and Mean Teacher arms"),
    "validate": ("re_agent.immuno.validate", "run the validation experiments"),
    "report": ("re_agent.immuno.report", "score a FASTA or raw sequence"),
    "figures": ("re_agent.immuno.figures", "render result figures"),
}


def usage() -> str:
    width = max(len(s) for s in STAGES)
    lines = [f"usage: {Path(sys.argv[0]).name} <stage> [args...]", "", "stages:"]
    lines += [f"  {name:<{width}}  {desc}" for name, (_, desc) in STAGES.items()]
    return "\n".join(lines)


def main() -> int:
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(usage())
        return 0 if len(sys.argv) > 1 else 2

    stage = sys.argv[1]
    if stage not in STAGES:
        print(f"unknown stage {stage!r}\n\n{usage()}", file=sys.stderr)
        return 2

    module_name, _ = STAGES[stage]
    module = importlib.import_module(module_name)
    # Rewrite argv so each stage's own argparse sees a sensible program name.
    sys.argv = [f"immuno {stage}", *sys.argv[2:]]
    return module.main() or 0


if __name__ == "__main__":
    raise SystemExit(main())
