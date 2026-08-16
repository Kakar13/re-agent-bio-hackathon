"""Retarget Proto's Modal GPU profile to tiers that need no payment method.

Proto ships Boltz-2 pinned to ``["H100:1", "H200:1", "A100-80GB:1"]``. Modal
gates all three behind a card even when the account holds credits, so the deploy
dies with "Please add a payment method to use H100 GPU functions."

Probing this workspace, T4 (16 GB), L4 (23 GB) and A10 (23 GB) all run; A100 in
either size, L40S, H100 and H200 do not. A 23 GB card is enough for the ~250
residue protease-plus-segment co-folds this pipeline submits, so the profile is
rewritten to prefer L4, then A10, then T4.

Idempotent. ``--restore`` puts the original list back.

    uv run python scripts/patch_proto_gpu.py [--restore]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import proto_tools

TARGET = Path(proto_tools.__file__).parent / "modal" / "gpu_profiles.py"
BACKUP = TARGET.with_suffix(".py.orig")

ORIGINAL = 'GPU_DEFAULT: Final[list[str]] = ["H100:1", "H200:1", "A100-80GB:1"]'
PATCHED = (
    '# PATCHED by scripts/patch_proto_gpu.py: H100/H200/A100 need a Modal\n'
    '# payment method; L4/A10/T4 do not and fit a ~250-residue co-fold.\n'
    'GPU_DEFAULT: Final[list[str]] = ["L4:1", "A10:1", "T4:1"]'
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--restore", action="store_true", help="undo the patch")
    args = ap.parse_args()

    text = TARGET.read_text()

    if args.restore:
        if BACKUP.exists():
            TARGET.write_text(BACKUP.read_text())
            print(f"restored {TARGET}")
        else:
            print("no backup found; nothing to restore")
        return 0

    if PATCHED.splitlines()[-1] in text:
        print(f"already patched: {TARGET}")
        return 0
    if ORIGINAL not in text:
        print(f"unexpected content in {TARGET}; refusing to patch", file=sys.stderr)
        print("looked for:", ORIGINAL, file=sys.stderr)
        return 1

    if not BACKUP.exists():
        BACKUP.write_text(text)
    TARGET.write_text(text.replace(ORIGINAL, PATCHED))
    print(f"patched {TARGET}\n  GPU_DEFAULT -> L4:1, A10:1, T4:1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
