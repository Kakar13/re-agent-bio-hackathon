"""Smoke test for the Boltz-2 on Modal structure stage.

Two modes:
  tiny  - minimal two-chain complex, no MSA, few steps. Validates that Proto can
          deploy to Modal and return a structure at all, and measures cold start.
  real  - mature cathepsin S (UniProt P25774, chain 115-331) co-folded with a CD74
          segment spanning the CLIP boundary (P04233, residues 105-135).

Usage: uv run python scripts/smoke_boltz_ctss.py [tiny|real]
"""

import json
import sys
import time
from pathlib import Path

from proto_tools.tools.structure_prediction.boltz2.boltz2 import (
    Boltz2Config,
    Boltz2Input,
    run_boltz2,
)

# UniProt P25774 residues 115-331: the mature, catalytically active enzyme.
# Catalytic triad in this numbering: Cys25, His164, Asn184.
CTSS_MATURE = (
    "LPDSVDWREKGCVTEVKYQGSCGACWAFSAVGALEAQLKLKTGKLVSLSAQNLVDCSTEKYGNKGCNGGFMTTAFQYII"
    "DNKGIDSDASYPYKAMDQKCQYDSKYRAATCSKYTELPYGREDVLKEAVANKGPVSVGVDARHPSFFLYRSGVYYEPSC"
    "TQNVNHGVLVVGYGDLNGKEYWLVKNSWGHNFGEEGYIRMARNKGNHCGIASFPSYPEI"
)

# UniProt P04233 residues 105-135. CLIP is annotated 97-120, so this spans the
# C-terminal CLIP boundary (M120 | G121) with room either side for the protease
# to engage P4-P4'.
CD74_SEGMENT = "SKMRMATPLLMQALPMGALPQGPMQNATKYG"
CD74_SEGMENT_START = 105

OUT_DIR = Path("results/smoke")


def extract_cif(structure) -> str | None:
    """Structure carries CIF text, but the attribute name varies by version."""
    for attr in ("cif", "content", "structure", "data", "text"):
        val = getattr(structure, attr, None)
        if isinstance(val, str) and len(val) > 100:
            return val
    for meth in ("to_cif", "as_cif", "get_cif"):
        fn = getattr(structure, meth, None)
        if callable(fn):
            try:
                val = fn()
                if isinstance(val, str):
                    return val
            except Exception:
                pass
    return None


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "tiny"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if mode == "tiny":
        chains = [
            {"sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", "entity_type": "protein"},
            {"sequence": "SKMRMATPLLMQALPM", "entity_type": "protein"},
        ]
        config = Boltz2Config(
            device="modal",
            use_msa=False,
            recycling_steps=1,
            sampling_steps=25,
            verbose=2,
            timeout=1800,
        )
    elif mode == "real":
        chains = [
            {"sequence": CTSS_MATURE, "entity_type": "protein"},
            {"sequence": CD74_SEGMENT, "entity_type": "protein"},
        ]
        config = Boltz2Config(device="modal", use_msa=True, verbose=2, timeout=3600)
    else:
        print(f"unknown mode {mode!r}; use tiny or real")
        return 2

    print(f"[{mode}] chain lengths: {[len(c['sequence']) for c in chains]}", flush=True)
    print(f"[{mode}] submitting to Modal (first call deploys the app; not a hang)", flush=True)

    started = time.time()
    out = run_boltz2(Boltz2Input(complexes=[{"chains": chains}]), config)
    elapsed = time.time() - started

    print(f"[{mode}] returned in {elapsed:.1f}s success={out.success} errors={out.errors}")
    if not out.structures:
        print(f"[{mode}] SMOKE-FAIL no structures returned")
        return 1

    struct = out.structures[0]
    m = struct.metrics
    summary = {
        "mode": mode,
        "elapsed_s": round(elapsed, 1),
        "chain_lengths": [len(c["sequence"]) for c in chains],
        "confidence_score": getattr(m, "confidence_score", None),
        "ptm": getattr(m, "ptm", None),
        "iptm": getattr(m, "iptm", None),
        "complex_plddt": getattr(m, "complex_plddt", None),
        "avg_pae": getattr(m, "avg_pae", None),
        "pair_chains_iptm": getattr(m, "pair_chains_iptm", None),
    }
    if mode == "real":
        summary["cd74_segment_start"] = CD74_SEGMENT_START
        summary["scissile_bond"] = "M120|G121 (CLIP C-terminal boundary)"

    print(json.dumps(summary, indent=2))
    (OUT_DIR / f"{mode}_metrics.json").write_text(json.dumps(summary, indent=2))

    cif = extract_cif(struct)
    if cif:
        path = OUT_DIR / f"{mode}_model.cif"
        path.write_text(cif)
        print(f"[{mode}] structure written to {path} ({len(cif)} chars)")
    else:
        print(f"[{mode}] could not locate CIF text; Structure attrs: {dir(struct)}")

    print(f"[{mode}] SMOKE-OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
