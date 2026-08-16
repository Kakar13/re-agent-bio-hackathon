#!/usr/bin/env python3
"""Join Cao affinity design names to sequences in scripts_and_main_pdbs.tar.gz."""

from __future__ import annotations

import json
import tarfile
from io import TextIOWrapper
from pathlib import Path

from eda_cao2022 import AFFINITY_PREFIX, ARCHIVE, SKIP_SUBDIR, parse_sc

SEQ_ARC = Path("data/raw/cao2022/scripts_and_main_pdbs.tar.gz")
OUT = Path("results/cao2022/sequence_join.json")


def load_sequences() -> dict[str, dict]:
    """name -> {binder, target, file, ntok}."""
    recs: dict[str, dict] = {}
    with tarfile.open(SEQ_ARC, "r:gz") as tar:
        for m in tar.getmembers():
            if not m.isfile() or not m.name.endswith(".seq"):
                continue
            if "design_models_sequence" not in m.name:
                continue
            src = Path(m.name).stem
            fh = tar.extractfile(m)
            if fh is None:
                continue
            text = TextIOWrapper(fh, encoding="utf-8", errors="replace").read()
            for line in text.splitlines():
                parts = line.split()
                if len(parts) < 2:
                    continue
                name = parts[-1]
                binder = parts[0]
                antigen = parts[1] if len(parts) >= 3 else ""
                recs[name] = {
                    "binder": binder,
                    "target_seq": antigen,
                    "file": src,
                    "ntok": len(parts),
                    "binder_len": len(binder),
                    "target_len": len(antigen),
                }
    return recs


def load_affinity_names() -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        for m in tar.getmembers():
            if not m.isfile() or not m.name.endswith(".sc"):
                continue
            if not m.name.startswith(AFFINITY_PREFIX) or SKIP_SUBDIR in m.name:
                continue
            extracted = tar.extractfile(m)
            if extracted is None:
                continue
            text = TextIOWrapper(extracted, encoding="utf-8", errors="replace").read()
            frame = parse_sc(text)
            if frame.empty:
                continue
            tgt = Path(m.name).stem
            for desc in frame["description"].astype(str):
                pairs.append((tgt, desc))
    return pairs


def main() -> None:
    print("loading sequences…")
    seqs = load_sequences()
    print(f"  {len(seqs):,} unique sequence names")
    print("loading affinity names…")
    aff = load_affinity_names()
    print(f"  {len(aff):,} affinity rows")

    exact = 0
    prefix = 0
    none = 0
    by_target: dict[str, dict] = {}
    samples = []

    # also index seq names for prefix: affinity desc is often a prefix of the seq name
    # or seq name is a prefix of affinity desc
    for tgt, desc in aff:
        slot = by_target.setdefault(
            tgt,
            {"n": 0, "exact": 0, "prefix": 0, "none": 0},
        )
        slot["n"] += 1
        hit = seqs.get(desc)
        how = "exact"
        if hit is None:
            # try: desc is prefix of a seq name, or vice versa
            # too slow to scan 1M names per row — use startswith via a few heuristics
            how = None
            if desc in seqs:
                hit = seqs[desc]
                how = "exact"
        if hit is None:
            none += 1
            slot["none"] += 1
        elif how == "exact":
            exact += 1
            slot["exact"] += 1
            if len(samples) < 8:
                samples.append(
                    {
                        "target": tgt,
                        "description": desc[:72],
                        "binder_len": hit["binder_len"],
                        "target_len": hit["target_len"],
                        "binder_head": hit["binder"][:18] + "…",
                        "match": "exact",
                    }
                )
        else:
            prefix += 1
            slot["prefix"] += 1

    # Second pass: unmatched affinity names vs seq names sharing a long prefix.
    # Build a map of stripped SSM parent (drop __pos__AA) and truncated names.
    unmatched_by_tgt = {t: s["none"] for t, s in by_target.items()}

    seq_files: dict[str, int] = {}
    binder_lens = []
    for rec in seqs.values():
        seq_files[rec["file"]] = seq_files.get(rec["file"], 0) + 1
        binder_lens.append(rec["binder_len"])

    summary = {
        "seq_archive": str(SEQ_ARC),
        "n_seq_names": len(seqs),
        "n_affinity_rows": len(aff),
        "exact_match": exact,
        "exact_pct": exact / len(aff) if aff else 0,
        "unmatched": none,
        "unmatched_pct": none / len(aff) if aff else 0,
        "seq_files": dict(sorted(seq_files.items())),
        "binder_len": {
            "min": min(binder_lens),
            "median": sorted(binder_lens)[len(binder_lens) // 2],
            "max": max(binder_lens),
        },
        "by_target": [
            {
                "target": t,
                **s,
                "exact_pct": s["exact"] / s["n"] if s["n"] else 0,
            }
            for t, s in sorted(by_target.items(), key=lambda kv: -kv[1]["n"])
        ],
        "samples": samples,
        "note": (
            "Sequences live in scripts_and_main_pdbs.tar.gz "
            "under design_models_sequence/*.seq (binder target name). "
            "The 68 GB design_models_pdb.tar.gz is 3D models, not required for the sequence join."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2))
    print(f"wrote {OUT}")
    print(
        f"exact {exact:,} ({exact/len(aff):.1%})  "
        f"unmatched {none:,} ({none/len(aff):.1%})"
    )
    print("worst unmatched targets:")
    for row in sorted(summary["by_target"], key=lambda r: -r["none"])[:8]:
        print(f"  {row['target']:22} n={row['n']:7,} exact={row['exact']:7,} none={row['none']:7,}")


if __name__ == "__main__":
    main()
