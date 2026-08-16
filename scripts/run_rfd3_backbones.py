#!/usr/bin/env python3
"""Generate RFdiffusion3 binder backbones for undruggable targets on Modal.

Baker-lab RFD3 binder settings: hotspots ori, is_non_loopy, step_scale=3, gamma_0=0.2.
Chunked + resumable (16 designs/chunk by default).

Example:
  uv run python scripts/run_rfd3_backbones.py --target pdl1 --n 16 --chunks 1
  uv run python scripts/run_rfd3_backbones.py --target all --n 100
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

from proto_tools.tools.structure_design.rfdiffusion3.rfdiffusion3_sample import (
    RFdiffusion3Config,
    RFdiffusion3DesignSpec,
    RFdiffusion3Input,
    run_rfdiffusion3,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "rfd3_targets.json"
RESULTS_ROOT = ROOT / "results" / "rfd3_binders"


def log(msg: str) -> None:
    print(msg, flush=True)


def load_registry() -> dict[str, Any]:
    if not CONFIG_PATH.exists():
        raise SystemExit(f"Missing {CONFIG_PATH}; run scripts/prep_binder_targets.py first")
    return json.loads(CONFIG_PATH.read_text())


def target_by_name(registry: dict[str, Any], name: str) -> dict[str, Any]:
    for t in registry["targets"]:
        if t["name"] == name:
            return t
    raise SystemExit(f"Unknown target {name!r}; have {[t['name'] for t in registry['targets']]}")


def chunk_dir(target: str, chunk_idx: int) -> Path:
    return RESULTS_ROOT / target / "chunks" / f"chunk_{chunk_idx:03d}"


def backbone_dir(target: str) -> Path:
    return RESULTS_ROOT / target / "backbones"


def manifest_path(target: str) -> Path:
    return RESULTS_ROOT / target / "manifest.jsonl"


def count_existing_backbones(target: str) -> int:
    d = backbone_dir(target)
    if not d.exists():
        return 0
    return len(list(d.glob("*.pdb")))


def chunk_done(target: str, chunk_idx: int) -> bool:
    return (chunk_dir(target, chunk_idx) / "done.json").exists()


def build_specs(target: dict[str, Any]) -> list[RFdiffusion3DesignSpec]:
    target_pdb = ROOT / target["target_pdb"]
    if not target_pdb.exists():
        raise SystemExit(f"Missing target PDB: {target_pdb}")
    specs: list[RFdiffusion3DesignSpec] = []
    for contig in target["contigs"].values():
        specs.append(
            RFdiffusion3DesignSpec(
                input_structure=str(target_pdb),
                contig=contig,
                select_hotspots=target["select_hotspots"],
                infer_ori_strategy="hotspots",
                is_non_loopy=True,
            )
        )
    return specs


def designs_per_chunk(target: dict[str, Any], batch_size: int) -> int:
    return len(target["contigs"]) * batch_size


def save_designs(out: Any, target_name: str, chunk_idx: int, global_start_idx: int) -> list[dict[str, Any]]:
    bb_dir = backbone_dir(target_name)
    bb_dir.mkdir(parents=True, exist_ok=True)
    cdir = chunk_dir(target_name, chunk_idx)
    cdir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []
    design_i = global_start_idx
    for bundle in out.designed_structures:
        for structure_wrap in bundle:
            design_i += 1
            fname = f"{target_name}_bb_{design_i:04d}.pdb"
            dest = bb_dir / fname
            structure_wrap.structure.write_pdb(dest)
            structure_wrap.structure.write_pdb(cdir / fname)
            meta = structure_wrap.metadata or {}
            try:
                rg = structure_wrap.structure.gyration_radius()
            except Exception:
                rg = None
            try:
                ss = structure_wrap.structure.secondary_structure_percentages()
            except Exception:
                ss = None
            records.append(
                {
                    "design_idx": design_i,
                    "pdb": str(dest.relative_to(ROOT)),
                    "spec_key": getattr(bundle, "spec_key", None),
                    "num_residues": structure_wrap.structure.num_residues,
                    "gyration_radius": rg,
                    "secondary_structure": ss,
                    "metadata_keys": sorted(meta.keys()),
                }
            )
    return records


def run_chunk(
    target: dict[str, Any],
    chunk_idx: int,
    *,
    seed: int,
    batch_size: int,
    step_scale: float,
    gamma_0: float,
    num_timesteps: int,
    timeout: int,
    device: str,
    dry_run: bool,
) -> dict[str, Any]:
    name = target["name"]
    if chunk_done(name, chunk_idx):
        log(f"[{name}] chunk {chunk_idx} already done — skipping")
        return json.loads((chunk_dir(name, chunk_idx) / "done.json").read_text())

    specs = build_specs(target)
    n_expected = designs_per_chunk(target, batch_size)
    existing = count_existing_backbones(name)
    log(f"[{name}] chunk {chunk_idx}: {len(specs)} specs x batch={batch_size} → {n_expected} designs")
    log(f"  hotspots={target['select_hotspots']}")
    log(f"  contigs={list(target['contigs'].values())}")
    log(f"  device={device} seed={seed}")

    config = RFdiffusion3Config(
        n_batches=1,
        diffusion_batch_size=batch_size,
        num_timesteps=num_timesteps,
        step_scale=step_scale,
        gamma_0=gamma_0,
        seed=seed,
        timeout=timeout,
        device=device,
        verbose=1,
        prevalidate_inputs=True,
    )
    inputs = RFdiffusion3Input(design_specs=specs)

    run_config = {
        "target": name,
        "chunk_idx": chunk_idx,
        "tool_key": "rfdiffusion3-design",
        "device": device,
        "seed": seed,
        "n_batches": 1,
        "diffusion_batch_size": batch_size,
        "num_timesteps": num_timesteps,
        "step_scale": step_scale,
        "gamma_0": gamma_0,
        "infer_ori_strategy": "hotspots",
        "is_non_loopy": True,
        "select_hotspots": target["select_hotspots"],
        "contigs": target["contigs"],
        "target_pdb": target["target_pdb"],
        "expected_designs": n_expected,
    }
    cdir = chunk_dir(name, chunk_idx)
    cdir.mkdir(parents=True, exist_ok=True)
    (cdir / "run_config.json").write_text(json.dumps(run_config, indent=2) + "\n")

    if dry_run:
        log(f"[{name}] dry-run — not calling Modal")
        return run_config

    t0 = time.perf_counter()
    out = run_rfdiffusion3(inputs, config)
    elapsed = time.perf_counter() - t0
    log(f"[{name}] chunk {chunk_idx} finished success={out.success} wall={elapsed:.1f}s")
    if not out.success:
        err = {"success": False, "errors": out.errors, "elapsed_s": elapsed, **run_config}
        (cdir / "error.json").write_text(json.dumps(err, indent=2) + "\n")
        raise RuntimeError(f"RFD3 failed for {name} chunk {chunk_idx}: {out.errors}")

    records = save_designs(out, name, chunk_idx, global_start_idx=existing)
    done = {
        **run_config,
        "success": True,
        "elapsed_s": elapsed,
        "n_written": len(records),
        "designs": records,
        "tool_id": getattr(out, "tool_id", None),
        "execution_time": getattr(out, "execution_time", None),
    }
    (cdir / "done.json").write_text(json.dumps(done, indent=2) + "\n")
    with open(manifest_path(name), "a") as fh:
        fh.write(json.dumps({"chunk_idx": chunk_idx, "n_written": len(records), "elapsed_s": elapsed}) + "\n")
    log(f"[{name}] wrote {len(records)} PDBs (total now {count_existing_backbones(name)})")
    return done


def chunks_needed(n: int, per_chunk: int) -> int:
    return (n + per_chunk - 1) // per_chunk


def run_target(
    target: dict[str, Any],
    n: int,
    *,
    max_chunks: int | None,
    base_seed: int,
    batch_size: int,
    step_scale: float,
    gamma_0: float,
    num_timesteps: int,
    timeout: int,
    device: str,
    dry_run: bool,
) -> None:
    name = target["name"]
    RESULTS_ROOT.joinpath(name).mkdir(parents=True, exist_ok=True)
    (RESULTS_ROOT / name / "run_config.json").write_text(
        json.dumps(
            {
                "target": target,
                "n_requested": n,
                "baker_protocol": {
                    "step_scale": step_scale,
                    "gamma_0": gamma_0,
                    "is_non_loopy": True,
                    "infer_ori_strategy": "hotspots",
                },
            },
            indent=2,
        )
        + "\n"
    )

    per_chunk = designs_per_chunk(target, batch_size)
    n_chunks = chunks_needed(n, per_chunk)
    if max_chunks is not None:
        n_chunks = min(n_chunks, max_chunks)

    log(f"\n######## {name}: need {n}, {per_chunk}/chunk → {n_chunks} chunks ########")
    for i in range(n_chunks):
        have = count_existing_backbones(name)
        if have >= n:
            log(f"[{name}] already have {have} ≥ {n}; stopping")
            break
        seed = base_seed + i * 1009 + hash(name) % 10007
        run_chunk(
            target,
            i,
            seed=seed,
            batch_size=batch_size,
            step_scale=step_scale,
            gamma_0=gamma_0,
            num_timesteps=num_timesteps,
            timeout=timeout,
            device=device,
            dry_run=dry_run,
        )

    bb = sorted(backbone_dir(name).glob("*.pdb"))
    if len(bb) > n:
        for extra in bb[n:]:
            extra.unlink()
        log(f"[{name}] trimmed to {n} backbones")
    log(f"[{name}] final count: {count_existing_backbones(name)}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", default="all", help="Target name or 'all'")
    p.add_argument("--n", type=int, default=100, help="Backbones per target")
    p.add_argument("--chunks", type=int, default=None, help="Max chunks (for control smoke)")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--step-scale", type=float, default=3.0)
    p.add_argument("--gamma-0", type=float, default=0.2)
    p.add_argument("--num-timesteps", type=int, default=200)
    p.add_argument("--timeout", type=int, default=7200)
    p.add_argument("--device", default="modal")
    p.add_argument("--seed", type=int, default=20260815)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    registry = load_registry()
    targets = registry["targets"] if args.target == "all" else [target_by_name(registry, args.target)]

    for t in targets:
        run_target(
            t,
            args.n,
            max_chunks=args.chunks,
            base_seed=args.seed,
            batch_size=args.batch_size,
            step_scale=args.step_scale,
            gamma_0=args.gamma_0,
            num_timesteps=args.num_timesteps,
            timeout=args.timeout,
            device=args.device,
            dry_run=args.dry_run,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
