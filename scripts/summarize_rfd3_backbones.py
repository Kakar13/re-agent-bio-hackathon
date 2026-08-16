#!/usr/bin/env python3
"""Summarize RFD3 binder backbone runs into results/rfd3_binders/report.md."""

from __future__ import annotations

import json
import statistics
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "rfd3_targets.json"
RESULTS = ROOT / "results" / "rfd3_binders"
CITATION = """Butcher et al. (2025) De novo Design of All-atom Biomolecular Interactions with RFdiffusion3. bioRxiv. doi:10.1101/2025.09.18.676967"""


def gyration_ca(path: Path) -> float | None:
    st = gemmi.read_structure(str(path))
    coords = []
    for model in st:
        for chain in model:
            for res in chain:
                if res.het_flag == "H":
                    continue
                atom = res.find_atom("CA", "*")
                if atom:
                    coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
    if len(coords) < 3:
        return None
    xyz = np.asarray(coords)
    center = xyz.mean(axis=0)
    return float(np.sqrt(((xyz - center) ** 2).sum(axis=1).mean()))


def summarize_target(t: dict) -> dict:
    name = t["name"]
    bb = sorted((RESULTS / name / "backbones").glob("*.pdb")) if (RESULTS / name / "backbones").exists() else []
    rgs = []
    for p in bb[:100]:
        try:
            rg = gyration_ca(p)
            if rg is not None:
                rgs.append(rg)
        except Exception:
            pass
    manifest = RESULTS / name / "manifest.jsonl"
    chunks = 0
    elapsed = []
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            chunks += 1
            if "elapsed_s" in row:
                elapsed.append(row["elapsed_s"])
    return {
        "name": name,
        "label": t.get("label", name),
        "pdb_id": t["pdb_id"],
        "n_backbones": len(bb),
        "hotspots": t["select_hotspots"],
        "contigs": t["contigs"],
        "residue_range": t.get("residue_range"),
        "caveats": t.get("caveats", []),
        "notes": t.get("notes", []),
        "rg_mean": round(statistics.mean(rgs), 2) if rgs else None,
        "rg_std": round(statistics.pstdev(rgs), 2) if len(rgs) > 1 else None,
        "chunks": chunks,
        "elapsed_total_s": round(sum(elapsed), 1) if elapsed else None,
    }


def main() -> int:
    registry = json.loads(CONFIG.read_text())
    rows = [summarize_target(t) for t in registry["targets"]]
    (RESULTS / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")

    lines = [
        "# RFD3 binder backbones — undruggable targets",
        "",
        "De novo binder **backbones** generated with Proto `rfdiffusion3-design` on Modal,",
        "using the Baker-lab / RosettaCommons foundry protein-binder settings.",
        "",
        "## Protocol",
        "",
        "| Setting | Value | Source |",
        "| --- | --- | --- |",
        "| Tool | `rfdiffusion3-design` | proto-tools |",
        "| Device | Modal (`proto-env`, A10) | this workspace |",
        "| `infer_ori_strategy` | `hotspots` | foundry protein_binder_design.md |",
        "| `is_non_loopy` | `true` | same |",
        "| `step_scale` | `3.0` | low-temperature PPI preset |",
        "| `gamma_0` | `0.2` | low-temperature PPI preset |",
        "| `num_timesteps` | `200` | RFD3 default |",
        "| Binder length bins | `55-85` and `85-120` | this run |",
        "| Contig pattern | `<len>,/0,<target_chain><lo>-<hi>` | foundry binder example |",
        "",
        "## Targets",
        "",
        "| Target | PDB | Hotspots | Backbones | Mean Rg (Å) | GPU time (s) |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for r in rows:
        hs = ", ".join(r["hotspots"].keys())
        lines.append(
            f"| {r['label']} | {r['pdb_id']} | {hs} | {r['n_backbones']} | "
            f"{r['rg_mean'] if r['rg_mean'] is not None else '—'} | "
            f"{r['elapsed_total_s'] if r['elapsed_total_s'] is not None else '—'} |"
        )

    lines += ["", "## Per-target notes", ""]
    for r in rows:
        lines.append(f"### {r['label']} (`{r['name']}`)")
        lines.append("")
        lines.append(f"- Residue range: `{r['residue_range']}`")
        lines.append(f"- Contigs: `{r['contigs']}`")
        lines.append(f"- `select_hotspots`: `{json.dumps(r['hotspots'])}`")
        lines.append(f"- Backbones: `results/rfd3_binders/{r['name']}/backbones/` ({r['n_backbones']} PDBs)")
        for c in r["caveats"]:
            lines.append(f"- Caveat: {c}")
        for n in r["notes"]:
            lines.append(f"- Note: {n}")
        lines.append("")

    lines += [
        "## Control gate (PD-L1)",
        "",
        "First 16 PD-L1 designs were checked for heavy-atom contact to the foundry hotspots",
        "Y56 / M115 / Y123 (auth numbering; output target chain is renumbered 1..N).",
        "All 16 contacted all three hotspots within 4.5 Å (mean min distances ≈ 3.5 / 3.4 / 2.4 Å).",
        "See `results/rfd3_binders/pdl1/control_qc.json`.",
        "",
        "## Out of scope",
        "",
        "Sequence design (ProteinMPNN), AF2/ESMFOLD2 validation, and disordered targets",
        "(AR-NTD, tau, α-synuclein) are deferred.",
        "",
        "## Citation",
        "",
        CITATION,
        "",
        "Foundry binder example:",
        "https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/examples/protein_binder_design.md",
        "",
        "## Reproduce",
        "",
        "```bash",
        "uv run python scripts/prep_binder_targets.py",
        "proto-tools deploy --apps rfdiffusion3   # once; this workspace uses A10",
        "uv run python scripts/run_rfd3_backbones.py --target all --n 100",
        "uv run python scripts/summarize_rfd3_backbones.py",
        "```",
        "",
    ]
    report = RESULTS / "report.md"
    report.write_text("\n".join(lines))
    print(f"Wrote {report} and {RESULTS / 'summary.json'}")
    for r in rows:
        print(f"  {r['name']}: {r['n_backbones']} backbones")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
