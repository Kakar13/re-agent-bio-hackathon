#!/usr/bin/env python3
"""Prepare undruggable binder targets and RFD3 hotspot specs.

Fetches PDB/assembly structures, crops to target(+partner) chains, computes
interface residues within 4.5 A, crops contigs to continuous residue segments
(RFD3 rejects gaps), and emits:
  - data/targets/<name>_target.pdb
  - data/targets/<name>_complex.pdb
  - data/targets/<name>_hotspots_report.json
  - configs/rfd3_targets.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import urllib.request
from pathlib import Path
from typing import Any

import gemmi

ROOT = Path(__file__).resolve().parents[1]
TARGETS_DIR = ROOT / "data" / "targets"
CONFIG_PATH = ROOT / "configs" / "rfd3_targets.json"

TIP_ATOMS: dict[str, str] = {
    "ALA": "CB",
    "ARG": "NH1,NH2",
    "ASN": "OD1,ND2",
    "ASP": "OD1,OD2",
    "CYS": "SG",
    "GLN": "OE1,NE2",
    "GLU": "OE1,OE2",
    "GLY": "CA",
    "HIS": "ND1,NE2",
    "ILE": "CD1",
    "LEU": "CD1,CD2",
    "LYS": "NZ",
    "MET": "CG,SD",
    "PHE": "CD2,CZ",
    "PRO": "CG",
    "SER": "OG",
    "THR": "OG1",
    "TRP": "CD2,CH2",
    "TYR": "CD2,OH",
    "VAL": "CG1,CG2",
}

LITERATURE_HOTSPOTS: dict[str, list[int]] = {
    "kras": [37, 38, 40],
    "stat3": [609, 611, 613],
    "beta_catenin": [312, 435],
}

PDL1_DOC_HOTSPOTS: dict[str, str] = {
    "A56": "CG,OH",
    "A115": "CG,SD",
    "A123": "CD2,OH",
}


def log(msg: str) -> None:
    print(msg, flush=True)


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not (dest.exists() and dest.stat().st_size > 0):
        log(f"  downloading {url}")
        with urllib.request.urlopen(url, timeout=120) as resp, open(dest, "wb") as out:
            shutil.copyfileobj(resp, out)
    if dest.suffix == ".gz" or url.endswith(".gz"):
        unzipped = dest.with_suffix("") if dest.suffix == ".gz" else Path(str(dest) + ".out")
        if not (unzipped.exists() and unzipped.stat().st_size > 0):
            with gzip.open(dest, "rb") as fin, open(unzipped, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        return unzipped
    return dest


def fetch_entry_pdb(pdb_id: str) -> Path:
    dest = TARGETS_DIR / "raw" / f"{pdb_id.upper()}.pdb"
    return download(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", dest)


def fetch_assembly_pdb(pdb_id: str, assembly: int = 1) -> Path:
    gz = TARGETS_DIR / "raw" / f"{pdb_id.upper()}.pdb{assembly}.gz"
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb{assembly}.gz"
    try:
        raw = download(url, gz)
    except Exception:
        dest = TARGETS_DIR / "raw" / f"{pdb_id.upper()}.pdb{assembly}"
        raw = download(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb{assembly}", dest)
    normalized = TARGETS_DIR / "raw" / f"{pdb_id.upper()}_assembly{assembly}.pdb"
    if not (normalized.exists() and normalized.stat().st_size > 0):
        shutil.copy(raw, normalized)
    return normalized


def gemmi_read(path: Path) -> gemmi.Structure:
    if path.suffix.lower() not in {".pdb", ".cif", ".mmcif", ".ent"}:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp.write(path.read_bytes())
            tmp_path = Path(tmp.name)
        try:
            return gemmi.read_structure(str(tmp_path))
        finally:
            tmp_path.unlink(missing_ok=True)
    return gemmi.read_structure(str(path))


def crop_structure(
    st: gemmi.Structure,
    keep_chains: list[str],
    residue_ranges: dict[str, tuple[int, int] | None] | None = None,
    strip_hetero: bool = True,
) -> gemmi.Structure:
    residue_ranges = residue_ranges or {}
    out = gemmi.Structure()
    out.name = st.name
    out.cell = st.cell
    out.spacegroup_hm = st.spacegroup_hm
    model = gemmi.Model("1")
    src = st[0]
    for chain_id in keep_chains:
        if chain_id not in [c.name for c in src]:
            raise ValueError(f"chain {chain_id!r} not in structure (have {[c.name for c in src]})")
        new_chain = gemmi.Chain(chain_id)
        rng = residue_ranges.get(chain_id)
        for res in src[chain_id]:
            if strip_hetero and res.het_flag == "H":
                continue
            if res.name not in TIP_ATOMS and res.het_flag == "H":
                continue
            if res.name not in TIP_ATOMS:
                continue
            seqnum = res.seqid.num
            if rng is not None:
                lo, hi = rng
                if seqnum < lo or seqnum > hi:
                    continue
            new_chain.add_residue(res)
        if len(new_chain) == 0:
            raise ValueError(f"no residues kept for chain {chain_id}")
        model.add_chain(new_chain)
    out.add_model(model)
    return out


def write_gemmi_pdb(st: gemmi.Structure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    st.write_minimal_pdb(str(path))


def residue_name_map(st: gemmi.Structure, chain_id: str) -> dict[int, str]:
    return {res.seqid.num: res.name for res in st[0][chain_id] if res.het_flag != "H"}


def tip_atoms_for(resname: str, available: set[str] | None = None) -> str:
    raw = TIP_ATOMS.get(resname, "CA")
    if available is None:
        return raw
    kept = [a for a in raw.split(",") if a in available]
    if kept:
        return ",".join(kept)
    for cand in ("CA", "CB", "N", "C"):
        if cand in available:
            return cand
    return "CA"


def atom_names_in_residue(st: gemmi.Structure, chain_id: str, resnum: int) -> set[str]:
    for res in st[0][chain_id]:
        if res.seqid.num == resnum and res.het_flag != "H":
            return {a.name for a in res}
    return set()


def count_interface_contacts(
    st: gemmi.Structure,
    target_chain: str,
    partner_chain: str,
    cutoff: float = 4.5,
) -> list[dict[str, Any]]:
    prt = st[0][partner_chain]
    partner_atoms = []
    for res in prt:
        if res.het_flag == "H":
            continue
        for atom in res:
            if atom.is_hydrogen():
                continue
            partner_atoms.append(atom.pos)
    if not partner_atoms:
        return []
    rows: list[dict[str, Any]] = []
    for res in st[0][target_chain]:
        if res.het_flag == "H":
            continue
        contacts = 0
        min_dist = 999.0
        for atom in res:
            if atom.is_hydrogen():
                continue
            for p in partner_atoms:
                d = atom.pos.dist(p)
                if d < min_dist:
                    min_dist = d
                if d <= cutoff:
                    contacts += 1
        if contacts > 0:
            rows.append(
                {
                    "chain": target_chain,
                    "resnum": res.seqid.num,
                    "resname": res.name,
                    "contacts": contacts,
                    "min_dist": round(min_dist, 3),
                    "key": f"{target_chain}{res.seqid.num}",
                }
            )
    rows.sort(key=lambda r: (-r["contacts"], r["min_dist"], r["resnum"]))
    return rows


def pick_hotspots(
    st: gemmi.Structure,
    target_chain: str,
    ranked: list[dict[str, Any]],
    n: int = 3,
    preferred: list[int] | None = None,
    forced: dict[str, str] | None = None,
) -> tuple[dict[str, str], list[str]]:
    notes: list[str] = []
    if forced:
        return forced, ["using forced/doc hotspots"]
    chosen_nums: list[int] = []
    if preferred:
        avail = {r["resnum"] for r in ranked}
        hit = [p for p in preferred if p in avail]
        miss = [p for p in preferred if p not in avail]
        if hit:
            notes.append(f"literature residues present at interface: {hit}")
            chosen_nums.extend(hit[:n])
        if miss:
            notes.append(f"literature residues NOT in computed interface: {miss}")
    for row in ranked:
        if len(chosen_nums) >= n:
            break
        if row["resnum"] not in chosen_nums:
            chosen_nums.append(row["resnum"])
    hotspots: dict[str, str] = {}
    name_map = residue_name_map(st, target_chain)
    for num in chosen_nums:
        resname = name_map.get(num, "UNK")
        atoms = tip_atoms_for(resname, atom_names_in_residue(st, target_chain, num))
        hotspots[f"{target_chain}{num}"] = atoms
    return hotspots, notes


def continuous_segments(nums: list[int]) -> list[tuple[int, int]]:
    if not nums:
        return []
    nums = sorted(nums)
    segs: list[tuple[int, int]] = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
        else:
            segs.append((start, prev))
            start = prev = n
    segs.append((start, prev))
    return segs


def hotspot_covering_segment(
    nums: list[int], hotspot_resnums: list[int]
) -> tuple[int, int] | None:
    segs = continuous_segments(nums)
    if not segs:
        return None
    best = None
    best_score = (-1, -1)
    for lo, hi in segs:
        covered = sum(1 for h in hotspot_resnums if lo <= h <= hi)
        score = (covered, hi - lo + 1)
        if score > best_score:
            best_score = score
            best = (lo, hi)
    return best


def observed_range(st: gemmi.Structure, chain_id: str) -> tuple[int, int]:
    nums = [res.seqid.num for res in st[0][chain_id] if res.het_flag != "H"]
    return min(nums), max(nums)


def observed_resnums(st: gemmi.Structure, chain_id: str) -> list[int]:
    return sorted({res.seqid.num for res in st[0][chain_id] if res.het_flag != "H"})


def build_contig(binder_lo: int, binder_hi: int, target_chain: str, res_lo: int, res_hi: int) -> str:
    return f"{binder_lo}-{binder_hi},/0,{target_chain}{res_lo}-{res_hi}"


def merge_multimodel_dimer(raw: Path, out_path: Path) -> Path | None:
    st = gemmi_read(raw)
    if len(st) < 2:
        return None
    polymers_per_model = []
    for mi, model in enumerate(st):
        for chain in model:
            n_aa = sum(1 for r in chain if r.het_flag != "H" and r.name in TIP_ATOMS)
            if n_aa > 50:
                polymers_per_model.append((mi, chain.name, n_aa))
    if len(polymers_per_model) < 2:
        log(f"  multi-model dimer not found: {polymers_per_model}")
        return None
    log(f"  merging multi-model dimer: {polymers_per_model[:2]}")
    (m1, c1, _), (m2, c2, _) = polymers_per_model[0], polymers_per_model[1]
    new = gemmi.Structure()
    new.name = st.name
    new.cell = st.cell
    new.spacegroup_hm = st.spacegroup_hm
    model = gemmi.Model("1")
    ch_a = gemmi.Chain("A")
    for res in st[m1][c1]:
        if res.het_flag == "H" or res.name not in TIP_ATOMS:
            continue
        ch_a.add_residue(res)
    ch_b = gemmi.Chain("B")
    for res in st[m2][c2]:
        if res.het_flag == "H" or res.name not in TIP_ATOMS:
            continue
        ch_b.add_residue(res)
    model.add_chain(ch_a)
    model.add_chain(ch_b)
    new.add_model(model)
    write_gemmi_pdb(new, out_path)
    return out_path


def prepare_target(spec: dict[str, Any]) -> dict[str, Any]:
    name = spec["name"]
    log(f"\n=== {name} ===")
    if spec.get("assembly"):
        raw = fetch_assembly_pdb(spec["pdb_id"], spec["assembly"])
    else:
        raw = fetch_entry_pdb(spec["pdb_id"])

    st = gemmi_read(raw)
    log(f"  chains in file: {[c.name for c in st[0]]}")
    keep = list(spec["keep_chains"])
    ranges = spec.get("residue_ranges") or {}
    cropped = crop_structure(st, keep, ranges, strip_hetero=True)

    target_chain = spec["target_chain"]
    partner_chain = spec.get("partner_chain")

    complex_path = TARGETS_DIR / f"{name}_complex.pdb"
    target_path = TARGETS_DIR / f"{name}_target.pdb"
    write_gemmi_pdb(cropped, complex_path)
    target_only = crop_structure(cropped, [target_chain], strip_hetero=True)
    write_gemmi_pdb(target_only, target_path)

    lo, hi = observed_range(target_only, target_chain)
    log(f"  target {target_chain} residues {lo}-{hi} ({hi - lo + 1} span, may have gaps)")

    ranked: list[dict[str, Any]] = []
    notes: list[str] = []
    if partner_chain and partner_chain in [c.name for c in cropped[0]]:
        ranked = count_interface_contacts(cropped, target_chain, partner_chain, cutoff=4.5)
        log(f"  interface residues (4.5 A): {len(ranked)}")
    else:
        notes.append("no partner chain — using forced/doc hotspots only")

    forced = spec.get("forced_hotspots")
    preferred = LITERATURE_HOTSPOTS.get(name)
    hotspots, pick_notes = pick_hotspots(
        cropped if partner_chain else target_only,
        target_chain,
        ranked,
        n=spec.get("n_hotspots", 3),
        preferred=preferred,
        forced=forced,
    )
    notes.extend(pick_notes)
    log(f"  select_hotspots: {hotspots}")

    hs_nums = [int("".join(ch for ch in k if ch.isdigit())) for k in hotspots]
    nums = observed_resnums(target_only, target_chain)
    seg = hotspot_covering_segment(nums, hs_nums)
    if seg is None:
        raise RuntimeError(f"{name}: no continuous segment for hotspots {hs_nums}")
    seg_lo, seg_hi = seg
    if (seg_lo, seg_hi) != (lo, hi) or len(nums) != (hi - lo + 1):
        notes.append(
            f"contig cropped to continuous segment {target_chain}{seg_lo}-{seg_hi} "
            f"(observed {len(nums)} residues in {lo}-{hi}; gaps omitted)"
        )
        log(f"  continuous contig window: {target_chain}{seg_lo}-{seg_hi}")
        target_only = crop_structure(target_only, [target_chain], {target_chain: (seg_lo, seg_hi)})
        write_gemmi_pdb(target_only, target_path)
        if partner_chain:
            cropped = crop_structure(
                cropped,
                [target_chain, partner_chain],
                {target_chain: (seg_lo, seg_hi), partner_chain: None},
            )
            write_gemmi_pdb(cropped, complex_path)
        lo, hi = seg_lo, seg_hi

    binder_bins = spec.get("binder_length_bins", [[55, 85], [85, 120]])
    contigs = {f"len_{a}_{b}": build_contig(a, b, target_chain, lo, hi) for a, b in binder_bins}

    report = {
        "name": name,
        "pdb_id": spec["pdb_id"],
        "assembly": spec.get("assembly"),
        "why_undruggable": spec.get("why_undruggable"),
        "caveats": spec.get("caveats", []),
        "target_chain": target_chain,
        "partner_chain": partner_chain,
        "residue_range": [lo, hi],
        "target_pdb": str(target_path.relative_to(ROOT)),
        "complex_pdb": str(complex_path.relative_to(ROOT)),
        "select_hotspots": hotspots,
        "contigs": contigs,
        "binder_length_bins": binder_bins,
        "notes": notes,
        "top_interface_residues": ranked[:25],
        "literature_hotspots": preferred,
    }
    report_path = TARGETS_DIR / f"{name}_hotspots_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    log(f"  wrote {report_path.relative_to(ROOT)}")

    return {
        "name": name,
        "pdb_id": spec["pdb_id"],
        "label": spec.get("label", name),
        "why_undruggable": spec.get("why_undruggable"),
        "caveats": spec.get("caveats", []),
        "target_chain": target_chain,
        "partner_chain": partner_chain,
        "residue_range": [lo, hi],
        "target_pdb": str(target_path.relative_to(ROOT)),
        "complex_pdb": str(complex_path.relative_to(ROOT)),
        "select_hotspots": hotspots,
        "infer_ori_strategy": "hotspots",
        "is_non_loopy": True,
        "binder_length_bins": binder_bins,
        "contigs": contigs,
        "hotspots_report": str(report_path.relative_to(ROOT)),
        "rfd3_sampler": {
            "step_scale": 3.0,
            "gamma_0": 0.2,
            "num_timesteps": 200,
            "diffusion_batch_size": 8,
        },
        "notes": notes,
    }


def target_specs() -> list[dict[str, Any]]:
    return [
        {
            "name": "pdl1",
            "label": "PD-L1 (positive control)",
            "pdb_id": "5O45",
            "keep_chains": ["A"],
            "target_chain": "A",
            "partner_chain": None,
            "residue_ranges": {"A": (17, 131)},
            "forced_hotspots": PDL1_DOC_HOTSPOTS,
            "why_undruggable": "Validated PPI control from RFD3 protein_binder_design example",
            "caveats": ["Uses foundry doc hotspots A56/A115/A123 verbatim"],
        },
        {
            "name": "kras",
            "label": "KRAS4b (effector face)",
            "pdb_id": "6VJJ",
            "keep_chains": ["A", "B"],
            "target_chain": "A",
            "partner_chain": "B",
            "residue_ranges": {"A": (1, 169), "B": None},
            "why_undruggable": "Flat shallow effector surface; picomolar GTP affinity",
            "caveats": [
                "6VJJ is wild-type KRAS4b–RAF1 RBD, not G12D; residue 12 is not part of the binder epitope"
            ],
            "n_hotspots": 3,
        },
        {
            "name": "myc_max",
            "label": "MYC leucine zipper (vs MAX)",
            "pdb_id": "1NKP",
            "keep_chains": ["A", "B"],
            "target_chain": "A",
            "partner_chain": "B",
            "residue_ranges": None,
            "why_undruggable": "Flexible shallow bHLH-LZ interface; only Omomyc reached clinic",
            "caveats": ["Cropped to MYC chain A + MAX chain B; DNA and duplicate copies removed"],
            "n_hotspots": 3,
        },
        {
            "name": "stat3",
            "label": "STAT3 SH2 dimer interface",
            "pdb_id": "1BG1",
            "assembly": 1,
            "keep_chains": ["A", "B"],
            "target_chain": "A",
            "partner_chain": "B",
            "residue_ranges": None,
            "why_undruggable": "SH2/pTyr-mediated dimerization rather than a classical pocket",
            "caveats": [
                "Uses biological assembly 1 (two models merged to A/B) to recover the SH2 dimer",
                "Forced literature SH2 pocket hotspots R609/S611/S613 — computed contacts favor the reciprocal pY705 face",
            ],
            "forced_hotspots": {"A609": "NH1,NH2", "A611": "OG", "A613": "OG"},
            "n_hotspots": 3,
        },
        {
            "name": "beta_catenin",
            "label": "beta-catenin ARM groove (TCF site)",
            "pdb_id": "1JDH",
            "keep_chains": ["A", "B"],
            "target_chain": "A",
            "partner_chain": "B",
            "residue_ranges": {"A": (134, 664), "B": None},
            "why_undruggable": "Extended shallow armadillo groove used by TCF/LEF/APC/cadherin",
            "caveats": ["Cropped ARM repeats; hotspots from beta-catenin–hTcf-4 contacts"],
            "n_hotspots": 3,
        },
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="*", help="Subset of target names to prepare")
    args = parser.parse_args()

    TARGETS_DIR.mkdir(parents=True, exist_ok=True)
    (TARGETS_DIR / "raw").mkdir(parents=True, exist_ok=True)
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)

    specs = target_specs()
    if args.only:
        specs = [s for s in specs if s["name"] in args.only]

    for s in specs:
        if s["name"] == "stat3":
            raw = fetch_assembly_pdb(s["pdb_id"], s.get("assembly", 1))
            remapped = TARGETS_DIR / "raw" / "1BG1_dimer_AB.pdb"
            merged = merge_multimodel_dimer(raw, remapped)
            if merged is None:
                raise SystemExit("STAT3 dimer merge failed — check assembly file")
            shutil.copy(merged, TARGETS_DIR / "raw" / "1BG1.pdb")
            s.pop("assembly", None)
            s["keep_chains"] = ["A", "B"]
            s["target_chain"] = "A"
            s["partner_chain"] = "B"
            s["residue_ranges"] = {"A": (500, 716), "B": (500, 716)}

    registry = [prepare_target(spec) for spec in specs]
    payload = {
        "tool_key": "rfdiffusion3-design",
        "citation": "proto-tools citation rfdiffusion3-design",
        "n_backbones_per_target": 100,
        "baker_protocol": {
            "infer_ori_strategy": "hotspots",
            "is_non_loopy": True,
            "step_scale": 3.0,
            "gamma_0": 0.2,
            "num_timesteps": 200,
            "source": "https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/examples/protein_binder_design.md",
        },
        "targets": registry,
    }
    CONFIG_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    log(f"\nWrote {CONFIG_PATH.relative_to(ROOT)} with {len(registry)} targets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
