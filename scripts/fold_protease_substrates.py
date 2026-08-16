"""Boltz-2 structural test of protease/substrate cleavage geometry.

Modes
-----
control
    Positive control: mature cathepsin S (UniProt P25774, 115-331) co-folded with
    the CD74 CLIP-boundary segment (P04233, 105-135). Scissile bond M120|G121.
shortlist
    Co-fold the top protease/design pairs from
    ``data/processed/immuno/pda/fold_shortlist.csv``.

Geometry metric (not eyeballing)
--------------------------------
Distance from the catalytic nucleophile (Cys25 SG in mature CatS numbering =
residue index 25 in the mature chain, which is Cys at position 25 of CTSS_MATURE)
to the scissile carbonyl carbon (P1 C), plus whether P1–P1' sit in the active-site
cleft, and interface iptm / PAE from Boltz.

Usage
-----
    uv run python scripts/fold_protease_substrates.py control
    uv run python scripts/fold_protease_substrates.py shortlist [--limit 20]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "immuno_risk" / "boltz_cleavage"
SHORTLIST = ROOT / "data" / "processed" / "immuno" / "pda" / "fold_shortlist.csv"
FASTA = ROOT / "data" / "processed" / "immuno" / "pda" / "pda_designed.fasta"
CHAINS_CSV = ROOT / "data" / "processed" / "immuno" / "pda" / "pda_designed_chains.csv"

# UniProt P25774 residues 115-331: mature catalytically active cathepsin S.
# Catalytic triad in mature numbering: Cys25, His164, Asn184.
CTSS_MATURE = (
    "LPDSVDWREKGCVTEVKYQGSCGACWAFSAVGALEAQLKLKTGKLVSLSAQNLVDCSTEKYGNKGCNGGFMTTAFQYII"
    "DNKGIDSDASYPYKAMDQKCQYDSKYRAATCSKYTELPYGREDVLKEAVANKGPVSVGVDARHPSFFLYRSGVYYEPSC"
    "TQNVNHGVLVVGYGDLNGKEYWLVKNSWGHNFGEEGYIRMARNKGNHCGIASFPSYPEI"
)
# Mature CatS Cys25 is 0-based index 24 in CTSS_MATURE (C of CGACWAF…).
CTSS_CYS25_INDEX = 24  # 0-based

# UniProt P04233 residues 105-135. CLIP C-terminal boundary M120|G121.
CD74_SEGMENT = "SKMRMATPLLMQALPMGALPQGPMQNATKYG"
CD74_SEGMENT_START = 105  # UniProt numbering of first residue
# Scissile P1 = M120 → index in segment = 120 - 105 = 15
CD74_P1_INDEX = 15

# Protease sequences keyed by site_id prefix for shortlist folding.
PROTEASE_SEQS: dict[str, tuple[str, int, str]] = {
    # site_id prefix → (sequence, nucleophile 0-based index, label)
    "cathepsin_s": (CTSS_MATURE, CTSS_CYS25_INDEX, "CTSS"),
    "cathepsin_l": (
        # Mature cathepsin L (UniProt P07711, approx 114-333) — catalytic Cys25
        "APRSVDWREKGYVTPVKNQGQCGSCWAFSATGALEGQMFRKTGRLISLSEQNLVDCSGPQGNEGCNGGLMDYAFQYVQ"
        "DNGGLDSEESYPYEATEESCKYNPKYSVANDTGFVDIPKQEKALMKAVATVGPISVAIDAGHESFLFYKEGIYFEPDCS"
        "SEDMDHGVLVVGYGFESTESDNNKYWLVKNSWGEEWGMGGYVKMAKDRRNHCGIASAASYPTV",
        24,
        "CTSL",
    ),
    "cathepsin_b": (
        # Mature cathepsin B heavy-chain region proxy (UniProt P07858) — Cys29 catalytic
        "LPASFDAREQWPQCPTIKEIRDQGSCGSCWAFGAVEAISDRICIHTNAHVSVEVSAEDLLTCCGSMCGDGCNGGYPA"
        "EAWNFWTRKGLVSGGLYESHVGCRPYSIPPCEHHVNGSRPPCTGEGDTPKCSKICEPGYSPTYKQDKHYGYNSYSVSN"
        "SEKDIMAEIYKNGPVEGAFSVYSDFLLYKSGVYQHVTGEMMGGHAIRILGWGVENGTPYWLVANSWNTDWGDNGFFKIL"
        "RGENHCGIESEIVAGIPRTDQYWE",
        29,
        "CTSB",
    ),
}


def extract_cif(structure) -> str | None:
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


def parse_cif_ca_coords(cif_text: str) -> dict[tuple[str, int], tuple[float, float, float]]:
    """Minimal mmCIF CA parser: (chain_id, seq_id) → (x, y, z).

    Also returns SG atoms under key (chain_id, seq_id, 'SG') via a parallel dict
    stored under the special chain key — see ``parse_cif_atoms``.
    """
    return {k: v for k, v in parse_cif_atoms(cif_text).items() if k[2] == "CA"}


def parse_cif_atoms(
    cif_text: str,
) -> dict[tuple[str, int, str], tuple[float, float, float]]:
    """Parse atom_site loop → (chain, seq_id, atom_name) → (x,y,z)."""
    lines = cif_text.splitlines()
    in_loop = False
    headers: list[str] = []
    atoms: dict[tuple[str, int, str], tuple[float, float, float]] = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("loop_"):
            in_loop = True
            headers = []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_"):
                headers.append(lines[i].strip())
                i += 1
            if not any("atom_site" in h for h in headers):
                in_loop = False
                continue
            # Map column indices
            def col(name: str) -> int | None:
                for hi, h in enumerate(headers):
                    if h.endswith(name) or h == name:
                        return hi
                return None

            idx_chain = col("_atom_site.auth_asym_id") or col("_atom_site.label_asym_id")
            idx_seq = col("_atom_site.auth_seq_id") or col("_atom_site.label_seq_id")
            idx_atom = col("_atom_site.label_atom_id") or col("_atom_site.auth_atom_id")
            idx_x = col("_atom_site.Cartn_x")
            idx_y = col("_atom_site.Cartn_y")
            idx_z = col("_atom_site.Cartn_z")
            if None in (idx_chain, idx_seq, idx_atom, idx_x, idx_y, idx_z):
                in_loop = False
                continue
            while i < len(lines):
                row = lines[i].strip()
                if not row or row.startswith("_") or row.startswith("loop_") or row.startswith("#"):
                    break
                if row.startswith("data_") or row.startswith(";"):
                    break
                parts = row.split()
                if len(parts) <= max(idx_chain, idx_seq, idx_atom, idx_x, idx_y, idx_z):
                    i += 1
                    continue
                try:
                    chain = parts[idx_chain]
                    seq_id = int(parts[idx_seq])
                    atom = parts[idx_atom].strip('"')
                    xyz = (float(parts[idx_x]), float(parts[idx_y]), float(parts[idx_z]))
                    atoms[(chain, seq_id, atom)] = xyz
                except (ValueError, IndexError):
                    pass
                i += 1
            in_loop = False
            continue
        i += 1
        _ = in_loop
    return atoms


def dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def geometry_report(
    cif_text: str,
    *,
    protease_chain: str = "A",
    substrate_chain: str = "B",
    nucleophile_seq_id: int,
    p1_seq_id: int,
    p1_prime_seq_id: int | None = None,
) -> dict:
    """Measure catalytic-nucleophile → scissile geometry.

    ``nucleophile_seq_id`` / ``p1_seq_id`` are 1-based auth_seq_id as written
    by Boltz (usually 1..N per chain).
    """
    atoms = parse_cif_atoms(cif_text)
    # Prefer SG for cysteine proteases; fall back to CA
    nuc = atoms.get((protease_chain, nucleophile_seq_id, "SG")) or atoms.get(
        (protease_chain, nucleophile_seq_id, "CA")
    )
    # Scissile carbonyl carbon ≈ backbone C of P1
    p1_c = atoms.get((substrate_chain, p1_seq_id, "C")) or atoms.get(
        (substrate_chain, p1_seq_id, "CA")
    )
    p1_ca = atoms.get((substrate_chain, p1_seq_id, "CA"))
    p1p_ca = None
    if p1_prime_seq_id is not None:
        p1p_ca = atoms.get((substrate_chain, p1_prime_seq_id, "CA"))

    report: dict = {
        "nucleophile_atom_found": nuc is not None,
        "p1_carbonyl_found": p1_c is not None,
        "nucleophile_to_p1_c_A": None,
        "nucleophile_to_p1_ca_A": None,
        "p1_p1prime_ca_A": None,
        "in_cleft_proxy": None,
        "caveat": (
            "Geometry from predicted coordinates; low pLDDT/iptm means the distance "
            "is not trustworthy. Cys SG→P1 C < ~4 Å with high iptm would support cleavage pose."
        ),
    }
    if nuc and p1_c:
        report["nucleophile_to_p1_c_A"] = round(dist(nuc, p1_c), 2)
    if nuc and p1_ca:
        report["nucleophile_to_p1_ca_A"] = round(dist(nuc, p1_ca), 2)
    if p1_ca and p1p_ca:
        report["p1_p1prime_ca_A"] = round(dist(p1_ca, p1p_ca), 2)
    # Crude cleft proxy: nucleophile within 8 Å of P1 CA
    if report["nucleophile_to_p1_ca_A"] is not None:
        report["in_cleft_proxy"] = report["nucleophile_to_p1_ca_A"] < 8.0
    return report


def run_boltz(chains: list[dict], *, timeout: int = 3600, use_msa: bool = False) -> tuple[object, float]:
    from proto_tools.tools.structure_prediction.boltz2.boltz2 import (
        Boltz2Config,
        Boltz2Input,
        run_boltz2,
    )

    # device lives on Boltz2Config in current proto-tools (not as a run_boltz2 kwarg).
    config = Boltz2Config(
        device="modal",
        use_msa=use_msa,
        recycling_steps=3 if use_msa else 1,
        sampling_steps=200 if use_msa else 50,
        verbose=2,
        timeout=timeout,
    )
    started = time.time()
    out = run_boltz2(Boltz2Input(complexes=[{"chains": chains}]), config)
    return out, time.time() - started


def load_design_sequences() -> dict[str, str]:
    seqs: dict[str, str] = {}
    if CHAINS_CSV.exists():
        with CHAINS_CSV.open() as f:
            for row in csv.DictReader(f):
                seqs[row["sequence_id"]] = row["sequence"]
        return seqs
    # Fallback: FASTA
    sid = None
    parts: list[str] = []
    for line in FASTA.read_text().splitlines():
        if line.startswith(">"):
            if sid:
                seqs[sid] = "".join(parts)
            sid = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if sid:
        seqs[sid] = "".join(parts)
    return seqs


def resolve_protease(site_id: str) -> tuple[str, int, str] | None:
    for prefix, val in PROTEASE_SEQS.items():
        if site_id.startswith(prefix):
            return val
    return None


def fold_control(*, use_msa: bool = False) -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    chains = [
        {"sequence": CTSS_MATURE, "entity_type": "protein"},
        {"sequence": CD74_SEGMENT, "entity_type": "protein"},
    ]
    print("[control] submitting CatS + CD74 CLIP segment to Modal…", flush=True)
    out, elapsed = run_boltz(chains, use_msa=use_msa)
    summary: dict = {
        "mode": "control",
        "elapsed_s": round(elapsed, 1),
        "success": getattr(out, "success", None),
        "errors": getattr(out, "errors", None),
        "scissile_bond": "M120|G121 (CLIP C-terminal boundary)",
        "protease": "CTSS_mature",
        "substrate": "CD74_105-135",
    }
    if not out.structures:
        summary["status"] = "SMOKE-FAIL"
        (OUT_DIR / "control_metrics.json").write_text(json.dumps(summary, indent=2))
        return summary

    struct = out.structures[0]
    m = struct.metrics
    summary.update(
        {
            "confidence_score": getattr(m, "confidence_score", None),
            "ptm": getattr(m, "ptm", None),
            "iptm": getattr(m, "iptm", None),
            "complex_plddt": getattr(m, "complex_plddt", None),
            "avg_pae": getattr(m, "avg_pae", None),
        }
    )
    cif = extract_cif(struct)
    if cif:
        cif_path = OUT_DIR / "control_model.cif"
        cif_path.write_text(cif)
        summary["cif"] = str(cif_path)
        # Boltz numbers chains A=protease, B=substrate; seq_id 1-based
        summary["geometry"] = geometry_report(
            cif,
            protease_chain="A",
            substrate_chain="B",
            nucleophile_seq_id=CTSS_CYS25_INDEX + 1,
            p1_seq_id=CD74_P1_INDEX + 1,
            p1_prime_seq_id=CD74_P1_INDEX + 2,
        )
    summary["status"] = "SMOKE-OK"
    (OUT_DIR / "control_metrics.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return summary


def fold_shortlist(*, limit: int = 20, use_msa: bool = False) -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not SHORTLIST.exists():
        raise SystemExit(f"shortlist not found: {SHORTLIST} — run join-epitope-protease first")
    designs = load_design_sequences()
    rows = list(csv.DictReader(SHORTLIST.open()))
    results = []
    n = 0
    for row in rows:
        if n >= limit:
            break
        site_id = row["site_id"]
        prot = resolve_protease(site_id)
        if prot is None:
            continue
        chain_id = row["chain_id"]
        design_seq = designs.get(chain_id)
        if not design_seq:
            continue
        # Keep substrate manageable for Boltz: ±20 aa around cleavage
        cut = int(row["cleavage_position"])
        window_start = max(0, cut - 20)
        window_end = min(len(design_seq), cut + 22)
        substrate = design_seq[window_start:window_end]
        p1_in_window = cut - window_start  # 0-based
        prot_seq, nuc_idx, prot_label = prot

        tag = f"{chain_id}__{site_id}__{n}"
        print(f"[shortlist] {tag} protease={prot_label} substrate_len={len(substrate)}", flush=True)
        chains = [
            {"sequence": prot_seq, "entity_type": "protein"},
            {"sequence": substrate, "entity_type": "protein"},
        ]
        try:
            out, elapsed = run_boltz(chains, use_msa=use_msa)
        except Exception as exc:  # noqa: BLE001
            results.append({"tag": tag, "status": "error", "error": str(exc)})
            continue

        entry: dict = {
            "tag": tag,
            "chain_id": chain_id,
            "site_id": site_id,
            "protease": prot_label,
            "peptide": row.get("peptide"),
            "allele": row.get("allele"),
            "percentile_rank": row.get("percentile_rank"),
            "relation": row.get("relation"),
            "cleavage_position": cut,
            "substrate_window": [window_start, window_end],
            "elapsed_s": round(elapsed, 1),
            "success": getattr(out, "success", None),
        }
        if out.structures:
            struct = out.structures[0]
            m = struct.metrics
            entry.update(
                {
                    "confidence_score": getattr(m, "confidence_score", None),
                    "ptm": getattr(m, "ptm", None),
                    "iptm": getattr(m, "iptm", None),
                    "complex_plddt": getattr(m, "complex_plddt", None),
                }
            )
            cif = extract_cif(struct)
            if cif:
                cif_path = OUT_DIR / f"{tag}.cif"
                cif_path.write_text(cif)
                entry["cif"] = str(cif_path)
                entry["geometry"] = geometry_report(
                    cif,
                    nucleophile_seq_id=nuc_idx + 1,
                    p1_seq_id=p1_in_window + 1,
                    p1_prime_seq_id=p1_in_window + 2,
                )
            entry["status"] = "ok"
        else:
            entry["status"] = "no_structure"
            entry["errors"] = getattr(out, "errors", None)
        results.append(entry)
        n += 1
        (OUT_DIR / "shortlist_metrics.json").write_text(json.dumps(results, indent=2))

    summary = {"mode": "shortlist", "n": len(results), "results": results}
    (OUT_DIR / "shortlist_metrics.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"n": len(results), "statuses": [r.get("status") for r in results]}, indent=2))
    return summary


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["control", "shortlist"])
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--use-msa", action="store_true")
    args = p.parse_args()
    if args.mode == "control":
        fold_control(use_msa=args.use_msa)
    else:
        fold_shortlist(limit=args.limit, use_msa=args.use_msa)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
