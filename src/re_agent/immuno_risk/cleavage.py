"""Subsite-aware protease cleavage prediction for antigen-processing diagnostics.

Python is the source of truth for cohort-scale scans. The TypeScript harness
mirrors the same catalytic-site catalog for interactive demos; a parity test
pins both to identical output on a fixture sequence.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

from re_agent.immuno_risk.batch import parse_fasta
from re_agent.immuno_risk.peptides import clean_sequence
from re_agent.immuno_risk.reference_data import ROOT

DEFAULT_OUT = ROOT / "data" / "processed" / "immuno" / "pda" / "cleavage_cohort.csv"


@dataclass(frozen=True)
class CatalyticSite:
    id: str
    name: str
    protease_class: str
    motif: str
    p1: tuple[str, ...] = ()
    p2: tuple[str, ...] = ()
    p3: tuple[str, ...] = ()
    p1_prime: tuple[str, ...] = ()
    blocked_p1_prime: tuple[str, ...] = ()
    notes: str = ""
    # Special pattern handlers: "furin", "mmp", "ctsb_dipeptidyl", or "" for P1/subsite scan
    pattern: str = ""


@dataclass
class CleavageEvent:
    site_id: str
    site_name: str
    protease_class: str
    position: int  # 0-based P1 index; cut is after this residue
    p1: str
    p1_prime: str
    p2: str = ""
    p3: str = ""
    score: float = 1.0
    n_terminal_product: str = ""
    c_terminal_product: str = ""


# ---------------------------------------------------------------------------
# Catalytic-site catalog (legacy demo sites + immunologically relevant set)
# ---------------------------------------------------------------------------

DEFAULT_CATALYTIC_SITES: list[CatalyticSite] = [
    CatalyticSite(
        id="trypsin_kr",
        name="Trypsin-like (K/R)",
        protease_class="serine",
        motif="P1 = K/R; blocked if P1' = P",
        p1=("K", "R"),
        blocked_p1_prime=("P",),
        notes="Classic tryptic cut; intracellular processing proxy.",
    ),
    CatalyticSite(
        id="chymotrypsin_fwy",
        name="Chymotrypsin-like (F/Y/W)",
        protease_class="serine",
        motif="P1 = F/Y/W",
        p1=("F", "Y", "W"),
        notes="Aromatic P1 preference.",
    ),
    CatalyticSite(
        id="caspase_asp",
        name="Caspase-like (D)",
        protease_class="cysteine",
        motif="P1 = D (simplified; real caspases need DXXD-like context)",
        p1=("D",),
        notes="Stub — tighten to DXXD when evaluating apoptosis-linked paths.",
    ),
    CatalyticSite(
        id="furin_rxkr",
        name="Furin-like (R-X-K/R-R)",
        protease_class="serine",
        motif="R-X-[KR]-R↓",
        pattern="furin",
        notes="Handled by dedicated pattern matcher, not plain P1 list.",
    ),
    CatalyticSite(
        id="thrombin_r",
        name="Thrombin-like (R)",
        protease_class="serine",
        motif="P1 = R (simplified)",
        p1=("R",),
        notes="Overlaps trypsin; kept as separate labeled site for demos.",
    ),
    CatalyticSite(
        id="pepsin_fl",
        name="Pepsin-like (F/L)",
        protease_class="aspartic",
        motif="P1 = F/L (simplified)",
        p1=("F", "L"),
    ),
    CatalyticSite(
        id="elastase_agv",
        name="Elastase-like (A/G/V)",
        protease_class="serine",
        motif="P1 = A/G/V",
        p1=("A", "G", "V"),
    ),
    CatalyticSite(
        id="proteasome_hydrophobic",
        name="Proteasome-like (hydrophobic C-term preference)",
        protease_class="threonine",
        motif="P1 = L/I/V/F/Y (cytosolic MHC I antigen processing proxy)",
        p1=("L", "I", "V", "F", "Y"),
        notes="Coarse stand-in for immunoproteasome cuts feeding MHC I.",
    ),
    CatalyticSite(
        id="mmp_gp",
        name="MMP-like (G-P soft site)",
        protease_class="metallo",
        motif="P1–P1' ≈ G-P / P-X (simplified GP motif scan)",
        pattern="mmp",
        notes="Pattern: GP soft sites via custom rule.",
    ),
    CatalyticSite(
        id="legumain_n",
        name="Legumain / AEP (N)",
        protease_class="cysteine",
        motif="P1 = N; prefer non-Pro P1'",
        p1=("N",),
        blocked_p1_prime=("P",),
        notes="Asn-specific endopeptidase (AEP/legumain) — MHC-II pathway relevant.",
    ),
    # --- Extended immunologically relevant proteases ---
    CatalyticSite(
        id="cathepsin_s",
        name="Cathepsin S (P2 hydrophobic)",
        protease_class="cysteine",
        motif="P2 = V/L/I/M/F; P1 broad (not P)",
        p1=("A", "G", "S", "T", "N", "Q", "K", "R", "H", "L", "I", "V", "M", "F", "Y", "W", "E", "D"),
        p2=("V", "L", "I", "M", "F"),
        blocked_p1_prime=("P",),
        notes=(
            "Endosomal MHC-II processing. CatS is P2-driven (bulky hydrophobic); "
            "P1 is relatively permissive. Cite: Riese et al.; Turk et al. MEROPS C01.034."
        ),
    ),
    CatalyticSite(
        id="cathepsin_l",
        name="Cathepsin L (P2 aromatic/hydrophobic)",
        protease_class="cysteine",
        motif="P2 = F/Y/W/L; P1 broad",
        p1=("A", "G", "S", "T", "N", "Q", "K", "R", "H", "L", "I", "V", "M", "F", "Y", "W", "E", "D"),
        p2=("F", "Y", "W", "L"),
        blocked_p1_prime=("P",),
        notes="Lysosomal endopeptidase; aromatic/Leu P2 preference. MEROPS C01.032.",
    ),
    CatalyticSite(
        id="cathepsin_b",
        name="Cathepsin B (Arg P1 + hydrophobic P2)",
        protease_class="cysteine",
        motif="P1 = R/K; P2 = hydrophobic; also dipeptidyl-carboxypeptidase",
        p1=("R", "K"),
        p2=("V", "L", "I", "M", "F", "A"),
        blocked_p1_prime=("P",),
        notes=(
            "Endopeptidase mode: Arg/Lys P1 with hydrophobic P2. "
            "Also has dipeptidyl-carboxypeptidase activity (pattern=ctsb_dipeptidyl)."
        ),
    ),
    CatalyticSite(
        id="cathepsin_b_cpx",
        name="Cathepsin B dipeptidyl-CPX",
        protease_class="cysteine",
        motif="Removes C-terminal dipeptides (P1–P1' at C-terminus −2)",
        pattern="ctsb_dipeptidyl",
        notes="Dipeptidyl carboxypeptidase mode; cuts two residues from the C-terminus.",
    ),
    CatalyticSite(
        id="immunoproteasome_b5i",
        name="Immunoproteasome β5i / LMP7 (chymotrypsin-like)",
        protease_class="threonine",
        motif="P1 = F/Y/W/L (hydrophobic/aromatic)",
        p1=("F", "Y", "W", "L"),
        notes="IFN-γ-induced immunoproteasome subunit; MHC-I epitope C-termini.",
    ),
    CatalyticSite(
        id="immunoproteasome_b1i",
        name="Immunoproteasome β1i / LMP2 (caspase-like)",
        protease_class="threonine",
        motif="P1 = D/E (acidic)",
        p1=("D", "E"),
        notes="Immunoproteasome caspase-like activity.",
    ),
    CatalyticSite(
        id="immunoproteasome_b2i",
        name="Immunoproteasome β2i / MECL-1 (trypsin-like)",
        protease_class="threonine",
        motif="P1 = K/R (basic)",
        p1=("K", "R"),
        blocked_p1_prime=("P",),
        notes="Immunoproteasome trypsin-like activity.",
    ),
]


def get_catalytic_sites(ids: Sequence[str] | None = None) -> list[CatalyticSite]:
    if not ids:
        return list(DEFAULT_CATALYTIC_SITES)
    want = set(ids)
    return [s for s in DEFAULT_CATALYTIC_SITES if s.id in want]


def sites_as_dicts(ids: Sequence[str] | None = None) -> list[dict[str, Any]]:
    """JSON-serializable catalog for CLI / harness parity."""
    out = []
    for s in get_catalytic_sites(ids):
        d = asdict(s)
        # Mirror TypeScript camelCase keys used by the harness
        d["proteaseClass"] = d.pop("protease_class")
        d["blockedP1Prime"] = list(d.pop("blocked_p1_prime"))
        d["p1Prime"] = list(d.pop("p1_prime"))
        d["p1"] = list(d["p1"])
        d["p2"] = list(d["p2"])
        d["p3"] = list(d["p3"])
        out.append(d)
    return out


def _aa_at(seq: str, i: int) -> str:
    if i < 0 or i >= len(seq):
        return ""
    return seq[i]


def _subsite_ok(preferred: tuple[str, ...], aa: str) -> bool:
    """Empty preference = unrestricted; otherwise aa must be in the set."""
    if not preferred:
        return True
    return bool(aa) and aa in preferred


def _emit(
    out: list[CleavageEvent],
    site: CatalyticSite,
    seq: str,
    p1_index: int,
    *,
    score: float = 1.0,
) -> None:
    if p1_index < 0 or p1_index >= len(seq) - 1:
        # Allow C-terminal dipeptidyl cuts where P1' is the last residue
        if not (site.pattern == "ctsb_dipeptidyl" and p1_index == len(seq) - 2 and len(seq) >= 2):
            if p1_index < 0 or p1_index >= len(seq) - 1:
                return
    p1 = _aa_at(seq, p1_index)
    p1p = _aa_at(seq, p1_index + 1)
    if site.blocked_p1_prime and p1p in site.blocked_p1_prime:
        return
    if not _subsite_ok(site.p1_prime, p1p) and site.p1_prime:
        return
    p2 = _aa_at(seq, p1_index - 1)
    p3 = _aa_at(seq, p1_index - 2)
    if site.p2 and not _subsite_ok(site.p2, p2):
        return
    if site.p3 and not _subsite_ok(site.p3, p3):
        return
    # Score: base 1.0; boost when required P2/P3 match (already required above)
    s = score
    if site.p2 and p2 in site.p2:
        s += 0.25
    if site.p3 and p3 in site.p3:
        s += 0.1
    out.append(
        CleavageEvent(
            site_id=site.id,
            site_name=site.name,
            protease_class=site.protease_class,
            position=p1_index,
            p1=p1,
            p1_prime=p1p,
            p2=p2,
            p3=p3,
            score=round(s, 3),
            n_terminal_product=seq[: p1_index + 1],
            c_terminal_product=seq[p1_index + 1 :],
        )
    )


def _match_furin(seq: str, site: CatalyticSite, out: list[CleavageEvent]) -> None:
    # R-X-[KR]-R↓  (P1 = final R)
    for i in range(0, len(seq) - 3):
        if seq[i] == "R" and seq[i + 2] in "KR" and seq[i + 3] == "R":
            _emit(out, site, seq, i + 3)


def _match_mmp(seq: str, site: CatalyticSite, out: list[CleavageEvent]) -> None:
    for i in range(len(seq) - 1):
        if seq[i] == "G" and seq[i + 1] == "P":
            _emit(out, site, seq, i)


def _match_ctsb_dipeptidyl(seq: str, site: CatalyticSite, out: list[CleavageEvent]) -> None:
    """Cathepsin B carboxypeptidase: cut after residue len-2 (remove C-terminal dipeptide)."""
    if len(seq) < 4:
        return
    # Cut between positions n-2 and n-1 (0-based P1 = len-2)
    p1_index = len(seq) - 2
    out.append(
        CleavageEvent(
            site_id=site.id,
            site_name=site.name,
            protease_class=site.protease_class,
            position=p1_index,
            p1=seq[p1_index],
            p1_prime=seq[p1_index + 1],
            p2=_aa_at(seq, p1_index - 1),
            p3=_aa_at(seq, p1_index - 2),
            score=0.8,
            n_terminal_product=seq[: p1_index + 1],
            c_terminal_product=seq[p1_index + 1 :],
        )
    )


def predict_cleavage(
    sequence: str,
    site_ids: Sequence[str] | None = None,
) -> list[CleavageEvent]:
    """Predict cleavage events against the curated catalytic-site catalog."""
    seq = clean_sequence(sequence)
    sites = get_catalytic_sites(site_ids)
    out: list[CleavageEvent] = []

    for site in sites:
        if site.pattern == "furin":
            _match_furin(seq, site, out)
            continue
        if site.pattern == "mmp":
            _match_mmp(seq, site, out)
            continue
        if site.pattern == "ctsb_dipeptidyl":
            _match_ctsb_dipeptidyl(seq, site, out)
            continue
        if not site.p1:
            continue
        p1set = set(site.p1)
        for i in range(len(seq) - 1):
            if seq[i] in p1set:
                _emit(out, site, seq, i)

    out.sort(key=lambda e: (e.position, e.site_id))
    return out


def peptide_pool(
    sequence: str,
    cleavages: Sequence[CleavageEvent],
    max_len: int = 25,
) -> list[str]:
    """Unique peptide fragments from cleavage (including full chain if no cuts)."""
    seq = clean_sequence(sequence)
    cuts = sorted({c.position + 1 for c in cleavages})
    frags: list[str] = []
    start = 0
    for end in cuts:
        if end > start:
            frags.append(seq[start:end])
        start = end
    if start < len(seq):
        frags.append(seq[start:])
    if not frags:
        frags.append(seq)

    sources = list({*frags, seq})
    windowed: set[str] = set()
    for f in sources:
        if 8 <= len(f) <= max_len:
            windowed.add(f)
        for L in range(8, 12):
            for i in range(0, len(f) - L + 1):
                windowed.add(f[i : i + L])
        for L in range(12, 19):
            for i in range(0, len(f) - L + 1):
                windowed.add(f[i : i + L])
    return sorted(windowed, key=lambda p: (len(p), p))


CLEAVAGE_FIELDS = [
    "chain_id",
    "site_id",
    "site_name",
    "protease_class",
    "position",
    "p1",
    "p1_prime",
    "p2",
    "p3",
    "score",
]


def scan_cohort_cleavage(
    fasta: Path,
    *,
    out_csv: Path | None = None,
    site_ids: Sequence[str] | None = None,
) -> dict:
    """Stream cleavage events for every chain in ``fasta`` to CSV."""
    out_csv = out_csv or DEFAULT_OUT
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    records = parse_fasta(fasta)
    n_events = 0
    by_site: dict[str, int] = {}

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CLEAVAGE_FIELDS)
        writer.writeheader()
        for chain_id, seq in records:
            events = predict_cleavage(seq, site_ids)
            for e in events:
                writer.writerow(
                    {
                        "chain_id": chain_id,
                        "site_id": e.site_id,
                        "site_name": e.site_name,
                        "protease_class": e.protease_class,
                        "position": e.position,
                        "p1": e.p1,
                        "p1_prime": e.p1_prime,
                        "p2": e.p2,
                        "p3": e.p3,
                        "score": e.score,
                    }
                )
                n_events += 1
                by_site[e.site_id] = by_site.get(e.site_id, 0) + 1
            if n_events % 100_000 < len(events):
                f.flush()

    summary = {
        "n_chains": len(records),
        "n_events": n_events,
        "events_by_site": by_site,
        "out_csv": str(out_csv),
        "site_ids": list(site_ids) if site_ids else [s.id for s in DEFAULT_CATALYTIC_SITES],
    }
    out_csv.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
    return summary


# ---------------------------------------------------------------------------
# Join: binders × cleavage → creating / destroying
# ---------------------------------------------------------------------------

JOIN_FIELDS = [
    "chain_id",
    "peptide",
    "allele",
    "epitope_start",
    "epitope_end",
    "percentile_rank",
    "presentation_score",
    "binder",
    "site_id",
    "site_name",
    "protease_class",
    "cleavage_position",
    "cleavage_score",
    "relation",  # creating_n|creating_c|destroying|flanking
]


def join_epitope_protease(
    mhc_csv: Path,
    cleavage_csv: Path,
    *,
    out_csv: Path | None = None,
    shortlist_csv: Path | None = None,
    max_rank: float = 2.0,
    shortlist_n: int = 40,
) -> dict:
    """Classify each nearby cleavage as epitope-creating or -destroying.

    Creating N-terminus: cut after position (start - 1)  → cleavage.position == start - 1
    Creating C-terminus: cut after position (end - 1)    → cleavage.position == end - 1
    Destroying: start <= cleavage.position < end - 1 (interior cut)
    """
    out_csv = out_csv or (
        ROOT / "data" / "processed" / "immuno" / "pda" / "epitope_protease_join.csv"
    )
    shortlist_csv = shortlist_csv or (
        ROOT / "data" / "processed" / "immuno" / "pda" / "fold_shortlist.csv"
    )
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # Index cleavages by chain
    cleavages_by_chain: dict[str, list[dict]] = {}
    with cleavage_csv.open() as f:
        for row in csv.DictReader(f):
            cleavages_by_chain.setdefault(row["chain_id"], []).append(row)

    n_creating = 0
    n_destroying = 0
    n_rows = 0
    shortlist_candidates: list[dict] = []

    with mhc_csv.open() as mf, out_csv.open("w", newline="") as of:
        reader = csv.DictReader(mf)
        writer = csv.DictWriter(of, fieldnames=JOIN_FIELDS)
        writer.writeheader()
        for row in reader:
            # Prefer binders / strong ranks
            try:
                rank = float(row["percentile_rank"]) if row.get("percentile_rank") not in (None, "") else 99.0
            except ValueError:
                rank = 99.0
            binder = str(row.get("binder", "0")) in {"1", "True", "true"}
            if not binder and rank > max_rank:
                continue
            chain = row["chain_id"]
            start = int(row["start"])
            end = int(row["end"])
            events = cleavages_by_chain.get(chain, [])
            for ev in events:
                pos = int(ev["position"])
                if pos == start - 1:
                    relation = "creating_n"
                elif pos == end - 1:
                    relation = "creating_c"
                elif start <= pos < end - 1:
                    relation = "destroying"
                elif abs(pos - start) <= 2 or abs(pos - (end - 1)) <= 2:
                    relation = "flanking"
                else:
                    continue
                if relation.startswith("creating"):
                    n_creating += 1
                elif relation == "destroying":
                    n_destroying += 1
                out_row = {
                    "chain_id": chain,
                    "peptide": row["peptide"],
                    "allele": row["allele"],
                    "epitope_start": start,
                    "epitope_end": end,
                    "percentile_rank": row.get("percentile_rank"),
                    "presentation_score": row.get("presentation_score"),
                    "binder": int(binder or rank <= max_rank),
                    "site_id": ev["site_id"],
                    "site_name": ev["site_name"],
                    "protease_class": ev["protease_class"],
                    "cleavage_position": pos,
                    "cleavage_score": ev.get("score", 1.0),
                    "relation": relation,
                }
                writer.writerow(out_row)
                n_rows += 1
                if relation.startswith("creating") and rank <= max_rank:
                    # Prefer immunologically relevant proteases for folding shortlist
                    if ev["site_id"].startswith(("cathepsin", "immuno", "legumain")):
                        shortlist_candidates.append({**out_row, "_rank": rank})

    # Rank shortlist: best (lowest) percentile rank, unique (chain, site_id)
    shortlist_candidates.sort(key=lambda r: (r["_rank"], -float(r.get("cleavage_score") or 0)))
    seen: set[tuple[str, str]] = set()
    shortlist: list[dict] = []
    for c in shortlist_candidates:
        key = (c["chain_id"], c["site_id"])
        if key in seen:
            continue
        seen.add(key)
        row = {k: v for k, v in c.items() if k != "_rank"}
        shortlist.append(row)
        if len(shortlist) >= shortlist_n:
            break

    with shortlist_csv.open("w", newline="") as sf:
        if shortlist:
            w = csv.DictWriter(sf, fieldnames=JOIN_FIELDS)
            w.writeheader()
            w.writerows(shortlist)
        else:
            w = csv.DictWriter(sf, fieldnames=JOIN_FIELDS)
            w.writeheader()

    summary = {
        "n_join_rows": n_rows,
        "n_creating": n_creating,
        "n_destroying": n_destroying,
        "n_shortlist": len(shortlist),
        "out_csv": str(out_csv),
        "shortlist_csv": str(shortlist_csv),
    }
    out_csv.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
    return summary
