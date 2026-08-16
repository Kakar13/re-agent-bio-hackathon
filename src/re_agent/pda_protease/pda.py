"""Protein Design Archive acquisition.

The PDA has no public REST API — the site is a catch-all single-page app — so the
curated dataset comes from the Zenodo snapshot, and the current PDB state comes
from re-running PDA's own documented selection query against RCSB.

Every PDA entry is backed by a solved structure, which is what lets the
accessibility gate use experimental coordinates instead of a prediction.
"""

from __future__ import annotations

import json
import logging
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path

import httpx

from .paths import PDA_CACHE, STRUCTURE_CACHE, UNIPROT_CACHE, ensure_caches

log = logging.getLogger(__name__)

ZENODO_RECORD = "13928951"
ZENODO_FILE = "chronowska-stam-wood-2024-protein-design-archive-1.2024.09.05.zip"
ZENODO_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD}/files/{ZENODO_FILE}/content"
DATA_JSON_SUFFIX = "Data_collection_and_processing/20240827_data.json"

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILE = "https://files.rcsb.org/download/{pdb_id}.cif"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# The accessibility gate reads B-factors, which only diffraction structures carry.
DIFFRACTION_METHODS = {"X-RAY DIFFRACTION", "NEUTRON DIFFRACTION", "ELECTRON CRYSTALLOGRAPHY"}

# Natural anchors with real-world immunogenicity expectations, plus the positive
# control. Sequences are pulled from UniProt rather than pasted, so they stay
# traceable to an accession.
NATURAL_ANCHORS = {
    "HSA": ("P02768", "Human serum albumin, clinically low immunogenicity"),
    "CD74": ("P04233", "Invariant chain; cathepsin S positive control"),
}


@dataclass
class Design:
    """One designed chain, with the structure that backs it."""

    design_id: str
    pdb_id: str
    chain_id: str
    sequence: str
    length: int
    exptl_method: str
    release_date: str
    subtitle: str
    tags: list[str] = field(default_factory=list)
    source: str = "pda"

    @property
    def has_diffraction_structure(self) -> bool:
        return self.exptl_method in DIFFRACTION_METHODS

    def to_dict(self) -> dict:
        return asdict(self)


def _client(timeout: float = 180.0) -> httpx.Client:
    return httpx.Client(timeout=timeout, follow_redirects=True)


def download_zenodo_snapshot(force: bool = False) -> Path:
    """Fetch the curated PDA archive once and cache it."""
    ensure_caches()
    dest = PDA_CACHE / "pda.zip"
    if dest.exists() and dest.stat().st_size > 1_000_000 and not force:
        log.info("PDA snapshot cached at %s", dest)
        return dest
    log.info("downloading PDA snapshot from Zenodo record %s", ZENODO_RECORD)
    with _client() as c:
        r = c.get(ZENODO_URL)
        r.raise_for_status()
        dest.write_bytes(r.content)
    log.info("wrote %s (%.1f MB)", dest, dest.stat().st_size / 1e6)
    return dest


def load_pda_records(force: bool = False) -> list[dict]:
    """Raw PDA records, one per PDB entry."""
    archive = download_zenodo_snapshot(force=force)
    with zipfile.ZipFile(archive) as z:
        name = next(n for n in z.namelist() if n.endswith(DATA_JSON_SUFFIX))
        records = json.loads(z.read(name))
    log.info("loaded %d PDA records", len(records))
    return records


def extract_designs(
    records: list[dict],
    *,
    designed_only: bool = True,
    min_length: int = 40,
    max_length: int = 400,
    require_standard_aa: bool = True,
) -> list[Design]:
    """Flatten PDA records into per-chain designs suitable for analysis.

    PDA labels chains 'D' (designed), 'U' (unknown), 'N' (natural) or 'M'. Only
    'D' is unambiguously de novo, so that is the default. Chains carrying
    non-standard residues are dropped because both the position weight matrices
    and IEDB assume the standard twenty.
    """
    designs: list[Design] = []
    for rec in records:
        pdb_id = rec.get("pdb", "").lower()
        method = (rec.get("exptl_method") or ["UNKNOWN"])[0]
        for chain in rec.get("chains", []):
            if designed_only and chain.get("chain_type") != "D":
                continue
            seq = (chain.get("chain_seq_fasta") or "").strip().upper()
            if not seq:
                continue
            if require_standard_aa and not set(seq) <= STANDARD_AA:
                continue
            if not (min_length <= len(seq) <= max_length):
                continue
            # PDA collapses identical chains into one record, so chain_id can be
            # "A,B,C". Any one of them has the same sequence and fold, so take
            # the first as the representative for structure parsing.
            raw_cid = (chain.get("chain_id") or "A").strip()
            cid = raw_cid.split(",")[0].strip() or "A"
            designs.append(
                Design(
                    design_id=f"{pdb_id}_{cid}",
                    pdb_id=pdb_id,
                    chain_id=cid,
                    sequence=seq,
                    length=len(seq),
                    exptl_method=method,
                    release_date=rec.get("release_date", ""),
                    subtitle=(rec.get("subtitle") or "")[:200],
                    tags=list(rec.get("tags") or []),
                )
            )
    log.info("extracted %d designed chains", len(designs))
    return designs


def deduplicate(designs: list[Design]) -> list[Design]:
    """Collapse exact sequence duplicates, keeping the earliest release."""
    by_seq: dict[str, Design] = {}
    for d in sorted(designs, key=lambda x: x.release_date or "9999"):
        by_seq.setdefault(d.sequence, d)
    out = list(by_seq.values())
    log.info("deduplicated %d -> %d unique sequences", len(designs), len(out))
    return out


def candidate_pool(
    records: list[dict] | None = None,
    *,
    require_diffraction: bool = True,
) -> list[Design]:
    """Designs eligible for the full pipeline including the structural gate."""
    records = records if records is not None else load_pda_records()
    designs = deduplicate(extract_designs(records))
    if require_diffraction:
        before = len(designs)
        designs = [d for d in designs if d.has_diffraction_structure]
        log.info("kept %d/%d with diffraction structures (B-factors)", len(designs), before)
    return sorted(designs, key=lambda d: d.design_id)


def rcsb_synthetic_construct_ids(limit: int | None = None) -> list[str]:
    """Re-run PDA's documented selection query against current RCSB.

    Used to report how far the Zenodo snapshot has drifted from live PDB, not to
    replace the curated set — PDA's manual exclusions are the valuable part.
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                        "operator": "exact_match",
                        "value": "synthetic construct",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater",
                        "value": 0,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit or 10000}},
    }
    with _client(timeout=120.0) as c:
        r = c.post(RCSB_SEARCH, json=query)
        r.raise_for_status()
        data = r.json()
    ids = [item["identifier"].lower() for item in data.get("result_set", [])]
    log.info("RCSB reports %d synthetic-construct protein entries", data.get("total_count", 0))
    return ids


def fetch_structure(pdb_id: str, *, force: bool = False) -> Path | None:
    """Download and cache an mmCIF from RCSB."""
    ensure_caches()
    pdb_id = pdb_id.lower()
    dest = STRUCTURE_CACHE / f"{pdb_id}.cif"
    if dest.exists() and dest.stat().st_size > 1000 and not force:
        return dest
    try:
        with _client(timeout=120.0) as c:
            r = c.get(RCSB_FILE.format(pdb_id=pdb_id))
            r.raise_for_status()
            dest.write_bytes(r.content)
        return dest
    except Exception as exc:  # noqa: BLE001 - a missing structure is not fatal
        log.warning("could not fetch structure %s: %s", pdb_id, exc)
        return None


def fetch_structures(designs: list[Design], *, max_workers: int = 8) -> dict[str, Path]:
    """Fetch structures for a set of designs, skipping failures."""
    from concurrent.futures import ThreadPoolExecutor

    unique = sorted({d.pdb_id for d in designs})
    out: dict[str, Path] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        for pdb_id, path in zip(unique, pool.map(fetch_structure, unique)):
            if path is not None:
                out[pdb_id] = path
    log.info("fetched %d/%d structures", len(out), len(unique))
    return out


def fetch_uniprot_sequence(acc: str, *, force: bool = False) -> str:
    """Cached UniProt sequence lookup."""
    ensure_caches()
    dest = UNIPROT_CACHE / f"{acc}.fasta"
    if not dest.exists() or force:
        with _client(timeout=60.0) as c:
            r = c.get(UNIPROT_FASTA.format(acc=acc))
            r.raise_for_status()
            dest.write_text(r.text)
    lines = dest.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def natural_anchor_designs() -> list[Design]:
    """Non-PDA reference proteins carried through the pipeline alongside designs."""
    out = []
    for name, (acc, note) in NATURAL_ANCHORS.items():
        try:
            seq = fetch_uniprot_sequence(acc)
        except Exception as exc:  # noqa: BLE001
            log.warning("could not fetch anchor %s (%s): %s", name, acc, exc)
            continue
        out.append(
            Design(
                design_id=f"anchor_{name}",
                pdb_id="",
                chain_id="A",
                sequence=seq,
                length=len(seq),
                exptl_method="UNIPROT",
                release_date="",
                subtitle=note,
                tags=["anchor"],
                source="uniprot",
            )
        )
    return out


def write_pool(designs: list[Design], dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps([d.to_dict() for d in designs], indent=2))
    log.info("wrote %d designs to %s", len(designs), dest)
    return dest


def read_pool(src: Path) -> list[Design]:
    return [Design(**rec) for rec in json.loads(src.read_text())]
