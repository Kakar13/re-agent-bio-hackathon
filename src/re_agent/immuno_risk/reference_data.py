"""IEDB and HLA Ligand Atlas adapters with reproducible manifests."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import httpx

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "processed" / "immuno"
FIXTURES = ROOT / "data" / "processed" / "immuno" / "fixtures"

IEDB_MHC_EXPORT = "https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_v3.zip"
IEDB_TCELL_EXPORT = "https://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip"
# Compact IQ-API style endpoints (may require auth / rate limits) — fixtures always work offline
IEDB_SEARCH = "https://query-api.iedb.org/tcell_search"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(
    name: str,
    *,
    source_url: str,
    query: str,
    release: str,
    path: Path,
    license_note: str,
    extra: dict[str, Any] | None = None,
) -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    rows = 0
    if path.exists() and path.suffix == ".csv":
        with path.open() as f:
            rows = max(sum(1 for _ in f) - 1, 0)
    manifest = {
        "name": name,
        "source_url": source_url,
        "query": query,
        "release": release,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "license": license_note,
        "path": str(path.relative_to(ROOT)) if path.exists() and path.is_relative_to(ROOT) else str(path),
        "row_count": rows,
        "checksum_sha256": _sha256(path) if path.exists() else None,
        **(extra or {}),
    }
    out = DATA / f"{name}_manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    return out


def ensure_fixtures() -> dict[str, Path]:
    """Commit-safe tiny fixtures for offline demos and tests."""
    FIXTURES.mkdir(parents=True, exist_ok=True)
    tcell = FIXTURES / "iedb_tcell_fixture.csv"
    if not tcell.exists():
        rows = [
            {
                "peptide": "GILGFVFTL",
                "allele": "HLA-A*02:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Positive",
                "source_protein": "Influenza M1",
                "reference_id": "PMID:fixture-flu",
                "publication_year": "1991",
            },
            {
                "peptide": "NLVPMVATV",
                "allele": "HLA-A*02:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Positive",
                "source_protein": "CMV pp65",
                "reference_id": "PMID:fixture-cmv",
                "publication_year": "1995",
            },
            {
                "peptide": "AAAAAAAAA",
                "allele": "HLA-A*02:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Negative",
                "source_protein": "Synthetic",
                "reference_id": "PMID:fixture-neg1",
                "publication_year": "2020",
            },
            {
                "peptide": "GGGGGGGGG",
                "allele": "HLA-A*02:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Negative",
                "source_protein": "Synthetic",
                "reference_id": "PMID:fixture-neg2",
                "publication_year": "2021",
            },
            {
                "peptide": "KTGGPIYKR",
                "allele": "HLA-A*03:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Positive",
                "source_protein": "Influenza NP",
                "reference_id": "PMID:fixture-np",
                "publication_year": "1993",
            },
            {
                "peptide": "ILKEPVHGV",
                "allele": "HLA-A*02:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Positive",
                "source_protein": "HIV RT",
                "reference_id": "PMID:fixture-hiv",
                "publication_year": "1996",
            },
            {
                "peptide": "SLLMWITQC",
                "allele": "HLA-A*02:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Positive",
                "source_protein": "NY-ESO-1",
                "reference_id": "PMID:fixture-nyeso",
                "publication_year": "2000",
            },
            {
                "peptide": "YMDGTMSQV",
                "allele": "HLA-A*02:01",
                "assay_type": "tcell",
                "qualitative_outcome": "Positive",
                "source_protein": "Tyrosinase",
                "reference_id": "PMID:fixture-tyr",
                "publication_year": "1995",
            },
        ]
        with tcell.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    write_manifest(
        "iedb_tcell_fixture",
        source_url="fixture://local",
        query="hand-curated MHC-I T-cell assays for tests",
        release="fixture-v1",
        path=tcell,
        license_note="Synthetic fixture for CI; not redistributed IEDB bulk data",
    )

    atlas = FIXTURES / "hla_ligand_atlas_fixture.csv"
    if not atlas.exists():
        rows = [
            {"peptide": "YLLPAIVHI", "allele": "HLA-A*02:01", "tissue": "spleen", "source": "benign"},
            {"peptide": "LLFGLALAV", "allele": "HLA-A*02:01", "tissue": "liver", "source": "benign"},
            {"peptide": "KVLEYVIKV", "allele": "HLA-A*02:01", "tissue": "lung", "source": "benign"},
            {"peptide": "ALWGFFPVL", "allele": "HLA-A*02:01", "tissue": "kidney", "source": "benign"},
            {"peptide": "RMFPNAPYL", "allele": "HLA-A*02:01", "tissue": "brain", "source": "benign"},
            {"peptide": "GLCTLVAML", "allele": "HLA-A*02:01", "tissue": "spleen", "source": "benign"},
        ]
        with atlas.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    write_manifest(
        "hla_ligand_atlas_fixture",
        source_url="fixture://local",
        query="tiny benign ligand set for tolerance join demos",
        release="fixture-v1",
        path=atlas,
        license_note="Synthetic fixture approximating Atlas schema",
    )

    # NetMHCpan overlap ledger placeholder (8-mers to exclude from eval)
    ledger = FIXTURES / "netmhcpan42_overlap_ledger.txt"
    if not ledger.exists():
        # Known epitopes often in training — mark for sanity-vs-holdout separation
        ledger.write_text("\n".join(["GILGFVFT", "NLVPMVAT", "ILKEPVHG", "SLLMWITQ"]) + "\n")

    write_manifest(
        "netmhcpan42_overlap_ledger",
        source_url="https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/",
        query="8-mer prefixes for leakage checks; replace with published partitions when downloaded",
        release="fixture-partial",
        path=ledger,
        license_note="Download full NetMHCpan_train/eval partitions from DTU for production eval",
    )

    return {"tcell": tcell, "atlas": atlas, "overlap_ledger": ledger}


def load_tcell_rows(path: Path | None = None) -> list[dict[str, str]]:
    ensure_fixtures()
    path = path or FIXTURES / "iedb_tcell_fixture.csv"
    with path.open() as f:
        return list(csv.DictReader(f))


def load_atlas_peptides(path: Path | None = None) -> set[str]:
    ensure_fixtures()
    path = path or FIXTURES / "hla_ligand_atlas_fixture.csv"
    with path.open() as f:
        return {row["peptide"].upper() for row in csv.DictReader(f)}


def load_overlap_8mers(path: Path | None = None) -> set[str]:
    ensure_fixtures()
    path = path or FIXTURES / "netmhcpan42_overlap_ledger.txt"
    return {line.strip().upper() for line in path.read_text().splitlines() if line.strip()}


def fetch_iedb_tcell_sample(limit: int = 200, out: Path | None = None) -> Path:
    """Best-effort live fetch; falls back to fixtures on failure."""
    ensure_fixtures()
    out = out or (DATA / "iedb_tcell_sample.csv")
    try:
        with httpx.Client(timeout=60.0) as client:
            # IQ-API shape varies; attempt a simple GET and parse JSON if present
            r = client.get(IEDB_SEARCH, params={"limit": limit})
            if r.status_code != 200:
                raise RuntimeError(f"IEDB HTTP {r.status_code}")
            data = r.json()
            rows = data if isinstance(data, list) else data.get("data") or data.get("results") or []
            if not rows:
                raise RuntimeError("empty IEDB response")
            # Normalize loosely
            norm = []
            for item in rows[:limit]:
                peptide = str(item.get("peptide") or item.get("linear_sequence") or "").upper()
                allele = str(item.get("allele") or item.get("mhc_allele") or "")
                outcome = str(item.get("qualitative_outcome") or item.get("outcome") or "")
                if not peptide or not allele:
                    continue
                norm.append(
                    {
                        "peptide": peptide,
                        "allele": allele,
                        "assay_type": "tcell",
                        "qualitative_outcome": outcome,
                        "source_protein": str(item.get("antigen") or item.get("source_protein") or ""),
                        "reference_id": str(item.get("pubmed_id") or item.get("reference_id") or ""),
                        "publication_year": str(item.get("year") or ""),
                    }
                )
            if not norm:
                raise RuntimeError("no normalized rows")
            out.parent.mkdir(parents=True, exist_ok=True)
            with out.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(norm[0].keys()))
                w.writeheader()
                w.writerows(norm)
            write_manifest(
                "iedb_tcell_sample",
                source_url=IEDB_SEARCH,
                query=f"limit={limit}",
                release="live-sample",
                path=out,
                license_note="IEDB data — cite IEDB; respect terms of use",
            )
            return out
    except Exception as exc:  # noqa: BLE001
        log.warning("IEDB live fetch failed (%s); using fixture", exc)
        return FIXTURES / "iedb_tcell_fixture.csv"


def fetch_atlas_stub(out: Path | None = None) -> Path:
    """Atlas full tables are large; ship fixture and document real download URL."""
    ensure_fixtures()
    src = FIXTURES / "hla_ligand_atlas_fixture.csv"
    out = out or (DATA / "hla_ligand_atlas.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        out.write_text(src.read_text())
    write_manifest(
        "hla_ligand_atlas",
        source_url="https://hla-ligand-atlas.org/",
        query="fixture copy — replace with release table from HLA Ligand Atlas",
        release="fixture-or-local",
        path=out,
        license_note="HLA Ligand Atlas — check release license before redistribution",
    )
    return out
