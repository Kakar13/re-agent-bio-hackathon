"""Ingest Wells Wood Protein Design Archive sequences for MHC scoring.

PDA entries are unlabeled de novo / designed proteins. They are a **de novo
evaluation cohort**, not NetMHCpan training labels. NetMHCpan-4.2e (the pinned
comparator in this repo) needs peptide–MHC BA/EL measurements to train; we only
emit FASTA / peptide files that the licensed binary can score with `-f` / `-p`.

Source: https://pragmaticproteindesign.bio.ed.ac.uk/pda/
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Any, Iterable

import httpx

from re_agent.immuno_risk.mhcflurry_backend import DEFAULT_ALLELES_I
from re_agent.immuno_risk.peptides import clean_sequence, sliding_windows
from re_agent.immuno_risk.reference_data import ROOT, write_manifest

log = logging.getLogger(__name__)

PDA_SITE = "https://pragmaticproteindesign.bio.ed.ac.uk/pda/"
PDA_API = "https://pragmaticproteindesign.bio.ed.ac.uk/pda-api"
GITHUB_REPO = "wells-wood-research/chronowska-stam-wood-2024-protein-design-archive"
GITHUB_DATA_API = f"https://api.github.com/repos/{GITHUB_REPO}/contents/data"
GITHUB_RAW = f"https://raw.githubusercontent.com/{GITHUB_REPO}/main/data"

RAW_DIR = ROOT / "data" / "raw" / "immuno" / "pda"
OUT_DIR = ROOT / "data" / "processed" / "immuno" / "pda"

HIS_TAG_RE = re.compile(r"(?:GS|LE)?H{6,}H?$", re.I)
FASTA_WRAP = 60
MIN_CLEAN_FRACTION = 0.8
DNA_LETTERS = set("ACGTN")


def netmhcpan_allele(allele: str) -> str:
    """HLA-A*02:01 -> HLA-A02:01 (NetMHCpan -a form)."""
    return allele.replace("*", "")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _first_chain_id(chain_id: str) -> str:
    token = (chain_id or "A").split(",")[0].strip()
    return token or "A"


def _strip_his_tag(seq: str) -> tuple[str, bool]:
    cleaned = HIS_TAG_RE.sub("", seq)
    return cleaned, cleaned != seq


def _looks_like_dna(seq: str) -> bool:
    s = seq.upper()
    return len(s) >= 12 and set(s) <= DNA_LETTERS


def _usable_sequence(raw: str, *, strip_his: bool) -> tuple[str, dict[str, Any]] | None:
    text = (raw or "").replace("\n", "").strip()
    if not text:
        return None
    if strip_his:
        text, stripped = _strip_his_tag(text)
    else:
        stripped = False
    letters = sum(c.isalpha() for c in text)
    seq = clean_sequence(text)
    if len(seq) < 8:
        return None
    if _looks_like_dna(seq):
        return None
    if letters and len(seq) / letters < MIN_CLEAN_FRACTION:
        return None
    return seq, {"his_tag_stripped": stripped, "raw_length": letters, "clean_length": len(seq)}


def extract_designed_chains(
    entries: Iterable[dict[str, Any]],
    *,
    designed_only: bool = True,
    include_unknown: bool = False,
    strip_his: bool = True,
    pdb_codes: set[str] | None = None,
    limit: int | None = None,
) -> list[dict[str, Any]]:
    """Keep designed (D) chains; optionally unknown (U). Skip natural (N)."""
    allowed = {"D"}
    if include_unknown or not designed_only:
        allowed.add("U")
    if not designed_only:
        allowed.update({"N", "M", "U", "D"})

    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for entry in entries:
        pdb = str(entry.get("pdb") or "").lower()
        if not pdb:
            continue
        if pdb_codes and pdb not in pdb_codes:
            continue
        tags = [str(t) for t in (entry.get("tags") or [])]
        subtitle = str(entry.get("subtitle") or "")
        release = str(entry.get("release_date") or "")
        for chain in entry.get("chains") or []:
            ctype = str(chain.get("chain_type") or "").upper() or "U"
            if designed_only and ctype not in allowed:
                continue
            chain_id = _first_chain_id(str(chain.get("chain_id") or "A"))
            key = (pdb, chain_id)
            if key in seen:
                continue
            raw = chain.get("chain_seq_nat") or chain.get("chain_seq_fasta") or chain.get("chain_seq_unnat") or ""
            parsed = _usable_sequence(str(raw), strip_his=strip_his)
            if parsed is None:
                continue
            seq, meta = parsed
            seen.add(key)
            seq_id = f"{pdb}_{chain_id}"
            rows.append(
                {
                    "sequence_id": seq_id,
                    "pdb": pdb,
                    "chain_id": chain_id,
                    "chain_type": ctype,
                    "chain_source": chain.get("chain_source") or "",
                    "sequence": seq,
                    "length": len(seq),
                    "release_date": release,
                    "subtitle": subtitle,
                    "tags": ";".join(tags),
                    "his_tag_stripped": meta["his_tag_stripped"],
                    "seq_max_sim_natural": (entry.get("seq_max_sim_natural") or {}).get("sim"),
                    "seq_max_sim_natural_partner": (entry.get("seq_max_sim_natural") or {}).get("partner"),
                }
            )
            if limit is not None and len(rows) >= limit:
                return rows
    return rows


def write_fasta(rows: list[dict[str, Any]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    chunks: list[str] = []
    for row in rows:
        chunks.append(f">{row['sequence_id']} pdb={row['pdb']} chain={row['chain_id']} type={row['chain_type']}")
        seq = row["sequence"]
        for i in range(0, len(seq), FASTA_WRAP):
            chunks.append(seq[i : i + FASTA_WRAP])
    path.write_text("\n".join(chunks) + ("\n" if chunks else ""))
    return path


def write_index_csv(rows: list[dict[str, Any]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sequence_id",
        "pdb",
        "chain_id",
        "chain_type",
        "chain_source",
        "length",
        "release_date",
        "his_tag_stripped",
        "seq_max_sim_natural",
        "seq_max_sim_natural_partner",
        "tags",
        "subtitle",
        "sequence",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    return path


def write_mhc1_peptides(
    rows: list[dict[str, Any]],
    csv_path: Path,
    peptide_path: Path,
    *,
    lengths: list[int] | None = None,
) -> int:
    lengths = lengths or list(range(8, 12))
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    unique: list[str] = []
    seen: set[str] = set()
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["peptide", "length", "sequence_id", "pdb", "start", "end"],
        )
        w.writeheader()
        for row in rows:
            for start, end, pep in sliding_windows(row["sequence"], lengths):
                w.writerow(
                    {
                        "peptide": pep,
                        "length": len(pep),
                        "sequence_id": row["sequence_id"],
                        "pdb": row["pdb"],
                        "start": start,
                        "end": end,
                    }
                )
                if pep not in seen:
                    seen.add(pep)
                    unique.append(pep)
    peptide_path.write_text("\n".join(unique) + ("\n" if unique else ""))
    return len(unique)


def write_allele_file(path: Path, alleles: list[str] | None = None) -> Path:
    alleles = alleles or DEFAULT_ALLELES_I
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(",".join(netmhcpan_allele(a) for a in alleles) + "\n")
    return path


def latest_github_curated_name(client: httpx.Client | None = None) -> str:
    http = client or httpx.Client(timeout=60.0, follow_redirects=True)
    r = http.get(GITHUB_DATA_API, headers={"User-Agent": "re-agent-immuno-risk"})
    r.raise_for_status()
    names = [
        item["name"]
        for item in r.json()
        if item.get("name", "").endswith("_data_curated.json")
    ]
    if not names:
        raise RuntimeError("No *_data_curated.json files found on GitHub PDA data/")
    return sorted(names)[-1]


def download_github_curated(
    dest: Path | None = None,
    *,
    filename: str | None = None,
    client: httpx.Client | None = None,
) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    http = client or httpx.Client(timeout=180.0, follow_redirects=True)
    filename = filename or latest_github_curated_name(http)
    dest = dest or (RAW_DIR / filename)
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log.info("Using cached PDA dump %s (%s bytes)", dest, dest.stat().st_size)
        return dest
    url = f"{GITHUB_RAW}/{filename}"
    log.info("Downloading %s", url)
    with http.stream("GET", url, headers={"User-Agent": "re-agent-immuno-risk"}) as r:
        r.raise_for_status()
        dest.write_bytes(r.read())
    return dest


def fetch_api_entries(
    pdb_codes: list[str],
    *,
    client: httpx.Client | None = None,
) -> list[dict[str, Any]]:
    http = client or httpx.Client(timeout=60.0, follow_redirects=True)
    out: list[dict[str, Any]] = []
    for pdb in pdb_codes:
        url = f"{PDA_API}/design-details/{pdb}"
        try:
            r = http.get(url, headers={"User-Agent": "re-agent-immuno-risk"})
            r.raise_for_status()
            out.append(r.json())
        except (httpx.HTTPError, json.JSONDecodeError) as exc:
            log.warning("PDA details failed for %s: %s", pdb, exc)
    return out


def fetch_api_stubs(*, client: httpx.Client | None = None) -> list[dict[str, Any]]:
    http = client or httpx.Client(timeout=60.0, follow_redirects=True)
    r = http.get(f"{PDA_API}/all-design-stubs", headers={"User-Agent": "re-agent-immuno-risk"})
    r.raise_for_status()
    data = r.json()
    if not isinstance(data, list):
        raise RuntimeError("Unexpected PDA stubs payload")
    return data


def load_entries(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text())
    if isinstance(data, dict):
        for key in ("designs", "entries", "data"):
            if isinstance(data.get(key), list):
                data = data[key]
                break
    if not isinstance(data, list):
        raise ValueError(f"PDA JSON at {path} is not a list of designs")
    return data


def ingest_pda(
    *,
    source: str = "github",
    json_path: Path | None = None,
    out_dir: Path | None = None,
    limit: int | None = None,
    pdb_codes: list[str] | None = None,
    include_unknown: bool = False,
    strip_his: bool = True,
    write_peptides: bool = False,
    alleles: list[str] | None = None,
) -> dict[str, Any]:
    """Download / parse PDA designs and write NetMHCpan-compatible FASTA."""
    out_dir = out_dir or OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    codes = {c.lower() for c in pdb_codes} if pdb_codes else None
    source_url = PDA_SITE
    release = "local"
    raw_path: Path | None = json_path

    with httpx.Client(timeout=180.0, follow_redirects=True) as http:
        if json_path:
            entries = load_entries(json_path)
            source_url = str(json_path)
            release = json_path.name
        elif source == "api":
            stubs = fetch_api_stubs(client=http)
            wanted = [s["pdb"] for s in stubs if s.get("pdb")]
            if codes:
                wanted = [p for p in wanted if p.lower() in codes]
            if limit is not None:
                wanted = wanted[:limit]
            entries = fetch_api_entries(wanted, client=http)
            source_url = f"{PDA_API}/design-details"
            release = "pda-api"
        else:
            raw_path = download_github_curated(client=http)
            entries = load_entries(raw_path)
            source_url = f"{GITHUB_RAW}/{raw_path.name}"
            release = raw_path.name.replace("_data_curated.json", "")

    rows = extract_designed_chains(
        entries,
        include_unknown=include_unknown,
        strip_his=strip_his,
        pdb_codes=codes,
        limit=limit if source != "api" else None,
    )
    fasta = write_fasta(rows, out_dir / "pda_designed.fasta")
    index = write_index_csv(rows, out_dir / "pda_designed_chains.csv")
    allele_path = write_allele_file(out_dir / "netmhcpan_alleles.txt", alleles)
    peptide_n = 0
    peptide_files: dict[str, str] = {}
    if write_peptides:
        n = write_mhc1_peptides(
            rows,
            out_dir / "pda_mhc1_peptides.csv",
            out_dir / "netmhcpan_peptides.txt",
        )
        peptide_n = n
        peptide_files = {
            "mhc1_peptides_csv": str(out_dir / "pda_mhc1_peptides.csv"),
            "netmhcpan_peptides": str(out_dir / "netmhcpan_peptides.txt"),
        }

    allele_str = allele_path.read_text().strip()
    cmd = (
        f"netMHCpan -f {fasta} -a {allele_str} -l 8,9,10,11"
    )
    (out_dir / "netmhcpan_cmd.txt").write_text(
        "\n".join(
            [
                "# Score PDA designed chains with licensed NetMHCpan-4.2e.",
                "# Do not retrain NetMHCpan on these sequences (no BA/EL labels).",
                cmd,
                "",
            ]
        )
    )

    extra = {
        "n_entries_parsed": len(entries),
        "n_designed_chains": len(rows),
        "n_unique_peptides": peptide_n,
        "netmhcpan_cmd": cmd,
        "caveat": (
            "PDA sequences have no MHC binding or T-cell labels. Use as a de novo "
            "scoring/evaluation cohort for NetMHCpan-4.2e / MHCflurry, not as training data. "
            "Retraining NetMHCpan is out of scope for this task."
        ),
        **peptide_files,
    }
    manifest = write_manifest(
        "pda_designed",
        source_url=source_url,
        query="chain_type=D designed chains; canonical AA; optional His-tag strip",
        release=str(release),
        path=index,
        license_note="PDA / Wells Wood Research Group; cite Chronowska et al. Nat Biotechnol 2025",
        extra=extra,
    )
    checksums = {
        "fasta": _sha256_bytes(fasta.read_bytes()),
        "index": _sha256_bytes(index.read_bytes()),
    }
    summary = {
        "n_designed_chains": len(rows),
        "n_entries_parsed": len(entries),
        "fasta": str(fasta),
        "index_csv": str(index),
        "alleles": str(allele_path),
        "manifest": str(manifest),
        "netmhcpan_cmd": cmd,
        "checksums_sha256": checksums,
        "raw_json": str(raw_path) if raw_path else None,
        **extra,
    }
    (out_dir / "ingest_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
