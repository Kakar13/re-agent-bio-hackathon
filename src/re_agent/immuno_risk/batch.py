"""Cohort-scale MHC-I scoring for PDA designed chains.

Writes flat CSV rows (not pydantic PeptideHit objects) so ~2.5M predictions
stay tractable. Backends are swappable: ``mhcflurry`` (default, no license) or
``netmhcpan`` (licensed DTU binary via NETMHCPAN_BIN).
"""

from __future__ import annotations

import csv
import logging
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Literal

from re_agent.immuno_risk.mhcflurry_backend import DEFAULT_ALLELES_I, mhcflurry_version
from re_agent.immuno_risk.peptides import clean_sequence, sliding_windows
from re_agent.immuno_risk.reference_data import ROOT

log = logging.getLogger(__name__)

Backend = Literal["mhcflurry", "netmhcpan"]

DEFAULT_OUT = ROOT / "data" / "processed" / "immuno" / "pda" / "mhc1_cohort.csv"
CHUNK_SIZE = 50_000
FIELDNAMES = [
    "chain_id",
    "start",
    "end",
    "peptide",
    "allele",
    "percentile_rank",
    "presentation_score",
    "affinity_nm",
    "processing_score",
    "el_score",
    "binder",
    "method",
]


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Return [(sequence_id, sequence), ...] from a multi-FASTA."""
    records: list[tuple[str, str]] = []
    sid: str | None = None
    parts: list[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                records.append((sid, clean_sequence("".join(parts))))
            sid = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if sid is not None:
        records.append((sid, clean_sequence("".join(parts))))
    return records


def build_peptide_index(
    records: Iterable[tuple[str, str]],
    *,
    lengths: list[int] | None = None,
) -> tuple[list[str], dict[str, list[tuple[str, int, int]]]]:
    """Deduplicate peptides cohort-wide.

    Returns
    -------
    peptides
        Unique peptide strings (stable order of first occurrence).
    index
        peptide -> [(chain_id, start, end), ...]
    """
    lengths = lengths or list(range(8, 12))
    index: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    peptides: list[str] = []
    seen: set[str] = set()
    for chain_id, seq in records:
        if len(seq) < 8:
            continue
        for start, end, pep in sliding_windows(seq, lengths):
            index[pep].append((chain_id, start, end))
            if pep not in seen:
                seen.add(pep)
                peptides.append(pep)
    return peptides, dict(index)


def _chunks(items: list[str], size: int) -> Iterator[list[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _score_mhcflurry_chunk(
    peptides: list[str],
    alleles: list[str],
    predictor=None,
) -> list[dict]:
    """Score one peptide chunk per allele; return flat prediction dicts."""
    from mhcflurry import Class1PresentationPredictor

    if predictor is None:
        predictor = Class1PresentationPredictor.load()
    version = mhcflurry_version()
    rows: list[dict] = []
    # PresentationPredictor collapses multi-allele lists to best_allele; score
    # each allele separately for a full panel.
    for allele in alleles:
        df = predictor.predict(
            peptides=peptides,
            alleles=[allele],
            verbose=0,
            include_affinity_percentile=True,
        )
        for rec in df.to_dict("records"):
            affinity = float(rec.get("affinity") or 0.0)
            presentation = float(rec.get("presentation_score") or 0.0)
            processing = float(rec.get("processing_score") or 0.0)
            pct = rec.get("affinity_percentile")
            if pct is None:
                pct = rec.get("presentation_percentile")
            percentile = float(pct) if pct is not None else None
            binder = (percentile is not None and percentile <= 2.0) or presentation >= 0.7
            rows.append(
                {
                    "peptide": str(rec["peptide"]),
                    "allele": str(rec.get("best_allele") or allele),
                    "percentile_rank": percentile,
                    "presentation_score": presentation,
                    "affinity_nm": affinity,
                    "processing_score": processing,
                    "el_score": None,
                    "binder": binder,
                    "method": f"mhcflurry_presentation:{version}",
                }
            )
    return rows


def _score_netmhcpan_chunk(
    peptides: list[str],
    alleles: list[str],
) -> list[dict]:
    """Score peptides via NetMHCpan ``-p`` peptide file mode."""
    from re_agent.immuno_risk.netmhcpan import (
        PINNED_VERSION,
        _parse_netmhcpan_table,
        netmhcpan_bin,
    )

    allele_str = ",".join(a.replace("*", "") for a in alleles)
    with tempfile.TemporaryDirectory() as tmp:
        pep_path = Path(tmp) / "peptides.txt"
        pep_path.write_text("\n".join(peptides) + "\n")
        cmd = [str(netmhcpan_bin()), "-p", str(pep_path), "-a", allele_str]
        import subprocess

        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"NetMHCpan exited {proc.returncode}: {(proc.stderr or proc.stdout)[:500]}"
            )
        hits = _parse_netmhcpan_table(proc.stdout)
    rows: list[dict] = []
    for h in hits:
        rows.append(
            {
                "peptide": h.peptide,
                "allele": h.allele,
                "percentile_rank": h.percentile_rank,
                "presentation_score": None,
                "affinity_nm": None,
                "processing_score": None,
                "el_score": h.el_score,
                "binder": h.binder,
                "method": f"netmhcpan:{PINNED_VERSION}",
            }
        )
    return rows


def _expand_and_write(
    pred_rows: list[dict],
    index: dict[str, list[tuple[str, int, int]]],
    writer: csv.DictWriter,
) -> int:
    n = 0
    for pred in pred_rows:
        pep = pred["peptide"]
        for chain_id, start, end in index.get(pep, []):
            writer.writerow(
                {
                    "chain_id": chain_id,
                    "start": start,
                    "end": end,
                    "peptide": pep,
                    "allele": pred["allele"],
                    "percentile_rank": pred["percentile_rank"],
                    "presentation_score": pred["presentation_score"],
                    "affinity_nm": pred["affinity_nm"],
                    "processing_score": pred["processing_score"],
                    "el_score": pred["el_score"],
                    "binder": int(bool(pred["binder"])),
                    "method": pred["method"],
                }
            )
            n += 1
    return n


def score_cohort(
    fasta: Path,
    *,
    backend: Backend = "mhcflurry",
    alleles: list[str] | None = None,
    out_csv: Path | None = None,
    chunk_size: int = CHUNK_SIZE,
    binder_only: bool = False,
) -> dict:
    """Score all unique 8–11mers in ``fasta`` and write ``out_csv``.

    Parameters
    ----------
    binder_only
        If True, only write rows where binder==True (keeps the CSV smaller
        for downstream joins; full predictions still run).
    """
    alleles = alleles or DEFAULT_ALLELES_I
    out_csv = out_csv or DEFAULT_OUT
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    records = parse_fasta(fasta)
    peptides, index = build_peptide_index(records)
    log.info(
        "cohort: %d chains, %d unique peptides, %d alleles, backend=%s",
        len(records),
        len(peptides),
        len(alleles),
        backend,
    )

    scorer = _score_mhcflurry_chunk if backend == "mhcflurry" else _score_netmhcpan_chunk
    mhcflurry_predictor = None
    if backend == "mhcflurry":
        from mhcflurry import Class1PresentationPredictor

        log.info("loading MHCflurry Class1PresentationPredictor…")
        mhcflurry_predictor = Class1PresentationPredictor.load()

    n_written = 0
    n_pred = 0
    n_chunks = (len(peptides) + chunk_size - 1) // chunk_size

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        for ci, chunk in enumerate(_chunks(peptides, chunk_size), start=1):
            log.info("scoring chunk %d/%d (%d peptides)", ci, n_chunks, len(chunk))
            if backend == "mhcflurry":
                pred_rows = _score_mhcflurry_chunk(chunk, alleles, predictor=mhcflurry_predictor)
            else:
                pred_rows = scorer(chunk, alleles)
            n_pred += len(pred_rows)
            if binder_only:
                pred_rows = [r for r in pred_rows if r["binder"]]
            n_written += _expand_and_write(pred_rows, index, writer)
            f.flush()

    summary = {
        "n_chains": len(records),
        "n_unique_peptides": len(peptides),
        "n_alleles": len(alleles),
        "n_predictions": n_pred,
        "n_rows_written": n_written,
        "backend": backend,
        "alleles": alleles,
        "out_csv": str(out_csv),
        "binder_only": binder_only,
    }
    summary_path = out_csv.with_suffix(".summary.json")
    import json

    summary_path.write_text(json.dumps(summary, indent=2))
    log.info("wrote %d rows → %s", n_written, out_csv)
    return summary
