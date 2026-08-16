"""MHCflurry MHC-I presentation backend (required open baseline)."""

from __future__ import annotations

import importlib.metadata
import logging
from typing import Any

from re_agent.immuno_risk.peptides import clean_sequence, sliding_windows
from re_agent.immuno_risk.schemas import PeptideHit

log = logging.getLogger(__name__)

DEFAULT_ALLELES_I = [
    "HLA-A*02:01",
    "HLA-A*01:01",
    "HLA-A*03:01",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-C*07:01",
]


def mhcflurry_version() -> str:
    try:
        return importlib.metadata.version("mhcflurry")
    except importlib.metadata.PackageNotFoundError:
        return "not_installed"


def _presentation_predictor():
    from mhcflurry import Class1PresentationPredictor

    return Class1PresentationPredictor.load()


def score_mhc_i(
    sequence: str,
    alleles: list[str] | None = None,
    *,
    peptide_lengths: list[int] | None = None,
    top_n: int = 50,
) -> list[PeptideHit]:
    """Score MHC-I presentation with MHCflurry. Fails closed if unavailable."""
    alleles = alleles or DEFAULT_ALLELES_I
    lengths = peptide_lengths or list(range(8, 12))
    seq = clean_sequence(sequence)
    if len(seq) < 8:
        raise ValueError("sequence too short for MHC-I windows")

    windows = sliding_windows(seq, lengths)
    try:
        predictor = _presentation_predictor()
    except Exception as exc:  # noqa: BLE001 — fail closed unless demo heuristic allowed
        import os

        if os.environ.get("IMMUNO_ALLOW_HEURISTIC_MHC") == "1":
            log.warning("MHCflurry unavailable (%s); IMMUNO_ALLOW_HEURISTIC_MHC=1 fallback", exc)
            return _heuristic_mhc_i(windows, alleles, top_n)
        raise RuntimeError(
            "MHCflurry presentation models unavailable. "
            "Install: uv sync --extra immuno && mhcflurry-downloads fetch models_class1_presentation "
            "(or set IMMUNO_ALLOW_HEURISTIC_MHC=1 for offline demo only)"
        ) from exc

    peptides = [w[2] for w in windows]
    starts = {w[2]: w[0] for w in windows}  # first occurrence for provenance

    # MHCflurry Class1PresentationPredictor returns best_allele per peptide when
    # given a multi-allele list. Score each allele separately for a full panel.
    version = mhcflurry_version()
    hits: list[PeptideHit] = []
    for allele in alleles:
        df = predictor.predict(
            peptides=peptides,
            alleles=[allele],
            verbose=0,
            include_affinity_percentile=True,
        )
        for rec in df.to_dict("records"):
            peptide = str(rec["peptide"])
            affinity = float(rec.get("affinity") or 0.0)
            presentation = float(rec.get("presentation_score") or 0.0)
            processing = float(rec.get("processing_score") or 0.0)
            pct = rec.get("affinity_percentile")
            if pct is None:
                pct = rec.get("presentation_percentile")
            percentile = float(pct) if pct is not None else None
            binder = (percentile is not None and percentile <= 2.0) or presentation >= 0.7
            start = starts.get(peptide)
            hits.append(
                PeptideHit(
                    peptide=peptide,
                    allele=str(rec.get("best_allele") or allele),
                    mhc_class="I",
                    start=start,
                    end=(start + len(peptide)) if start is not None else None,
                    length=len(peptide),
                    affinity_nm=affinity,
                    presentation_score=presentation,
                    processing_score=processing,
                    percentile_rank=percentile,
                    binder=binder,
                    method="mhcflurry_presentation",
                    version=version,
                    provenance={"predictor": "mhcflurry", "models": "models_class1_presentation"},
                    caveat="Presentation propensity — not immunogenicity probability.",
                )
            )

    # Sort: lower percentile better; fallback higher presentation
    def sort_key(h: PeptideHit) -> tuple[float, float]:
        pct = h.percentile_rank if h.percentile_rank is not None else 50.0
        pres = -(h.presentation_score or 0.0)
        return (pct, pres)

    hits.sort(key=sort_key)
    return hits[:top_n]


def _heuristic_mhc_i(
    windows: list[tuple[int, int, str]],
    alleles: list[str],
    top_n: int,
) -> list[PeptideHit]:
    hits: list[PeptideHit] = []
    for start, end, peptide in windows:
        for allele in alleles:
            h = sum((i + 1) * ord(c) for i, c in enumerate(peptide + allele)) % 5000
            rank = round(h / 100.0, 2)
            hits.append(
                PeptideHit(
                    peptide=peptide,
                    allele=allele,
                    mhc_class="I",
                    start=start,
                    end=end,
                    length=len(peptide),
                    percentile_rank=rank,
                    presentation_score=max(0.0, 1.0 - rank / 50.0),
                    processing_score=0.5,
                    affinity_nm=500.0 * (1.0 + rank),
                    binder=rank <= 2.0,
                    method="mhc_i_heuristic_demo",
                    version="demo",
                    provenance={"predictor": "heuristic", "mhcflurry": False},
                    caveat="DEMO ONLY — MHCflurry not installed. Not a binder predictor.",
                )
            )
    hits.sort(key=lambda h: h.percentile_rank if h.percentile_rank is not None else 99.0)
    return hits[:top_n]


def score_peptide_allele_pairs(
    pairs: list[tuple[str, str]],
) -> list[dict[str, Any]]:
    """Score explicit (peptide, allele) pairs — used for IEDB feature extraction."""
    if not pairs:
        return []
    try:
        predictor = _presentation_predictor()
    except Exception as exc:  # noqa: BLE001
        import os

        if os.environ.get("IMMUNO_ALLOW_HEURISTIC_MHC") == "1":
            out = []
            for peptide, allele in pairs:
                h = sum((i + 1) * ord(c) for i, c in enumerate(peptide + allele)) % 5000
                rank = round(h / 100.0, 2)
                out.append(
                    {
                        "peptide": peptide,
                        "allele": allele,
                        "affinity_nm": 500.0 * (1.0 + rank),
                        "presentation_score": max(0.0, 1.0 - rank / 50.0),
                        "processing_score": 0.5,
                        "percentile_rank": rank,
                    }
                )
            return out
        raise RuntimeError("MHCflurry unavailable for feature extraction") from exc

    peptides = [p for p, _ in pairs]
    # Score each unique allele separately (PresentationPredictor collapses panels).
    alleles_unique = sorted({a for _, a in pairs})
    by_key: dict[tuple[str, str], dict] = {}
    for allele in alleles_unique:
        allele_peptides = [p for p, a in pairs if a == allele]
        if not allele_peptides:
            continue
        df = predictor.predict(
            peptides=allele_peptides,
            alleles=[allele],
            verbose=0,
            include_affinity_percentile=True,
        )
        for rec in df.to_dict("records"):
            by_key[(str(rec["peptide"]), allele)] = rec

    out: list[dict[str, Any]] = []
    for peptide, allele in pairs:
        row = by_key.get((peptide, allele))
        if row is None:
            out.append(
                {
                    "peptide": peptide,
                    "allele": allele,
                    "affinity_nm": None,
                    "presentation_score": None,
                    "processing_score": None,
                    "percentile_rank": None,
                }
            )
            continue
        pct = row.get("affinity_percentile")
        if pct is None:
            pct = row.get("presentation_percentile")
        out.append(
            {
                "peptide": peptide,
                "allele": allele,
                "affinity_nm": float(row.get("affinity") or 0.0),
                "presentation_score": float(row.get("presentation_score") or 0.0),
                "processing_score": float(row.get("processing_score") or 0.0),
                "percentile_rank": float(pct) if pct is not None else None,
            }
        )
    return out
