"""Thin MHC-II presentation: HLAIIPred (required) ± optional NetMHCIIpan-4.3k."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from re_agent.immuno_risk.peptides import clean_sequence, sliding_windows
from re_agent.immuno_risk.schemas import PeptideHit

log = logging.getLogger(__name__)

DEFAULT_ALLELES_II = [
    "DRB1*01:01",
    "DRB1*03:01",
    "DRB1*04:01",
    "DRB1*07:01",
    "DRB1*11:01",
    "DRB1*15:01",
]

PINNED_NETMHCII = "4.3k"


def hlaipred_available() -> bool:
    try:
        import hlaipred  # noqa: F401

        return True
    except ImportError:
        try:
            import HLAIIPred  # noqa: F401

            return True
        except ImportError:
            return False


def _heuristic_allowed() -> bool:
    return os.environ.get("IMMUNO_ALLOW_HEURISTIC_MHC", "").lower() in {"1", "true", "yes"}


def score_mhc_ii(
    sequence: str,
    alleles: list[str] | None = None,
    *,
    peptide_lengths: list[int] | None = None,
    top_n: int = 40,
) -> list[PeptideHit]:
    """Thin MHC-II presentation baseline."""
    alleles = alleles or DEFAULT_ALLELES_II
    lengths = peptide_lengths or list(range(13, 22))
    seq = clean_sequence(sequence)
    if len(seq) < 13:
        return []

    if hlaipred_available():
        return _score_hlaipred(seq, alleles, lengths, top_n)
    if netmhciipan_available():
        return score_netmhciipan(seq, alleles, top_n=top_n)

    if _heuristic_allowed():
        log.warning("No MHC-II predictor installed; using explicit demo-only heuristic ranks")
        return _heuristic_ii(seq, alleles, lengths, top_n)
    raise RuntimeError(
        "No MHC-II predictor available. Install HLAIIPred or set NETMHCIIPAN_BIN. "
        "For offline demos only, set IMMUNO_ALLOW_HEURISTIC_MHC=1."
    )


def _score_hlaipred(
    seq: str,
    alleles: list[str],
    lengths: list[int],
    top_n: int,
) -> list[PeptideHit]:
    """Call HLAIIPred if API shape allows; otherwise fall back gracefully."""
    try:
        # HLAIIPred packaging varies; try common entry points
        try:
            from hlaipred import predict as hla_predict  # type: ignore
        except ImportError:
            from HLAIIPred import predict as hla_predict  # type: ignore

        windows = sliding_windows(seq, lengths)
        hits: list[PeptideHit] = []
        for start, end, peptide in windows:
            for allele in alleles:
                try:
                    out = hla_predict(peptide, allele)
                except TypeError:
                    out = hla_predict(peptides=[peptide], alleles=[allele])
                score = float(out[0]) if isinstance(out, (list, tuple)) else float(out)
                # Interpret as presentation probability-ish; convert to rank proxy
                rank = max(0.01, (1.0 - score) * 100.0)
                hits.append(
                    PeptideHit(
                        peptide=peptide,
                        allele=allele,
                        mhc_class="II",
                        start=start,
                        end=end,
                        length=len(peptide),
                        presentation_score=score,
                        percentile_rank=rank,
                        binder=rank <= 10.0,
                        method="hlaipred",
                        version="installed",
                        provenance={"predictor": "HLAIIPred"},
                        caveat="MHC-II presentation only — not ADA / T-cell risk.",
                    )
                )
        hits.sort(key=lambda h: h.percentile_rank if h.percentile_rank is not None else 99.0)
        return hits[:top_n]
    except Exception as exc:  # noqa: BLE001
        if _heuristic_allowed():
            log.warning("HLAIIPred call failed (%s); using explicit demo-only heuristic", exc)
            return _heuristic_ii(seq, alleles, lengths, top_n)
        raise RuntimeError(f"HLAIIPred call failed: {exc}") from exc


def _heuristic_ii(
    seq: str,
    alleles: list[str],
    lengths: list[int],
    top_n: int,
) -> list[PeptideHit]:
    hits: list[PeptideHit] = []
    for start, end, peptide in sliding_windows(seq, lengths):
        for allele in alleles:
            # Deterministic pseudo-rank from composition (inspectable, not biology)
            h = sum((i + 1) * ord(c) for i, c in enumerate(peptide + allele)) % 5000
            rank = round(h / 100.0, 2)  # 0–49.99
            hits.append(
                PeptideHit(
                    peptide=peptide,
                    allele=allele,
                    mhc_class="II",
                    start=start,
                    end=end,
                    length=len(peptide),
                    percentile_rank=rank,
                    binder=rank <= 10.0,
                    method="mhc_ii_heuristic_v0",
                    version="fallback",
                    provenance={"predictor": "heuristic_fallback"},
                    caveat=(
                        "HLAIIPred not available — heuristic rank only. "
                        "Install HLAIIPred or provide NetMHCIIpan for real presentation scores."
                    ),
                )
            )
    hits.sort(key=lambda h: h.percentile_rank if h.percentile_rank is not None else 99.0)
    return hits[:top_n]


def netmhciipan_available() -> bool:
    path = os.environ.get("NETMHCIIPAN_BIN") or shutil.which("netMHCIIpan")
    return bool(path and Path(path).exists())


def score_netmhciipan(
    sequence: str,
    alleles: list[str],
    *,
    top_n: int = 40,
) -> list[PeptideHit]:
    env = os.environ.get("NETMHCIIPAN_BIN") or shutil.which("netMHCIIpan")
    if not env:
        raise FileNotFoundError(
            "NetMHCIIpan not found. Set NETMHCIIPAN_BIN to licensed 4.3k binary."
        )
    seq = clean_sequence(sequence)
    allele_str = ",".join(a.replace("*", "") for a in alleles)
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "query.fsa"
        fasta.write_text(f">query\n{seq}\n")
        proc = subprocess.run(
            [env, "-f", str(fasta), "-a", allele_str],
            capture_output=True,
            text=True,
            timeout=600,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"NetMHCIIpan failed: {(proc.stderr or proc.stdout)[:400]}")
        hits: list[PeptideHit] = []
        for line in proc.stdout.splitlines():
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                pos = int(parts[0])
            except ValueError:
                continue
            allele, peptide = parts[1], parts[2]
            floats = []
            for tok in parts[3:]:
                try:
                    floats.append(float(tok))
                except ValueError:
                    pass
            if len(floats) < 2:
                continue
            rank = floats[-1]
            hits.append(
                PeptideHit(
                    peptide=peptide,
                    allele=allele,
                    mhc_class="II",
                    start=pos,
                    end=pos + len(peptide),
                    length=len(peptide),
                    el_score=floats[-2],
                    percentile_rank=rank,
                    binder=rank <= 10.0,
                    method="netmhciipan",
                    version=PINNED_NETMHCII,
                    provenance={"predictor": "netMHCIIpan", "pinned": PINNED_NETMHCII},
                    caveat="Optional licensed MHC-II comparator — presentation only.",
                )
            )
        hits.sort(key=lambda h: h.percentile_rank if h.percentile_rank is not None else 99.0)
        return hits[:top_n]
