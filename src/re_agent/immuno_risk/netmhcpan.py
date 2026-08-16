"""Optional NetMHCpan-4.2e local adapter (academic license — user-provided binary)."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from re_agent.immuno_risk.peptides import clean_sequence
from re_agent.immuno_risk.schemas import PeptideHit

PINNED_VERSION = "4.2e"


def netmhcpan_available() -> bool:
    path = os.environ.get("NETMHCPAN_BIN") or shutil.which("netMHCpan")
    return bool(path and Path(path).exists())


def netmhcpan_bin() -> Path:
    env = os.environ.get("NETMHCPAN_BIN")
    if env:
        return Path(env)
    which = shutil.which("netMHCpan")
    if which:
        return Path(which)
    raise FileNotFoundError(
        "NetMHCpan not found. Set NETMHCPAN_BIN to a licensed 4.2e binary "
        "(academic download from DTU). Do not automate the public web form."
    )


def _parse_netmhcpan_table(text: str) -> list[PeptideHit]:
    """Parse NetMHCpan stdout, retaining the Identity (FASTA id) column.

    Typical columns (whitespace-separated)::

        Pos HLA Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL [BindLevel]

    Identity is the last non-numeric token before the trailing Score_EL / %Rank_EL floats.
    """
    hits: list[PeptideHit] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("-"):
            continue
        if "Peptide" in line and ("%Rank" in line or "Rnk_EL" in line):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        try:
            pos = int(parts[0])
        except ValueError:
            continue
        allele = parts[1]
        peptide = parts[2]
        if not peptide.isalpha():
            continue

        # Walk tokens after peptide; collect trailing floats (Score_EL, %Rank_EL, …)
        # and remember the last alphabetic/underscore token before those floats = Identity.
        binder = False
        # Strip bind-level tags from the right
        toks = list(parts[3:])
        while toks and (toks[-1] in {"SB", "WB", "<="} or toks[-1].startswith("<=")):
            binder = True
            if toks[-1] in {"SB", "WB"}:
                toks.pop()
            elif toks[-1].startswith("<="):
                toks.pop()
            if toks and toks[-1] == "<=":
                toks.pop()

        # Trailing floats
        float_idxs: list[int] = []
        for i in range(len(toks) - 1, -1, -1):
            try:
                float(toks[i])
                float_idxs.append(i)
            except ValueError:
                break
        float_idxs.reverse()
        if len(float_idxs) < 2:
            continue
        el_score = float(toks[float_idxs[-2]])
        rank = float(toks[float_idxs[-1]])
        # Identity = last non-float token before the float block
        first_float = float_idxs[0]
        identity: str | None = None
        for i in range(first_float - 1, -1, -1):
            try:
                float(toks[i])
                continue
            except ValueError:
                identity = toks[i]
                break

        binder = binder or rank <= 2.0
        provenance: dict = {
            "predictor": "netMHCpan",
            "pinned": PINNED_VERSION,
        }
        if identity is not None:
            provenance["identity"] = identity
        hits.append(
            PeptideHit(
                peptide=peptide,
                allele=allele,
                mhc_class="I",
                start=pos,
                end=pos + len(peptide),
                length=len(peptide),
                el_score=el_score,
                percentile_rank=rank,
                binder=binder,
                method="netmhcpan",
                version=PINNED_VERSION,
                provenance=provenance,
                caveat=(
                    "Optional licensed comparator. EL/BA ≠ immunogenicity. "
                    "Pathogen/neo heads are domain-mismatched for de novo."
                ),
            )
        )
    return hits


def score_netmhcpan_i(
    sequence: str,
    alleles: list[str],
    *,
    lengths: list[int] | None = None,
    include_pathogen: bool = False,
    include_neo: bool = False,
    top_n: int = 50,
) -> list[PeptideHit]:
    """Run local NetMHCpan-4.2e. Raises if binary missing (fail closed)."""
    bin_path = netmhcpan_bin()
    seq = clean_sequence(sequence)
    lengths = lengths or list(range(8, 12))
    allele_str = ",".join(a.replace("*", "") for a in alleles)
    length_str = ",".join(str(L) for L in lengths)

    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "query.fsa"
        fasta.write_text(f">query\n{seq}\n")
        cmd = [
            str(bin_path),
            "-f",
            str(fasta),
            "-a",
            allele_str,
            "-l",
            length_str,
        ]
        if include_pathogen:
            cmd.append("-pathogen")
        if include_neo:
            cmd.append("-neo")
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise RuntimeError(f"NetMHCpan execution failed: {exc}") from exc
        if proc.returncode != 0:
            raise RuntimeError(
                f"NetMHCpan exited {proc.returncode}: {(proc.stderr or proc.stdout)[:500]}"
            )
        hits = _parse_netmhcpan_table(proc.stdout)
        # Prefer %Rank sort
        hits.sort(key=lambda h: h.percentile_rank if h.percentile_rank is not None else 99.0)
        return hits[:top_n]
