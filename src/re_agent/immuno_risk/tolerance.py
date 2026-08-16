"""Tolerance / self-likeness vs HLA Ligand Atlas-style references."""

from __future__ import annotations

from re_agent.immuno_risk.reference_data import load_atlas_peptides
from re_agent.immuno_risk.schemas import ToleranceEvidence


def _identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / max(len(a), len(b))


def check_tolerance(
    peptides: list[str],
    *,
    allele: str | None = None,
) -> list[ToleranceEvidence]:
    atlas = load_atlas_peptides()
    hits: list[ToleranceEvidence] = []
    for peptide in peptides:
        p = peptide.upper()
        if p in atlas:
            hits.append(
                ToleranceEvidence(
                    peptide=p,
                    allele=allele,
                    status="self_like",
                    nearest_self=p,
                    identity=1.0,
                    atlas_hit=True,
                    method="hla_ligand_atlas_exact",
                    caveat="Observed benign-tissue ligand — supports self presentation, not absolute tolerance.",
                )
            )
            continue
        best, best_id = None, 0.0
        for self_p in atlas:
            if abs(len(self_p) - len(p)) > 2:
                continue
            ident = _identity(p, self_p)
            if ident > best_id:
                best_id, best = ident, self_p
        if best_id >= 0.85:
            status = "self_like"
        elif best_id <= 0.4:
            status = "foreign_like"
        else:
            status = "unknown"
        hits.append(
            ToleranceEvidence(
                peptide=p,
                allele=allele,
                status=status,  # type: ignore[arg-type]
                nearest_self=best,
                identity=round(best_id, 3),
                atlas_hit=False,
                method="atlas_identity_v1",
                caveat="Sequence identity to Atlas ligands is a proxy — not JanusMatrix TCR-face logic.",
            )
        )
    return hits
