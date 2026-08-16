"""Lightweight aggregation / self-association report (separate from epitope risk)."""

from __future__ import annotations

from re_agent.immuno_risk.peptides import HYDROPHOBIC, clean_sequence, net_charge_ph74
from re_agent.immuno_risk.schemas import AggregationReport


def aggregation_report(sequence_id: str, sequence: str) -> AggregationReport:
    seq = clean_sequence(sequence)
    n = max(len(seq), 1)
    hydrophobic_fraction = sum(1 for c in seq if c in HYDROPHOBIC) / n
    charge = net_charge_ph74(seq)
    free_cys = seq.count("C")  # reducing cytosol → disulfides unreliable
    # β-edge proxy: runs of β-branched / aromatic suggestive of edge stacking risk
    betaish = set("VIFWY")
    runs = 0
    i = 0
    while i < len(seq):
        if seq[i] in betaish:
            j = i
            while j < len(seq) and seq[j] in betaish:
                j += 1
            if j - i >= 4:
                runs += 1
            i = j
        else:
            i += 1
    beta_edge_proxy = min(1.0, runs / 3.0)
    # Solubility proxy: prefer charge away from zero + lower hydrophobics
    solubility_proxy = min(
        1.0,
        (abs(charge) / 10.0) * 0.5 + (1.0 - hydrophobic_fraction) * 0.5,
    )

    factors: list[dict] = []
    score = 0.0
    if hydrophobic_fraction > 0.45:
        c = min(35.0, (hydrophobic_fraction - 0.45) * 200)
        factors.append(
            {
                "name": "hydrophobic_surface_proxy",
                "contribution": round(c, 1),
                "note": f"Hydrophobic fraction {hydrophobic_fraction:.2f}",
            }
        )
        score += c
    if abs(charge) < 2:
        factors.append(
            {
                "name": "near_neutral_charge",
                "contribution": 15.0,
                "note": f"Net charge at pH~7.4 ≈ {charge:.0f} (colloidal risk)",
            }
        )
        score += 15.0
    if free_cys > 0:
        c = min(20.0, free_cys * 5.0)
        factors.append(
            {
                "name": "free_cysteines_reducing_cytosol",
                "contribution": c,
                "note": f"{free_cys} Cys — disulfide-dependent stability unreliable intracellularly",
            }
        )
        score += c
    if beta_edge_proxy > 0:
        c = beta_edge_proxy * 20.0
        factors.append(
            {
                "name": "beta_edge_proxy",
                "contribution": round(c, 1),
                "note": f"{runs} hydrophobic/β-branched runs ≥4",
            }
        )
        score += c
    if solubility_proxy < 0.4:
        factors.append(
            {
                "name": "low_solubility_proxy",
                "contribution": 10.0,
                "note": f"Solubility proxy {solubility_proxy:.2f}",
            }
        )
        score += 10.0

    score = min(100.0, score)
    overall = "high" if score >= 55 else "moderate" if score >= 30 else "low"
    return AggregationReport(
        sequence_id=sequence_id,
        overall=overall,  # type: ignore[arg-type]
        score0to100=round(score, 1),
        factors=factors,
        net_charge_ph74=charge,
        hydrophobic_fraction=round(hydrophobic_fraction, 3),
        free_cysteine_count=free_cys,
        beta_edge_proxy=round(beta_edge_proxy, 3),
        solubility_proxy=round(solubility_proxy, 3),
        persistence_caveat=(
            "Protease resistance may increase exposure time but does not imply aggregation; "
            "report persistence separately from self-association."
        ),
    )
