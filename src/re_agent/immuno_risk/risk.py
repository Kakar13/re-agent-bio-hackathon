"""Risk aggregation, confidence, and residue projection."""

from __future__ import annotations

from re_agent.immuno_risk.schemas import (
    ArmSummary,
    ConfidenceReport,
    PeptideEvidence,
    PeptideHit,
    ResidueRisk,
    RiskBreakdown,
    ToleranceEvidence,
)


def project_residue_risk(
    sequence: str,
    evidence: list[PeptideEvidence],
) -> list[ResidueRisk]:
    scores = [0.0] * len(sequence)
    counts = [0] * len(sequence)
    peptides_at: list[set[str]] = [set() for _ in sequence]
    for ev in evidence:
        if ev.start is None or ev.end is None:
            # try locate first occurrence
            idx = sequence.find(ev.peptide)
            if idx < 0:
                continue
            start, end = idx, idx + len(ev.peptide)
        else:
            start, end = ev.start, ev.end
        weight = ev.contribution or (ev.point_score / 5.0)
        for i in range(start, min(end, len(sequence))):
            scores[i] += weight
            counts[i] += 1
            peptides_at[i].add(ev.peptide)
    out: list[ResidueRisk] = []
    for i, aa in enumerate(sequence):
        out.append(
            ResidueRisk(
                position=i + 1,
                residue=aa,
                risk=round(scores[i], 4),
                peptide_count=counts[i],
                peptides=sorted(peptides_at[i])[:20],
            )
        )
    return out


def build_confidence(
    *,
    mhc_i_hits: list[PeptideHit],
    mhc_ii_hits: list[PeptideHit],
    atlas_coverage: float,
    predictors_used: list[str],
    allele_count_i: int,
    allele_count_ii: int,
) -> ConfidenceReport:
    factors = []
    score = 0.2
    real_i = any("heuristic" not in h.method for h in mhc_i_hits)
    real_ii = any("heuristic" not in h.method for h in mhc_ii_hits)
    if real_i:
        score += 0.3
        factors.append({"name": "existing_mhc_i_predictor", "contribution": 0.3})
    else:
        factors.append({"name": "mhc_i_predictor_missing_or_demo", "contribution": 0.0})
    if real_ii:
        score += 0.2
        factors.append({"name": "existing_mhc_ii_predictor", "contribution": 0.2})
    else:
        factors.append({"name": "mhc_ii_predictor_missing_or_demo", "contribution": 0.0})
    score += min(0.1, atlas_coverage * 0.1)
    factors.append({"name": "atlas_coverage", "contribution": round(min(0.1, atlas_coverage * 0.1), 3)})
    allele_coverage = min(0.1, (allele_count_i + allele_count_ii) / 120.0)
    score += allele_coverage
    factors.append({"name": "allele_panel_coverage", "contribution": round(allele_coverage, 3)})
    if len(predictors_used) >= 2:
        score += 0.1
        factors.append({"name": "both_mhc_arms", "contribution": 0.1})
    score = min(max(score, 0.0), 1.0)
    return ConfidenceReport(score0to1=round(score, 3), factors=factors)


def score_risk(
    sequence_id: str,
    *,
    evidence: list[PeptideEvidence],
    mhc_i: list[PeptideHit],
    mhc_ii: list[PeptideHit],
    aggregation_score: float,
) -> RiskBreakdown:
    factors: list[dict] = []
    i_evidence = sorted(
        (e for e in evidence if e.mhc_class == "I"),
        key=lambda e: e.point_score,
        reverse=True,
    )
    ii_evidence = sorted(
        (e for e in evidence if e.mhc_class == "II"),
        key=lambda e: e.point_score,
        reverse=True,
    )
    top_i = i_evidence[:5]
    top_ii = ii_evidence[:5]

    # Transparent point system: each peptide earns 0–3 presentation points
    # plus 0–2 tolerance/novelty points. MHC-I is weighted 2x because
    # intracellular/plasmid delivery is the primary route.
    i_points = sum(e.point_score for e in top_i) * 2
    ii_points = sum(e.point_score for e in top_ii)
    i_binders = [e for e in i_evidence if e.mhc.binder]
    foreign_i = [e for e in i_binders if e.tolerance and e.tolerance.status == "foreign_like"]
    factors.append(
        {
            "name": "mhc_i_top5_points_x2",
            "contribution": i_points,
            "note": f"Top five MHC-I peptide points weighted 2x; {len(i_binders)} binders",
        }
    )
    factors.append(
        {
            "name": "mhc_ii_top5_points",
            "contribution": ii_points,
            "note": "Top five MHC-II presentation/tolerance points; not an ADA prediction",
        }
    )

    allele_breadth = len(
        {
            e.allele
            for e in evidence
            if e.mhc.binder and e.point_score >= 2
        }
    )
    breadth_points = min(10, allele_breadth)
    factors.append(
        {
            "name": "hla_breadth_points",
            "contribution": breadth_points,
            "note": f"Flagged peptides span {allele_breadth} HLA alleles",
        }
    )

    total_points = min(85, i_points + ii_points + breadth_points)
    score = round((total_points / 85.0) * 100.0, 1)

    # Aggregation deliberately stays outside the immuno point total.
    if aggregation_score >= 55:
        factors.append(
            {
                "name": "aggregation_separate",
                "contribution": 0,
                "note": "High aggregation flag reported separately; it does not change immuno points",
            }
        )
    overall = "high" if score >= 55 else "moderate" if score >= 30 else "low"

    def arm(mhc_class: str, hits: list[PeptideHit], evid: list[PeptideEvidence]) -> ArmSummary:
        top = hits[:10]
        return ArmSummary(
            mhc_class=mhc_class,  # type: ignore[arg-type]
            binder_count=sum(1 for h in hits if h.binder),
            top_hits=top,
            total_points=sum(e.point_score for e in sorted(evid, key=lambda e: e.point_score, reverse=True)[:5]),
            max_peptide_points=max((e.point_score for e in evid), default=0),
            method="existing_predictor_points_v1",
            caveat="Transparent points over existing predictor ranks; no new model is trained.",
        )

    flagged = list(
        dict.fromkeys(
            [e.peptide for e in foreign_i]
            + [e.peptide for e in (top_i + top_ii) if e.point_score >= 2]
        )
    )[:40]

    return RiskBreakdown(
        sequence_id=sequence_id,
        overall=overall,  # type: ignore[arg-type]
        score0to100=round(score, 1),
        mhc_i=arm("I", mhc_i, [e for e in evidence if e.mhc_class == "I"]),
        mhc_ii=arm("II", mhc_ii, [e for e in evidence if e.mhc_class == "II"]),
        total_points=total_points,
        max_points=85,
        factors=factors,
        peptides_flagged=flagged,
        method="dual_arm_existing_predictor_points_v1",
        caveat=(
            "Rule-based screening points over existing MHC-I and MHC-II predictors. "
            "No model is trained; this is not a clinical immunogenicity probability."
        ),
    )


def join_evidence(
    hits: list[PeptideHit],
    tolerance_by_peptide: dict[str, ToleranceEvidence],
) -> list[PeptideEvidence]:
    out: list[PeptideEvidence] = []
    for h in hits:
        rank = h.percentile_rank
        if rank is not None and rank <= 0.5:
            presentation_points = 3
        elif rank is not None and rank <= 2.0:
            presentation_points = 2
        elif rank is not None and rank <= (10.0 if h.mhc_class == "II" else 5.0):
            presentation_points = 1
        elif h.presentation_score is not None and h.presentation_score >= 0.9:
            presentation_points = 3
        elif h.presentation_score is not None and h.presentation_score >= 0.7:
            presentation_points = 2
        elif h.presentation_score is not None and h.presentation_score >= 0.5:
            presentation_points = 1
        else:
            presentation_points = 0
        tol = tolerance_by_peptide.get(h.peptide)
        tolerance_points = (
            2 if tol and tol.status == "foreign_like"
            else 1 if tol is None or tol.status == "unknown"
            else 0
        )
        point_score = presentation_points + tolerance_points
        out.append(
            PeptideEvidence(
                peptide=h.peptide,
                allele=h.allele,
                mhc_class=h.mhc_class,
                start=h.start,
                end=h.end,
                mhc=h,
                tolerance=tol,
                presentation_points=presentation_points,
                tolerance_points=tolerance_points,
                point_score=point_score,
                contribution=round(point_score / 5.0, 4),
            )
        )
    return out
