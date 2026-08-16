"""End-to-end immuno-risk pipeline with versioned artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import uuid
from datetime import datetime, timezone
from pathlib import Path

from re_agent.immuno_risk.aggregation import aggregation_report
from re_agent.immuno_risk.mhc_ii import DEFAULT_ALLELES_II, score_mhc_ii
from re_agent.immuno_risk.mhcflurry_backend import DEFAULT_ALLELES_I, mhcflurry_version, score_mhc_i
from re_agent.immuno_risk.peptides import clean_sequence
from re_agent.immuno_risk.reference_data import ROOT, ensure_fixtures
from re_agent.immuno_risk.risk import build_confidence, join_evidence, project_residue_risk, score_risk
from re_agent.immuno_risk.schemas import ImmunoRunResult
from re_agent.immuno_risk.tolerance import check_tolerance

log = logging.getLogger(__name__)


def _run_id(sequence_id: str) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}_{sequence_id[:40]}_{uuid.uuid4().hex[:8]}"


def run_immuno_risk(
    sequence: str,
    *,
    sequence_id: str = "query",
    delivery_mode: str = "intracellular_plasmid",
    alleles_i: list[str] | None = None,
    alleles_ii: list[str] | None = None,
    mhc_class: str = "both",
    write: bool = True,
    results_root: Path | None = None,
    use_netmhcpan: bool = False,
    use_netmhciipan: bool = False,
) -> ImmunoRunResult:
    ensure_fixtures()
    seq = clean_sequence(sequence)
    if len(seq) < 8:
        raise ValueError(f"sequence too short ({len(seq)}); need ≥8 AA")

    alleles_i = alleles_i or DEFAULT_ALLELES_I
    alleles_ii = alleles_ii or DEFAULT_ALLELES_II
    caveats: list[str] = []
    versions: dict[str, str] = {}
    predictors_used: list[str] = []

    mhc_i = []
    if mhc_class in {"I", "both"}:
        try:
            mhc_i = score_mhc_i(seq, alleles_i)
            primary_i = mhc_i[0].method if mhc_i else "none"
            versions["mhc_i"] = (
                mhcflurry_version() if primary_i == "mhcflurry_presentation" else primary_i
            )
            predictors_used.append(primary_i)
        except Exception as exc:  # noqa: BLE001
            caveats.append(f"MHCflurry failed: {exc}")
            log.error("MHCflurry failed: %s", exc)
            # Fail closed on MHC-I for honesty — empty hits, high caveat
            mhc_i = []

        if use_netmhcpan:
            try:
                from re_agent.immuno_risk.netmhcpan import score_netmhcpan_i

                nm = score_netmhcpan_i(seq, alleles_i)
                versions["netmhcpan"] = "4.2e"
                predictors_used.append("netmhcpan")
                # Keep as additional evidence rows (do not replace MHCflurry)
                mhc_i = mhc_i + nm
            except Exception as exc:  # noqa: BLE001
                caveats.append(f"NetMHCpan skipped: {exc}")

    mhc_ii = []
    if mhc_class in {"II", "both"}:
        try:
            mhc_ii = score_mhc_ii(seq, alleles_ii)
            versions["mhc_ii"] = mhc_ii[0].method if mhc_ii else "none"
            predictors_used.append(versions["mhc_ii"])
        except Exception as exc:  # noqa: BLE001
            caveats.append(f"MHC-II predictor failed: {exc}")
            log.error("MHC-II predictor failed: %s", exc)
            mhc_ii = []
        if use_netmhciipan:
            try:
                from re_agent.immuno_risk.mhc_ii import score_netmhciipan

                mhc_ii = mhc_ii + score_netmhciipan(seq, alleles_ii)
                versions["netmhciipan"] = "4.3k"
            except Exception as exc:  # noqa: BLE001
                caveats.append(f"NetMHCIIpan skipped: {exc}")

    # Tolerance on unique peptides from top hits
    top_peptides = list(
        dict.fromkeys([h.peptide for h in mhc_i[:40]] + [h.peptide for h in mhc_ii[:20]])
    )
    tol_list = check_tolerance(top_peptides)
    tol_map = {t.peptide: t for t in tol_list}

    evidence = join_evidence(
        (
            [h for h in mhc_i if h.method.startswith("mhcflurry")]
            or [h for h in mhc_i if h.method == "netmhcpan"]
            or mhc_i
        ),
        tol_map,
    )
    evidence.extend(join_evidence(mhc_ii[:30], tol_map))

    agg = aggregation_report(sequence_id, seq)
    risk = score_risk(
        sequence_id,
        evidence=evidence,
        mhc_i=mhc_i,
        mhc_ii=mhc_ii,
        aggregation_score=agg.score0to100,
    )
    atlas_hits = sum(1 for t in tol_list if t.atlas_hit)
    confidence = build_confidence(
        mhc_i_hits=mhc_i,
        mhc_ii_hits=mhc_ii,
        atlas_coverage=(atlas_hits / max(len(tol_list), 1)),
        predictors_used=predictors_used,
        allele_count_i=len(alleles_i),
        allele_count_ii=len(alleles_ii),
    )
    residue_risk = project_residue_risk(seq, [e for e in evidence if e.mhc_class == "I"])

    run_id = _run_id(sequence_id)
    result = ImmunoRunResult(
        run_id=run_id,
        sequence_id=sequence_id,
        sequence=seq,
        delivery_mode=delivery_mode,
        alleles_i=alleles_i,
        alleles_ii=alleles_ii,
        peptides=evidence,
        risk=risk,
        confidence=confidence,
        residue_risk=residue_risk,
        aggregation=agg,
        predictor_versions=versions,
        caveats=caveats,
    )

    if write:
        root = results_root or (ROOT / "results" / "immuno_risk")
        out_dir = root / run_id
        out_dir.mkdir(parents=True, exist_ok=True)
        _write_artifacts(out_dir, result)
        result.artifact_dir = str(out_dir)

    return result


def _write_artifacts(out_dir: Path, result: ImmunoRunResult) -> None:
    summary = result.model_dump()
    # Keep residue list but peptides can be large — still write full JSON
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    (out_dir / "residue-risk.json").write_text(
        json.dumps([r.model_dump() for r in result.residue_risk], indent=2)
    )
    (out_dir / "aggregation.json").write_text(json.dumps(result.aggregation.model_dump(), indent=2))

    pep_path = out_dir / "peptides.csv"
    with pep_path.open("w", newline="") as f:
        fields = [
            "peptide",
            "allele",
            "mhc_class",
            "start",
            "end",
            "percentile_rank",
            "presentation_score",
            "affinity_nm",
            "presentation_points",
            "tolerance_points",
            "point_score",
            "contribution",
            "tolerance_status",
            "method",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for e in result.peptides:
            w.writerow(
                {
                    "peptide": e.peptide,
                    "allele": e.allele,
                    "mhc_class": e.mhc_class,
                    "start": e.start,
                    "end": e.end,
                    "percentile_rank": e.mhc.percentile_rank,
                    "presentation_score": e.mhc.presentation_score,
                    "affinity_nm": e.mhc.affinity_nm,
                    "presentation_points": e.presentation_points,
                    "tolerance_points": e.tolerance_points,
                    "point_score": e.point_score,
                    "contribution": e.contribution,
                    "tolerance_status": e.tolerance.status if e.tolerance else "",
                    "method": e.mhc.method,
                }
            )

    checksums = {}
    for name in ["summary.json", "peptides.csv", "residue-risk.json", "aggregation.json"]:
        p = out_dir / name
        checksums[name] = hashlib.sha256(p.read_bytes()).hexdigest()

    manifest = {
        "run_id": result.run_id,
        "sequence_id": result.sequence_id,
        "delivery_mode": result.delivery_mode,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "predictor_versions": result.predictor_versions,
        "overall_risk": result.risk.overall,
        "score0to100": result.risk.score0to100,
        "total_points": result.risk.total_points,
        "max_points": result.risk.max_points,
        "confidence": result.confidence.score0to1,
        "aggregation_overall": result.aggregation.overall,
        "checksums_sha256": checksums,
        "caveats": result.caveats,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    md = [
        f"# Immuno-risk: {result.sequence_id}",
        "",
        f"- Run ID: `{result.run_id}`",
        f"- Delivery: **{result.delivery_mode}**",
        f"- Overall epitope risk: **{result.risk.overall}** "
        f"({result.risk.total_points}/{result.risk.max_points} points; "
        f"{result.risk.score0to100}/100 normalized)",
        f"- Confidence: {result.confidence.score0to1}",
        f"- Aggregation (separate): **{result.aggregation.overall}** ({result.aggregation.score0to100}/100)",
        f"- Predictors: {result.predictor_versions}",
        "",
        "## Factors",
        *[f"- {f['name']} (+{f['contribution']}): {f['note']}" for f in result.risk.factors],
        "",
        "## Caveats",
        *[f"- {c}" for c in (result.caveats or ["none"])],
        "",
        f"Artifacts in `{out_dir}`",
        "",
    ]
    (out_dir / "report.md").write_text("\n".join(md))
