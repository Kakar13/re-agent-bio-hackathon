"""Pipeline entry point.

    uv run python -m re_agent.pda_protease.cli run --help

Stages run in pathway order: resolve cleavage, digest, present, then co-fold as
an independent check. Each stage caches, so an interrupted run resumes cheaply.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path

from . import accessibility as acc
from . import cleavage, controls, digest, geometry, merops, pda, report, structure
from .iedb import IEDBClient
from .paths import run_dir
from .pilot import select_pilot

log = logging.getLogger("pda_protease")

# netmhcpan_el / netmhciipan_el percentile below which a peptide is called a
# binder. 2.0 is a strong-binder cutoff for class II; 1.0 is conventional for
# class I but 2.0 keeps the two arms comparable.
BINDER_PERCENTILE = 2.0


def setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )
    for noisy in ("httpx", "httpcore", "urllib3"):
        logging.getLogger(noisy).setLevel(logging.WARNING)


def _binder_rows(rows: list[dict], percentile_key: str) -> list[dict]:
    out = []
    for r in rows:
        v = r.get(percentile_key)
        if v is None:
            continue
        try:
            if float(v) <= BINDER_PERCENTILE:
                out.append(r)
        except (TypeError, ValueError):
            continue
    return out


def cmd_run(args: argparse.Namespace) -> int:
    run_id = args.run_id or report.new_run_id()
    out_dir = run_dir(run_id)
    log.info("run %s -> %s", run_id, out_dir)

    manifest: dict = {
        "run_id": run_id,
        "started_at": report.new_run_id(),
        "args": vars(args),
    }

    # ---------------------------------------------------------------- corpus
    records = pda.load_pda_records()
    pool = pda.candidate_pool(records)
    background = merops.background_from_sequences([d.sequence for d in pool])
    matrices = merops.load_matrices(background)
    manifest["pool_size"] = len(pool)
    manifest["amino_acid_background"] = {k: round(v, 5) for k, v in background.items()}
    manifest["proteases"] = {
        n: {"merops": m.merops_id, "n_cleavages": m.n_cleavages, "family": m.family}
        for n, m in matrices.items()
    }

    # ----------------------------------------------------------------- pilot
    anchors = pda.natural_anchor_designs()
    entries = select_pilot(
        pool,
        matrices,
        anchors,
        n_high=args.n_high,
        n_clean=args.n_clean,
        n_mid=args.n_mid,
        n_shuffles=args.screen_shuffles,
    )
    designs = [e.design for e in entries]
    sequences = {d.design_id: d.sequence for d in designs}
    manifest["pilot"] = [
        {
            "design_id": e.design.design_id,
            "pole": e.pole,
            "length": e.design.length,
            "pdb_id": e.design.pdb_id,
            "site_density": round(e.site_density, 3),
            "n_significant": e.n_significant,
            "subtitle": e.design.subtitle,
        }
        for e in entries
    ]

    # ------------------------------------------------------------ structures
    struct_paths = pda.fetch_structures([d for d in designs if d.pdb_id])
    features: dict[str, list] = {}
    for d in designs:
        p = struct_paths.get(d.pdb_id)
        if p is None:
            continue
        features[d.design_id] = acc.compute_structure_features(p, d.chain_id, d.sequence)
    log.info("structural features for %d/%d designs", len(features), len(designs))

    # --------------------------------------------------------------- netchop
    netchop: dict[str, list[float]] = {}
    if not args.skip_iedb:
        client = IEDBClient()
        try:
            netchop = client.netchop_scores(sequences)
        except Exception as exc:  # noqa: BLE001
            log.warning("netchop failed, continuing without class I arm: %s", exc)

    # -------------------------------------------------------------- cleavage
    all_sites: dict[str, list[cleavage.CleavageSite]] = {}
    for d in designs:
        all_sites[d.design_id] = cleavage.scan_design(
            d.design_id,
            d.sequence,
            matrices,
            features=features.get(d.design_id),
            netchop=netchop.get(d.design_id),
            n_shuffles=args.shuffles,
        )
    flat_sites = [s for v in all_sites.values() for s in v]
    significant = [s for s in flat_sites if s.p_value <= args.p_threshold]
    log.info("%d bonds scored, %d significant at p<=%.3g",
             len(flat_sites), len(significant), args.p_threshold)

    # ---------------------------------------------------------------- digest
    peptides: list[digest.Peptide] = []
    unconstrained: list[digest.Peptide] = []
    for d in designs:
        sites = all_sites[d.design_id]
        cuts = cleavage.high_confidence_cuts(sites, p_threshold=args.p_threshold)
        peptides += digest.digest_mhcii(d.design_id, d.sequence, cuts, sites)
        unconstrained += digest.unconstrained_peptides(d.design_id, d.sequence)
        nc = netchop.get(d.design_id)
        if nc:
            peptides += digest.digest_mhci(
                d.design_id, d.sequence, cleavage.netchop_cut_positions(nc)
            )
    peptides = digest.dedupe(peptides)
    log.info("digest produced %d peptides (%d unconstrained baseline)",
             len(peptides), len(unconstrained))

    # ------------------------------------------------------------ presentation
    peptide_summary: dict = {}
    scored_by_peptide: dict[str, float] = {}
    if not args.skip_iedb:
        client = IEDBClient()
        ii = [p for p in peptides if p.arm == digest.ARM_MHCII]
        i1 = [p for p in peptides if p.arm == digest.ARM_MHCI]
        try:
            if ii:
                rows = client.score_peptides(
                    [p.peptide for p in ii][: args.max_peptides], tool_group="mhcii"
                )
                binders = _binder_rows(rows, "netmhciipan_el_percentile")
                peptide_summary["mhcii_rows"] = len(rows)
                peptide_summary["constrained_binders"] = len({b["peptide"] for b in binders})
                for r in rows:
                    v = r.get("netmhciipan_el_percentile")
                    if v is not None:
                        k = r["peptide"]
                        scored_by_peptide[k] = min(scored_by_peptide.get(k, 1e9), float(v))
            if i1:
                rows = client.score_peptides(
                    [p.peptide for p in i1][: args.max_peptides], tool_group="mhci"
                )
                peptide_summary["mhci_rows"] = len(rows)
                peptide_summary["mhci_binders"] = len(
                    {b["peptide"] for b in _binder_rows(rows, "netmhcpan_el_percentile")}
                )
            if args.baseline:
                base = sorted({p.peptide for p in unconstrained})[: args.max_peptides]
                rows = client.score_peptides(base, tool_group="mhcii")
                peptide_summary["unconstrained_binders"] = len(
                    {b["peptide"] for b in _binder_rows(rows, "netmhciipan_el_percentile")}
                )
        except Exception as exc:  # noqa: BLE001
            log.warning("IEDB presentation scoring failed: %s", exc)

    # ------------------------------------------------------------- structure
    geometries: list[geometry.Geometry] = []
    positive_control: list[dict] = []
    if not args.skip_structure:
        real_sites = []
        for d in designs:
            if d.design_id == "anchor_CD74":
                continue
            best = [
                s
                for s in all_sites[d.design_id]
                if s.p_value <= args.p_threshold and s.protease in structure.COFOLDABLE
            ]
            best.sort(key=lambda s: s.p_value)
            for s in best[: args.sites_per_design]:
                real_sites.append((d.design_id, s.protease, s))
        real_sites = real_sites[: args.max_cofolds]

        plan = controls.build_control_plan(real_sites, all_sites, sequences)
        jobs = plan.jobs

        cd74 = next((d for d in designs if d.design_id == "anchor_CD74"), None)
        if cd74 is not None:
            jobs = controls.positive_control_jobs(
                cd74.sequence, all_sites[cd74.design_id], n=args.positive_controls
            ) + jobs

        if args.max_jobs:
            jobs = jobs[: args.max_jobs]
        manifest["cofold_jobs"] = len(jobs)
        manifest["cofold_arms"] = plan.counts

        results = structure.run_wave(
            jobs,
            use_msa=not args.no_msa,
            recycling_steps=args.recycling_steps,
            sampling_steps=args.sampling_steps,
        )
        geometries = geometry.analyse_all(results)
        by_job = {j.job_id: j for j in jobs}
        for g in geometries:
            j = by_job.get(g.job_id)
            if j and j.arm == controls.ARM_POSITIVE:
                positive_control.append(
                    g.to_dict() | {"in_clip_window": j.meta.get("in_clip_window", False)}
                )

    # ------------------------------------------------------------ clean gate
    clean = _evaluate_clean_designs(entries, all_sites, args.p_threshold)

    # ---------------------------------------------------------------- report
    gates = report.evaluate_gates(
        geometries=[g.to_dict() for g in geometries],
        peptide_summary=peptide_summary,
        positive_control=positive_control,
        clean_designs=clean,
    )
    summary = {
        "designs_in_pilot": len(designs),
        "bonds_scored": len(flat_sites),
        "significant_sites": len(significant),
        "peptides_after_digest": len(peptides),
        "unconstrained_baseline_peptides": len(unconstrained),
        "cofolds_completed": len(geometries),
        **peptide_summary,
    }
    report.write_report(
        out_dir,
        manifest=manifest,
        sites=[s.to_dict() for s in sorted(significant, key=lambda s: s.p_value)],
        peptides=[p.to_dict() for p in peptides],
        geometries=[g.to_dict() for g in geometries],
        gates=gates,
        summary=summary,
    )

    print()
    print(f"run {run_id} complete -> {out_dir}")
    for g in gates:
        print(f"  [{g.verdict:>13}] {g.gate_id}: {g.detail}")
    return 0


def _evaluate_clean_designs(entries, all_sites, p_threshold: float) -> dict:
    """Do the designs picked as predicted-clean actually stay clean?"""
    clean = [e for e in entries if e.pole == "predicted-clean"]
    high = [e for e in entries if e.pole == "high-signal"]
    if not clean or not high:
        return {"evaluated": False, "detail": "pilot lacks both poles"}

    def density(e) -> float:
        sites = all_sites.get(e.design.design_id, [])
        bonds = {s.cut_index for s in sites if s.p_value <= p_threshold}
        return 100.0 * len(bonds) / max(e.design.length, 1)

    dc = [density(e) for e in clean]
    dh = [density(e) for e in high]
    mc, mh = sum(dc) / len(dc), sum(dh) / len(dh)
    return {
        "evaluated": True,
        "passed": mc < mh,
        "detail": (
            f"predicted-clean designs carry {mc:.2f} significant sites per 100 aa "
            f"against {mh:.2f} for high-signal designs, on the full structure-aware scan"
        ),
        "value": {"clean_density": mc, "high_density": mh},
    }


def cmd_smoke(args: argparse.Namespace) -> int:
    """One co-fold, to prove the Modal path works before committing to a wave."""
    seq = pda.fetch_uniprot_sequence("P04233")
    job = structure.build_job(
        design_id="CD74",
        sequence=seq,
        cut_index=120,
        protease="cathepsin_S",
        arm=controls.ARM_POSITIVE,
    )
    rec = structure.run_cofold(
        job,
        use_msa=not args.no_msa,
        recycling_steps=args.recycling_steps,
        sampling_steps=args.sampling_steps,
    )
    if not rec:
        print("SMOKE-FAIL: no structure returned")
        return 1
    g = geometry.analyse(rec)
    print(json.dumps({"metrics": rec["metrics"], "geometry": asdict(g) if g else None}, indent=2))
    print("SMOKE-OK")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="pda_protease", description=__doc__)
    ap.add_argument("-v", "--verbose", action="store_true")
    sub = ap.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run the pipeline")
    r.add_argument("--run-id")
    r.add_argument("--n-high", type=int, default=6)
    r.add_argument("--n-clean", type=int, default=4)
    r.add_argument("--n-mid", type=int, default=3)
    r.add_argument("--shuffles", type=int, default=200, help="null shuffles in the full scan")
    r.add_argument("--screen-shuffles", type=int, default=50, help="null shuffles when ranking")
    r.add_argument("--p-threshold", type=float, default=0.01)
    r.add_argument("--max-peptides", type=int, default=800)
    r.add_argument("--sites-per-design", type=int, default=1)
    r.add_argument("--max-cofolds", type=int, default=8, help="real sites to co-fold")
    r.add_argument("--max-jobs", type=int, default=0, help="hard cap on all co-folds, 0 = no cap")
    r.add_argument("--positive-controls", type=int, default=2)
    r.add_argument("--recycling-steps", type=int, default=None)
    r.add_argument("--sampling-steps", type=int, default=None)
    r.add_argument("--no-msa", action="store_true")
    r.add_argument("--skip-structure", action="store_true")
    r.add_argument("--skip-iedb", action="store_true")
    r.add_argument("--baseline", action="store_true", help="also score the unconstrained set")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("smoke", help="single co-fold sanity check")
    s.add_argument("--no-msa", action="store_true")
    s.add_argument("--recycling-steps", type=int, default=None)
    s.add_argument("--sampling-steps", type=int, default=None)
    s.set_defaults(func=cmd_smoke)

    args = ap.parse_args(argv)
    setup_logging(args.verbose)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
