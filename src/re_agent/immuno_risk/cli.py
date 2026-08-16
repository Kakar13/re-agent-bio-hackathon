"""CLI: uv run python -m re_agent.immuno_risk.cli ..."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from re_agent.immuno_risk.peptides import clean_sequence
from re_agent.immuno_risk.pipeline import run_immuno_risk
from re_agent.immuno_risk.reference_data import ensure_fixtures, fetch_atlas_stub, fetch_iedb_tcell_sample


def _read_sequence(raw: str | None, fasta: Path | None) -> tuple[str, str]:
    if fasta:
        text = fasta.read_text()
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        sid = "query"
        seq_parts: list[str] = []
        for ln in lines:
            if ln.startswith(">"):
                if seq_parts:
                    break
                sid = ln[1:].split()[0]
            else:
                seq_parts.append(ln)
        return sid, clean_sequence("".join(seq_parts))
    if not raw:
        raise SystemExit("provide --sequence or --fasta")
    if raw.startswith(">"):
        lines = raw.splitlines()
        sid = lines[0][1:].split()[0]
        return sid, clean_sequence("".join(lines[1:]))
    return "query", clean_sequence(raw)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Immuno-risk screening CLI")
    sub = p.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser("run", help="Run dual-arm immuno-risk pipeline")
    run_p.add_argument("--sequence", type=str, default=None)
    run_p.add_argument("--fasta", type=Path, default=None)
    run_p.add_argument("--sequence-id", type=str, default=None)
    run_p.add_argument("--delivery-mode", default="intracellular_plasmid")
    run_p.add_argument("--mhc-class", choices=["I", "II", "both"], default="both")
    run_p.add_argument("--alleles-i", type=str, default=None, help="comma-separated HLA-I")
    run_p.add_argument("--alleles-ii", type=str, default=None, help="comma-separated HLA-II")
    run_p.add_argument("--no-write", action="store_true")
    run_p.add_argument("--netmhcpan", action="store_true")
    run_p.add_argument("--netmhciipan", action="store_true")
    run_p.add_argument("--json", action="store_true", help="print full JSON to stdout")

    sub.add_parser("ensure-fixtures", help="Create offline IEDB/Atlas fixtures + manifests")
    fetch_p = sub.add_parser("fetch-refs", help="Best-effort fetch IEDB sample + Atlas stub")
    fetch_p.add_argument("--limit", type=int, default=200)

    sub.add_parser(
        "train",
        help="Train the IEDB logistic MHC-I risk head (not NetMHCpan)",
    )

    pda_p = sub.add_parser(
        "pda-ingest",
        help="Pull Protein Design Archive designed chains into NetMHCpan-4.2e FASTA",
    )
    pda_p.add_argument("--source", choices=["github", "api", "json"], default="github")
    pda_p.add_argument("--json", type=Path, default=None, dest="pda_json")
    pda_p.add_argument("--out", type=Path, default=None)
    pda_p.add_argument("--limit", type=int, default=None)
    pda_p.add_argument("--pdb-codes", type=str, default=None)
    pda_p.add_argument("--include-unknown", action="store_true")
    pda_p.add_argument("--keep-his-tags", action="store_true")
    pda_p.add_argument("--peptides", action="store_true", help="also tile 8–11mer MHC-I peptides")

    cohort_p = sub.add_parser(
        "score-cohort",
        help="Batch MHC-I score all unique peptides in a FASTA (mhcflurry|netmhcpan)",
    )
    cohort_p.add_argument(
        "--fasta",
        type=Path,
        default=Path("data/processed/immuno/pda/pda_designed.fasta"),
    )
    cohort_p.add_argument("--backend", choices=["mhcflurry", "netmhcpan"], default="mhcflurry")
    cohort_p.add_argument("--alleles-i", type=str, default=None, help="comma-separated HLA-I")
    cohort_p.add_argument("--out", type=Path, default=None)
    cohort_p.add_argument("--chunk-size", type=int, default=50_000)
    cohort_p.add_argument(
        "--binder-only",
        action="store_true",
        help="write only binder rows (predictions still run on all peptides)",
    )

    cleave_p = sub.add_parser(
        "cleave-cohort",
        help="Scan FASTA for protease cleavage sites (subsite-aware catalog)",
    )
    cleave_p.add_argument(
        "--fasta",
        type=Path,
        default=Path("data/processed/immuno/pda/pda_designed.fasta"),
    )
    cleave_p.add_argument("--out", type=Path, default=None)
    cleave_p.add_argument("--site-ids", type=str, default=None, help="comma-separated site ids")

    join_p = sub.add_parser(
        "join-epitope-protease",
        help="Join MHC binders to cleavage events (creating/destroying)",
    )
    join_p.add_argument(
        "--mhc",
        type=Path,
        default=Path("data/processed/immuno/pda/mhc1_cohort.csv"),
    )
    join_p.add_argument(
        "--cleavage",
        type=Path,
        default=Path("data/processed/immuno/pda/cleavage_cohort.csv"),
    )
    join_p.add_argument("--out", type=Path, default=None)
    join_p.add_argument("--shortlist", type=Path, default=None)
    join_p.add_argument("--max-rank", type=float, default=2.0)
    join_p.add_argument("--shortlist-n", type=int, default=40)

    bench_pull = sub.add_parser("benchling-pull", help="Pull AA sequences from Benchling")
    bench_pull.add_argument("--ids", type=str, default=None, help="comma-separated AA sequence IDs")
    bench_pull.add_argument("--name-includes", type=str, default=None)
    bench_pull.add_argument("--dry-run", action="store_true")

    bench_pub = sub.add_parser("benchling-publish", help="Publish a run summary to Benchling")
    bench_pub.add_argument("--run-dir", type=Path, required=True)
    bench_pub.add_argument("--dry-run", action="store_true")

    args = p.parse_args(argv)

    if args.cmd == "ensure-fixtures":
        paths = ensure_fixtures()
        print(json.dumps({k: str(v) for k, v in paths.items()}, indent=2))
        return 0

    if args.cmd == "fetch-refs":
        t = fetch_iedb_tcell_sample(limit=args.limit)
        a = fetch_atlas_stub()
        print(json.dumps({"tcell": str(t), "atlas": str(a)}, indent=2))
        return 0

    if args.cmd == "train":
        from re_agent.immuno_risk.model import train_and_save

        out = train_and_save()
        print(json.dumps({"model": str(out)}, indent=2))
        return 0

    if args.cmd == "pda-ingest":
        from re_agent.immuno_risk.pda import ingest_pda

        codes = [x.strip() for x in args.pdb_codes.split(",")] if args.pdb_codes else None
        summary = ingest_pda(
            source=args.source if not args.pda_json else "json",
            json_path=args.pda_json,
            out_dir=args.out,
            limit=args.limit,
            pdb_codes=codes,
            include_unknown=args.include_unknown,
            strip_his=not args.keep_his_tags,
            write_peptides=args.peptides,
        )
        print(json.dumps(summary, indent=2))
        return 0

    if args.cmd == "score-cohort":
        from re_agent.immuno_risk.batch import score_cohort

        alleles = [a.strip() for a in args.alleles_i.split(",")] if args.alleles_i else None
        summary = score_cohort(
            args.fasta,
            backend=args.backend,
            alleles=alleles,
            out_csv=args.out,
            chunk_size=args.chunk_size,
            binder_only=args.binder_only,
        )
        print(json.dumps(summary, indent=2))
        return 0

    if args.cmd == "cleave-cohort":
        from re_agent.immuno_risk.cleavage import scan_cohort_cleavage

        site_ids = [x.strip() for x in args.site_ids.split(",")] if args.site_ids else None
        summary = scan_cohort_cleavage(args.fasta, out_csv=args.out, site_ids=site_ids)
        print(json.dumps(summary, indent=2))
        return 0

    if args.cmd == "join-epitope-protease":
        from re_agent.immuno_risk.cleavage import join_epitope_protease

        summary = join_epitope_protease(
            args.mhc,
            args.cleavage,
            out_csv=args.out,
            shortlist_csv=args.shortlist,
            max_rank=args.max_rank,
            shortlist_n=args.shortlist_n,
        )
        print(json.dumps(summary, indent=2))
        return 0

    if args.cmd == "benchling-pull":
        from re_agent.immuno_risk.benchling import pull_candidates

        ids = [x.strip() for x in args.ids.split(",")] if args.ids else None
        rows = pull_candidates(ids=ids, name_includes=args.name_includes, dry_run=args.dry_run)
        print(json.dumps(rows, indent=2))
        return 0

    if args.cmd == "benchling-publish":
        from re_agent.immuno_risk.benchling import publish_run

        result = publish_run(args.run_dir, dry_run=args.dry_run)
        print(json.dumps(result, indent=2))
        return 0

    if args.cmd == "run":
        sid, seq = _read_sequence(args.sequence, args.fasta)
        if args.sequence_id:
            sid = args.sequence_id
        alleles_i = [a.strip() for a in args.alleles_i.split(",")] if args.alleles_i else None
        alleles_ii = [a.strip() for a in args.alleles_ii.split(",")] if args.alleles_ii else None
        result = run_immuno_risk(
            seq,
            sequence_id=sid,
            delivery_mode=args.delivery_mode,
            alleles_i=alleles_i,
            alleles_ii=alleles_ii,
            mhc_class=args.mhc_class,
            write=not args.no_write,
            use_netmhcpan=args.netmhcpan,
            use_netmhciipan=args.netmhciipan,
        )
        if args.json:
            print(result.model_dump_json(indent=2))
        else:
            print(
                json.dumps(
                    {
                        "run_id": result.run_id,
                        "overall": result.risk.overall,
                        "score0to100": result.risk.score0to100,
                        "confidence": result.confidence.score0to1,
                        "aggregation": result.aggregation.overall,
                        "artifact_dir": result.artifact_dir,
                        "caveats": result.caveats,
                        "versions": result.predictor_versions,
                    },
                    indent=2,
                )
            )
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
