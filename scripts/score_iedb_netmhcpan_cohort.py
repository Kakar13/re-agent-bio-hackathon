"""Score the full PDA unique-peptide cohort via IEDB NetMHCpan-4.1 EL.

Polite single-flight client: batches same-length peptides, retries every
``--poll-s`` seconds on failure, resumes from a checkpoint JSONL.

Output CSV always includes the peptide sequence.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

from score_iedb_netmhcpan_shortlist import DEFAULT_METHOD, parse_tsv, request_with_poll

DEFAULT_ALLELES = [
    "HLA-A*02:01",
    "HLA-A*01:01",
    "HLA-A*03:01",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-C*07:01",
]


def load_peptides(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def chunked(items: list[str], size: int):
    for i in range(0, len(items), size):
        yield i, items[i : i + size]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--peptides",
        type=Path,
        default=Path("data/processed/immuno/pda/netmhcpan_peptides.txt"),
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("results/immuno_risk/iedb_netmhcpan_cohort.csv"),
    )
    p.add_argument("--batch-size", type=int, default=100)
    p.add_argument("--poll-s", type=float, default=3.0)
    p.add_argument("--max-wait-s", type=float, default=900.0)
    p.add_argument("--gap-s", type=float, default=1.5, help="pause after each successful batch")
    p.add_argument("--alleles", type=str, default=",".join(DEFAULT_ALLELES))
    p.add_argument("--limit-peptides", type=int, default=None, help="optional cap for smoke")
    p.add_argument("--shuffle", action="store_true", help="shuffle peptides before limit (reproducible with --seed)")
    p.add_argument("--seed", type=int, default=42, help="RNG seed used when --shuffle is set")
    p.add_argument("--max-batches", type=int, default=None, help="stop after N batches (debug)")
    args = p.parse_args()

    alleles = [a.strip() for a in args.alleles.split(",") if a.strip()]
    peptides = load_peptides(args.peptides)
    if args.shuffle:
        import random

        rng = random.Random(args.seed)
        peptides = list(peptides)
        rng.shuffle(peptides)
        print(f"shuffled {len(peptides)} peptides with seed={args.seed}", flush=True)
    if args.limit_peptides:
        peptides = peptides[: args.limit_peptides]
        print(f"limited to {len(peptides)} peptides", flush=True)
    by_len: dict[int, list[str]] = {}
    for pep in peptides:
        by_len.setdefault(len(pep), []).append(pep)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    ckpt = args.out.with_suffix(".checkpoint.json")
    done: set[str] = set()
    if ckpt.exists():
        done = set(json.loads(ckpt.read_text()).get("done_keys", []))
        print(f"resuming; {len(done)} batches already done", flush=True)

    write_header = not args.out.exists() or args.out.stat().st_size == 0
    fields = [
        "peptide",
        "allele",
        "length",
        "score",
        "percentile_rank",
        "core",
        "icore",
        "start",
        "end",
        "method",
    ]

    n_rows = 0
    n_batches = 0
    n_fail = 0
    t0 = time.time()

    with args.out.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if write_header:
            w.writeheader()

        for allele in alleles:
            for length, pep_list in sorted(by_len.items()):
                for offset, batch in chunked(pep_list, args.batch_size):
                    key = f"{allele}|L{length}|off{offset}|n{len(batch)}"
                    if key in done:
                        continue
                    n_batches += 1
                    print(
                        f"[{n_batches}] {allele} len={length} offset={offset} "
                        f"batch={len(batch)} …",
                        flush=True,
                    )
                    r = request_with_poll(
                        peptides=batch,
                        allele=allele,
                        length=str(length),
                        method=DEFAULT_METHOD,
                        poll_s=args.poll_s,
                        max_wait_s=args.max_wait_s,
                        timeout=300.0,
                    )
                    if not r["ok"]:
                        n_fail += 1
                        print(f"  FAIL status={r.get('status')} attempts={r.get('attempts')}", flush=True)
                        # persist progress and stop so we can resume later
                        ckpt.write_text(
                            json.dumps(
                                {
                                    "done_keys": sorted(done),
                                    "last_fail": key,
                                    "n_rows": n_rows,
                                    "n_fail": n_fail,
                                },
                                indent=2,
                            )
                        )
                        summary = {
                            "status": "paused_on_fail",
                            "n_rows_written": n_rows,
                            "n_batches_ok": len(done),
                            "n_fail": n_fail,
                            "elapsed_s": round(time.time() - t0, 1),
                            "out": str(args.out),
                            "checkpoint": str(ckpt),
                        }
                        args.out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
                        print(json.dumps(summary, indent=2), flush=True)
                        return 2

                    parsed = parse_tsv(r["body"])
                    # Map by peptide; IEDB returns one row per peptide for fixed length
                    for row in parsed:
                        pep = row.get("peptide") or ""
                        if not pep:
                            continue
                        w.writerow(
                            {
                                "peptide": pep,
                                "allele": row.get("allele") or allele,
                                "length": len(pep),
                                "score": row.get("score"),
                                "percentile_rank": row.get("percentile_rank"),
                                "core": row.get("core"),
                                "icore": row.get("icore"),
                                "start": row.get("start"),
                                "end": row.get("end"),
                                "method": DEFAULT_METHOD,
                            }
                        )
                        n_rows += 1
                    f.flush()
                    done.add(key)
                    ckpt.write_text(
                        json.dumps({"done_keys": sorted(done), "n_rows": n_rows}, indent=2)
                    )
                    if args.max_batches and n_batches >= args.max_batches:
                        print("hit --max-batches; stopping for resume", flush=True)
                        break
                    if args.gap_s > 0:
                        time.sleep(args.gap_s)
                else:
                    continue
                break
            else:
                continue
            break

    summary = {
        "status": "complete" if n_fail == 0 and args.max_batches is None else "partial",
        "n_peptides": len(peptides),
        "n_alleles": len(alleles),
        "n_rows_written": n_rows,
        "n_batches_ok": len(done),
        "n_fail": n_fail,
        "elapsed_s": round(time.time() - t0, 1),
        "out": str(args.out),
        "batch_size": args.batch_size,
        "method": DEFAULT_METHOD,
        "note": "CSV peptide column is the epitope sequence scored by NetMHCpan-4.1 EL via IEDB.",
    }
    args.out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
