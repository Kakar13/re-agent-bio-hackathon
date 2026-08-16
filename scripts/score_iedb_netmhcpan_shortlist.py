"""Polite IEDB NetMHCpan-4.1 EL client with fixed-interval polling on failure.

On 403/429/5xx/network errors, wait ``poll_s`` (default 3s) and retry until
success or ``max_wait_s`` is exceeded. Single-flight only — no concurrency.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
DEFAULT_METHOD = "netmhcpan_el-4.1"


def _post(
    *,
    peptides: list[str],
    allele: str,
    length: str | None = None,
    method: str = DEFAULT_METHOD,
    timeout: float = 120.0,
) -> dict:
    if not peptides:
        raise ValueError("peptides required")
    if len(peptides) == 1:
        seq = peptides[0]
    else:
        seq = "\n".join(f">p{i}\n{p}" for i, p in enumerate(peptides))
    length = length or str(len(peptides[0]))
    data = urllib.parse.urlencode(
        {
            "method": method,
            "sequence_text": seq,
            "allele": allele,
            "length": length,
        }
    ).encode()
    req = urllib.request.Request(
        URL,
        data=data,
        method="POST",
        headers={"User-Agent": "re-agent-immuno-risk/0.1 (hackathon; polite single-flight)"},
    )
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode("utf-8", errors="replace")
            return {
                "ok": True,
                "status": resp.status,
                "elapsed_s": round(time.time() - t0, 3),
                "body": body,
                "n_peptides": len(peptides),
                "attempts": 1,
            }
    except urllib.error.HTTPError as e:
        err = e.read().decode("utf-8", errors="replace")[:500]
        return {
            "ok": False,
            "status": e.code,
            "elapsed_s": round(time.time() - t0, 3),
            "error": err,
            "n_peptides": len(peptides),
            "attempts": 1,
        }
    except Exception as e:  # noqa: BLE001
        return {
            "ok": False,
            "status": None,
            "elapsed_s": round(time.time() - t0, 3),
            "error": f"{type(e).__name__}: {e}",
            "n_peptides": len(peptides),
            "attempts": 1,
        }


def request_with_poll(
    *,
    peptides: list[str],
    allele: str,
    length: str | None = None,
    method: str = DEFAULT_METHOD,
    poll_s: float = 3.0,
    max_wait_s: float = 600.0,
    timeout: float = 120.0,
    verbose: bool = True,
) -> dict:
    """POST once; on failure, sleep ``poll_s`` and retry until ``max_wait_s``."""
    deadline = time.time() + max_wait_s
    attempt = 0
    last: dict = {}
    while True:
        attempt += 1
        last = _post(
            peptides=peptides,
            allele=allele,
            length=length,
            method=method,
            timeout=timeout,
        )
        last["attempts"] = attempt
        if last["ok"]:
            if verbose:
                print(
                    f"  ok allele={allele} n={len(peptides)} attempt={attempt} "
                    f"{last['elapsed_s']}s",
                    flush=True,
                )
            return last
        remaining = deadline - time.time()
        if remaining <= 0:
            if verbose:
                print(
                    f"  give up allele={allele} after {attempt} attempts "
                    f"status={last.get('status')} err={(last.get('error') or '')[:80]}",
                    flush=True,
                )
            return last
        if verbose:
            print(
                f"  retry in {poll_s:.0f}s (attempt={attempt} status={last.get('status')} "
                f"remaining={remaining:.0f}s)",
                flush=True,
            )
        time.sleep(min(poll_s, remaining))


def parse_tsv(body: str) -> list[dict]:
    lines = [ln for ln in body.splitlines() if ln.strip()]
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    return list(reader)


def score_pairs(
    pairs: list[tuple[str, str]],
    *,
    poll_s: float = 3.0,
    max_wait_s: float = 600.0,
    gap_s: float = 1.0,
) -> list[dict]:
    """Score (peptide, allele) pairs one at a time with polling on failure."""
    out: list[dict] = []
    for i, (peptide, allele) in enumerate(pairs, 1):
        print(f"[{i}/{len(pairs)}] {peptide} {allele}", flush=True)
        r = request_with_poll(
            peptides=[peptide],
            allele=allele,
            length=str(len(peptide)),
            poll_s=poll_s,
            max_wait_s=max_wait_s,
        )
        row: dict = {
            "peptide": peptide,
            "allele": allele,
            "ok": r["ok"],
            "status": r.get("status"),
            "attempts": r.get("attempts"),
            "elapsed_s": r.get("elapsed_s"),
        }
        if r["ok"]:
            parsed = parse_tsv(r["body"])
            # Prefer exact peptide match
            hit = next((p for p in parsed if p.get("peptide") == peptide), parsed[0] if parsed else {})
            row["score"] = hit.get("score")
            row["percentile_rank"] = hit.get("percentile_rank")
            row["core"] = hit.get("core")
        else:
            row["error"] = (r.get("error") or "")[:200]
        out.append(row)
        if i < len(pairs) and gap_s > 0:
            time.sleep(gap_s)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--shortlist",
        type=Path,
        default=Path("data/processed/immuno/pda/fold_shortlist.csv"),
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("results/immuno_risk/iedb_shortlist_netmhcpan.csv"),
    )
    p.add_argument("--poll-s", type=float, default=3.0, help="seconds between retries on failure")
    p.add_argument("--max-wait-s", type=float, default=600.0, help="max seconds to keep polling one request")
    p.add_argument("--gap-s", type=float, default=1.0, help="pause between successful requests")
    p.add_argument("--limit", type=int, default=None, help="max unique (peptide,allele) pairs")
    p.add_argument("--smoke", action="store_true", help="single GILGFVFTL / HLA-A*02:01 probe with polling")
    args = p.parse_args()

    if args.smoke:
        print("smoke: GILGFVFTL HLA-A*02:01 with 3s poll…", flush=True)
        r = request_with_poll(
            peptides=["GILGFVFTL"],
            allele="HLA-A*02:01",
            poll_s=args.poll_s,
            max_wait_s=args.max_wait_s,
        )
        print(json.dumps({k: v for k, v in r.items() if k != "body"}, indent=2))
        if r["ok"]:
            print(r["body"])
        return 0 if r["ok"] else 1

    if not args.shortlist.exists():
        raise SystemExit(f"shortlist not found: {args.shortlist}")

    pairs: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    with args.shortlist.open() as f:
        for row in csv.DictReader(f):
            key = (row["peptide"], row["allele"])
            if key in seen:
                continue
            seen.add(key)
            pairs.append(key)
            if args.limit and len(pairs) >= args.limit:
                break

    print(f"scoring {len(pairs)} unique (peptide, allele) pairs via IEDB…", flush=True)
    rows = score_pairs(
        pairs,
        poll_s=args.poll_s,
        max_wait_s=args.max_wait_s,
        gap_s=args.gap_s,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "peptide",
        "allele",
        "ok",
        "status",
        "attempts",
        "elapsed_s",
        "score",
        "percentile_rank",
        "core",
        "error",
    ]
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n_ok = sum(1 for r in rows if r["ok"])
    summary = {
        "n_pairs": len(rows),
        "n_ok": n_ok,
        "n_fail": len(rows) - n_ok,
        "out": str(args.out),
        "poll_s": args.poll_s,
    }
    args.out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0 if n_ok == len(rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
