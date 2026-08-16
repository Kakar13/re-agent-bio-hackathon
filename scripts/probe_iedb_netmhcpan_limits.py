"""Probe IEDB tools_api/mhci NetMHCpan rate / size limits.

Stops on hard errors (429/503/timeouts). Does NOT attempt the full PDA cohort.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
METHOD = "netmhcpan_el-4.1"
ALLELE = "HLA-A*02:01"
LENGTH = "9"
OUT = Path("results/immuno_risk/iedb_rate_probe.json")

# Small synthetic 9mers (valid AA only)
PEPTIDES = [
    "GILGFVFTL",
    "NLVPMVATV",
    "GLCTLVAML",
    "ILKEPVHGV",
    "FLPSDFFPSV",
    "KTWGQYWQV",
    "RMFPNAPYL",
    "SLYNTVATL",
    "YLVAYQATV",
    "CLGGLLTMV",
]


def one_request(peptides: list[str], timeout: float = 120.0) -> dict:
    if len(peptides) == 1:
        seq = peptides[0]
    else:
        parts = []
        for i, p in enumerate(peptides):
            parts.append(f">p{i}\n{p}")
        seq = "\n".join(parts)
    data = urllib.parse.urlencode(
        {
            "method": METHOD,
            "sequence_text": seq,
            "allele": ALLELE,
            "length": LENGTH,
        }
    ).encode()
    req = urllib.request.Request(URL, data=data, method="POST")
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode("utf-8", errors="replace")
            return {
                "ok": True,
                "status": resp.status,
                "elapsed_s": round(time.time() - t0, 3),
                "bytes": len(body),
                "n_peptides": len(peptides),
                "n_lines": body.count("\n"),
                "snippet": body.splitlines()[:3],
            }
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")[:500]
        return {
            "ok": False,
            "status": e.code,
            "elapsed_s": round(time.time() - t0, 3),
            "n_peptides": len(peptides),
            "error": body,
        }
    except Exception as e:  # noqa: BLE001
        return {
            "ok": False,
            "status": None,
            "elapsed_s": round(time.time() - t0, 3),
            "n_peptides": len(peptides),
            "error": f"{type(e).__name__}: {e}",
        }


def request_with_poll(
    peptides: list[str],
    *,
    poll_s: float = 3.0,
    max_wait_s: float = 120.0,
    timeout: float = 120.0,
) -> dict:
    """Retry ``one_request`` every ``poll_s`` seconds until success or timeout."""
    deadline = time.time() + max_wait_s
    attempt = 0
    last: dict = {}
    while True:
        attempt += 1
        last = one_request(peptides, timeout=timeout)
        last["attempts"] = attempt
        if last["ok"]:
            return last
        remaining = deadline - time.time()
        if remaining <= 0:
            return last
        time.sleep(min(poll_s, remaining))


def phase_sequential(n: int = 20, sleep_s: float = 0.0, *, poll: bool = False) -> dict:
    """Fire n single-peptide requests sequentially (optional sleep / poll)."""
    results = []
    for i in range(n):
        pep = PEPTIDES[i % len(PEPTIDES)]
        r = request_with_poll([pep], poll_s=3.0, max_wait_s=60.0) if poll else one_request([pep])
        results.append(r)
        print(
            f"  seq[{i+1}/{n}] status={r.get('status')} ok={r['ok']} "
            f"{r['elapsed_s']}s attempts={r.get('attempts', 1)}",
            flush=True,
        )
        if not r["ok"] and r.get("status") in {429, 503, 502, 500} and not poll:
            break
        if sleep_s:
            time.sleep(sleep_s)
    return {"phase": "sequential", "sleep_s": sleep_s, "poll": poll, "n_attempted": n, "results": results}


def phase_concurrent(concurrency: int, n: int) -> dict:
    """Fire n single-peptide requests with given concurrency."""
    results = []
    print(f"  concurrent workers={concurrency} jobs={n}", flush=True)
    with ThreadPoolExecutor(max_workers=concurrency) as ex:
        futs = [ex.submit(one_request, [PEPTIDES[i % len(PEPTIDES)]]) for i in range(n)]
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            results.append(r)
            print(
                f"  conc[{i}/{n}] status={r.get('status')} ok={r['ok']} {r['elapsed_s']}s",
                flush=True,
            )
    fails = [r for r in results if not r["ok"]]
    return {
        "phase": "concurrent",
        "concurrency": concurrency,
        "n": n,
        "n_ok": sum(1 for r in results if r["ok"]),
        "n_fail": len(fails),
        "fail_statuses": sorted({r.get("status") for r in fails}),
        "results": results,
    }


def phase_batch_sizes(sizes: list[int]) -> dict:
    """Grow peptides-per-request until failure."""
    results = []
    for size in sizes:
        # Repeat base peptides to reach size
        peptides = [PEPTIDES[i % len(PEPTIDES)] for i in range(size)]
        # Make them unique-ish by mutating last AA cyclically so server doesn't dedupe
        aa = "ACDEFGHIKLMNPQRSTVWY"
        peptides = [p[:8] + aa[i % len(aa)] for i, p in enumerate(peptides)]
        print(f"  batch size={size} …", flush=True)
        r = one_request(peptides, timeout=300.0)
        results.append(r)
        print(
            f"  batch[{size}] status={r.get('status')} ok={r['ok']} {r['elapsed_s']}s "
            f"lines={r.get('n_lines')}",
            flush=True,
        )
        if not r["ok"]:
            break
    return {"phase": "batch_size", "results": results}


def main() -> int:
    report: dict = {"url": URL, "method": METHOD, "phases": []}
    OUT.parent.mkdir(parents=True, exist_ok=True)

    print("=== phase 1: sequential ×20, no sleep ===", flush=True)
    report["phases"].append(phase_sequential(20, sleep_s=0.0))

    print("=== phase 2: concurrent ramp ===", flush=True)
    for conc in (2, 4, 8, 16):
        report["phases"].append(phase_concurrent(conc, n=conc * 2))
        # brief pause between ramps to be polite
        time.sleep(2)

    print("=== phase 3: batch size ramp ===", flush=True)
    report["phases"].append(phase_batch_sizes([1, 10, 50, 100, 250, 500, 1000, 2000, 5000]))

    # Summary
    summary = {"hard_errors": [], "batch_max_ok": None, "concurrent_max_all_ok": None}
    for phase in report["phases"]:
        if phase["phase"] == "batch_size":
            oks = [r for r in phase["results"] if r["ok"]]
            if oks:
                summary["batch_max_ok"] = max(r["n_peptides"] for r in oks)
            fails = [r for r in phase["results"] if not r["ok"]]
            for f in fails:
                summary["hard_errors"].append({"phase": "batch", **{k: f.get(k) for k in ("status", "error", "n_peptides")}})
        if phase["phase"] == "concurrent":
            if phase["n_fail"] == 0:
                summary["concurrent_max_all_ok"] = phase["concurrency"]
            else:
                summary["hard_errors"].append(
                    {
                        "phase": "concurrent",
                        "concurrency": phase["concurrency"],
                        "n_fail": phase["n_fail"],
                        "fail_statuses": phase["fail_statuses"],
                    }
                )
        if phase["phase"] == "sequential":
            fails = [r for r in phase["results"] if not r["ok"]]
            for f in fails:
                summary["hard_errors"].append({"phase": "sequential", **{k: f.get(k) for k in ("status", "error")}})

    report["summary"] = summary
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    print(f"wrote {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
