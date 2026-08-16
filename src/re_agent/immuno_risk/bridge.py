"""JSON stdin/stdout bridge for TypeScript harness tools."""

from __future__ import annotations

import json
import sys
from typing import Any

from re_agent.immuno_risk.pipeline import run_immuno_risk


def main() -> int:
    raw = sys.stdin.read()
    if not raw.strip():
        print(json.dumps({"error": "expected JSON on stdin"}))
        return 1
    req: dict[str, Any] = json.loads(raw)
    action = req.get("action", "run")
    if action == "run":
        result = run_immuno_risk(
            req["sequence"],
            sequence_id=req.get("sequence_id", "query"),
            delivery_mode=req.get("delivery_mode", "intracellular_plasmid"),
            alleles_i=req.get("alleles_i"),
            alleles_ii=req.get("alleles_ii"),
            mhc_class=req.get("mhc_class", "both"),
            write=req.get("write", True),
            use_netmhcpan=bool(req.get("use_netmhcpan")),
            use_netmhciipan=bool(req.get("use_netmhciipan")),
        )
        print(result.model_dump_json())
        return 0
    if action == "aggregation":
        from re_agent.immuno_risk.aggregation import aggregation_report

        print(aggregation_report(req.get("sequence_id", "query"), req["sequence"]).model_dump_json())
        return 0
    print(json.dumps({"error": f"unknown action {action}"}))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
