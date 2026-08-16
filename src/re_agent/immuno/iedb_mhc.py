"""Real NetMHCIIpan EL/BA provider through the documented IEDB Tools API."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import time
from pathlib import Path

import httpx

from re_agent.immuno.adapters import sequence_sha256
from re_agent.immuno.contracts import (
    MHCHit,
    MHCProviderResult,
    MHCRequest,
    Provenance,
)

IEDB_MHCII_API = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"


class IEDBNetMHCIIpanProvider:
    """Call NetMHCIIpan 4.3 EL and BA separately, then cache their joined rows."""

    provider_id = "netmhciipan"

    def __init__(
        self,
        cache_dir: Path,
        version: str = "4.3",
        endpoint: str = IEDB_MHCII_API,
        timeout_seconds: float = 180.0,
    ) -> None:
        self.cache_dir = cache_dir
        self.version = version
        self.endpoint = endpoint
        self.timeout_seconds = timeout_seconds

    def predict(self, request: MHCRequest) -> MHCProviderResult:
        started = time.perf_counter()
        key = self._cache_key(request)
        cache_path = self.cache_dir / f"{key}.json"
        cached = cache_path.exists()
        warnings = [
            "EL presentation and BA binding are independent channels; raw values are not averaged."
        ]
        try:
            if cached:
                payload = json.loads(cache_path.read_text())
            else:
                payload = self._fetch(request)
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                cache_path.write_text(json.dumps(payload, indent=2) + "\n")
            hits = [MHCHit.model_validate(row) for row in payload["hits"]]
            supported = sorted({hit.allele for hit in hits})
            missing = sorted(set(request.alleles) - set(supported))
            if missing:
                warnings.append(f"provider returned no rows for {len(missing)} requested alleles")
            return MHCProviderResult(
                provider_id=self.provider_id,
                status="ok" if hits else "unsupported",
                hits=hits,
                supported_alleles=supported,
                missing_alleles=missing,
                provenance=self._provenance(
                    request,
                    runtime_seconds=time.perf_counter() - started,
                    cached=cached,
                ),
                warnings=warnings,
            )
        except Exception as exc:  # external provider boundary
            return MHCProviderResult(
                provider_id=self.provider_id,
                status="error",
                missing_alleles=request.alleles,
                provenance=self._provenance(
                    request,
                    runtime_seconds=time.perf_counter() - started,
                    cached=cached,
                ),
                warnings=warnings,
                error=f"{type(exc).__name__}: {exc}",
            )

    def _fetch(self, request: MHCRequest) -> dict:
        allele_csv = ",".join(request.alleles)
        el = self._call("netmhciipan_el", request.parent_sequence, allele_csv)
        ba = self._call("netmhciipan_ba", request.parent_sequence, allele_csv)
        joined: dict[tuple[str, int, int, str, str], dict] = {}
        for row in el:
            key = self._row_key(row)
            joined[key] = {
                "allele": row["allele"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "peptide": row["peptide"],
                "core": row["core_peptide"],
                "el_score": float(row["score"]),
                "el_rank": float(row["rank"]),
            }
        for row in ba:
            key = self._row_key(row)
            record = joined.setdefault(
                key,
                {
                    "allele": row["allele"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "peptide": row["peptide"],
                    "core": row["core_peptide"],
                },
            )
            record["ba_ic50_nm"] = float(row["ic50"])
            record["ba_rank"] = float(row["rank"])
        return {
            "schema_version": "1.0.0",
            "provider": self.provider_id,
            "version": self.version,
            "input_sha256": sequence_sha256(request.parent_sequence),
            "alleles": request.alleles,
            "hits": list(joined.values()),
        }

    def _call(self, method: str, sequence: str, alleles: str) -> list[dict[str, str]]:
        response = httpx.post(
            self.endpoint,
            data={
                "method": f"{method}-{self.version}",
                "sequence_text": sequence,
                "allele": alleles,
            },
            timeout=self.timeout_seconds,
        )
        response.raise_for_status()
        if response.text.lstrip().startswith("<"):
            raise RuntimeError("IEDB returned HTML instead of a prediction table")
        rows = list(csv.DictReader(io.StringIO(response.text), delimiter="\t"))
        if not rows:
            raise RuntimeError(f"IEDB returned no rows for {method}")
        return rows

    @staticmethod
    def _row_key(row: dict[str, str]) -> tuple[str, int, int, str, str]:
        return (
            row["allele"],
            int(row["start"]),
            int(row["end"]),
            row["peptide"],
            row["core_peptide"],
        )

    def _cache_key(self, request: MHCRequest) -> str:
        payload = json.dumps(
            {
                "provider": self.provider_id,
                "version": self.version,
                "sequence": request.parent_sequence,
                "alleles": sorted(request.alleles),
                "channels": ["el", "ba"],
            },
            sort_keys=True,
        )
        return hashlib.sha256(payload.encode()).hexdigest()

    def _provenance(
        self,
        request: MHCRequest,
        runtime_seconds: float,
        cached: bool,
    ) -> Provenance:
        return Provenance(
            provider=self.provider_id,
            version=f"{self.version} via IEDB Tools API",
            capability="mhc_evidence",
            source=self.endpoint,
            parameters={
                "el_method": f"netmhciipan_el-{self.version}",
                "ba_method": f"netmhciipan_ba-{self.version}",
                "alleles": request.alleles,
            },
            input_sha256=sequence_sha256(request.parent_sequence),
            runtime_seconds=runtime_seconds,
            cached=cached,
        )
