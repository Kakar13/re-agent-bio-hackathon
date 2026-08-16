"""IEDB next-generation Tools API client.

Two distinct uses, and the distinction is the point of the whole pipeline:

* ``netchop_scores`` runs the proteasome cleavage predictor on a *sequence*, with
  no binding predictor attached. This resolves cleavage sites.
* ``score_peptides`` submits an *explicit peptide list* produced by the digest.
  IEDB therefore only ever scores peptides the digest says a protease could make,
  rather than every window in the sequence.

Every response is cached on disk under a hash of the request, so re-runs cost
nothing and any number in the report traces back to a stored payload.
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import httpx

from .paths import IEDB_CACHE, ensure_caches

log = logging.getLogger(__name__)

API_ROOT = "https://api-nextgen-tools.iedb.org/api/v1"
PIPELINE_URL = f"{API_ROOT}/pipeline"
RESULTS_URL = f"{API_ROOT}/results/{{result_id}}"

# Common HLA class I supertype representatives; broad population coverage
# without paying for a large allele panel on every call.
DEFAULT_MHCI_ALLELES = (
    "HLA-A*01:01,HLA-A*02:01,HLA-A*03:01,HLA-A*24:02,"
    "HLA-B*07:02,HLA-B*08:01,HLA-B*15:01,HLA-B*40:01"
)

# Standard reference DRB1 panel.
DEFAULT_MHCII_ALLELES = (
    "HLA-DRB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:01,"
    "HLA-DRB1*07:01,HLA-DRB1*11:01,HLA-DRB1*15:01"
)


class IEDBError(RuntimeError):
    pass


@dataclass
class IEDBTable:
    """A single result table plus the column names needed to read it."""

    table_type: str
    columns: list[str]
    rows: list[list[Any]]

    def as_dicts(self) -> list[dict]:
        return [dict(zip(self.columns, r)) for r in self.rows]


def _digest(payload: dict) -> str:
    blob = json.dumps(payload, sort_keys=True).encode()
    return hashlib.sha1(blob).hexdigest()[:16]


class IEDBClient:
    def __init__(
        self,
        *,
        cache_dir: Path | None = None,
        timeout: float = 120.0,
        poll_interval: float = 4.0,
        max_wait: float = 900.0,
        min_interval: float = 1.0,
    ) -> None:
        ensure_caches()
        self.cache_dir = cache_dir or IEDB_CACHE
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.poll_interval = poll_interval
        self.max_wait = max_wait
        self.min_interval = min_interval
        self._last_call = 0.0

    def _throttle(self) -> None:
        delta = time.time() - self._last_call
        if delta < self.min_interval:
            time.sleep(self.min_interval - delta)
        self._last_call = time.time()

    def _submit(self, payload: dict) -> str:
        backoff = 5.0
        for attempt in range(5):
            self._throttle()
            with httpx.Client(timeout=self.timeout) as c:
                r = c.post(PIPELINE_URL, json=payload)
            if r.status_code == 429 or r.status_code >= 500:
                log.warning("IEDB %s on submit, backing off %.0fs", r.status_code, backoff)
                time.sleep(backoff)
                backoff *= 2
                continue
            if r.status_code >= 400:
                raise IEDBError(f"IEDB submit failed {r.status_code}: {r.text[:300]}")
            data = r.json()
            rid = data.get("result_id")
            if not rid:
                raise IEDBError(f"no result_id in response: {json.dumps(data)[:300]}")
            for w in data.get("warnings") or []:
                log.warning("IEDB warning: %s", w)
            return rid
        raise IEDBError(f"IEDB submit failed after {attempt + 1} attempts")

    def _fetch(self, result_id: str) -> dict:
        """Poll until the job actually finishes.

        A pending job still returns HTTP 200 with a ``results`` array holding the
        echoed ``input_sequence_table``, so presence of results is not a
        completion signal. Only ``status`` is.
        """
        deadline = time.time() + self.max_wait
        while True:
            self._throttle()
            with httpx.Client(timeout=self.timeout) as c:
                r = c.get(RESULTS_URL.format(result_id=result_id), params={"format": "json"})
            if r.status_code in (202, 404):
                pass  # not ready yet
            elif r.status_code >= 400:
                raise IEDBError(f"IEDB results {r.status_code}: {r.text[:300]}")
            else:
                data = r.json()
                payload = data.get("data") or {}
                status = (data.get("status") or payload.get("status") or "").lower()
                if status in {"error", "failed"}:
                    errs = payload.get("errors") or payload.get("warnings") or []
                    raise IEDBError(f"IEDB job {result_id} failed: {json.dumps(errs)[:400]}")
                if status in {"done", "complete", "completed", "success", "finished"}:
                    return data
                if not status and self._has_prediction_table(payload):
                    return data
            if time.time() > deadline:
                raise IEDBError(f"IEDB job {result_id} timed out after {self.max_wait}s")
            time.sleep(self.poll_interval)

    @staticmethod
    def _has_prediction_table(payload: dict) -> bool:
        return any(
            res.get("type") not in (None, "input_sequence_table")
            for res in payload.get("results", [])
        )

    def run(self, payload: dict, *, tag: str, force: bool = False) -> dict:
        """Submit a pipeline payload, honouring the on-disk cache."""
        key = f"{tag}_{_digest(payload)}"
        cached = self.cache_dir / f"{key}.json"
        if cached.exists() and not force:
            log.debug("IEDB cache hit %s", cached.name)
            return json.loads(cached.read_text())

        log.info("IEDB submit [%s]", tag)
        rid = self._submit(payload)
        data = self._fetch(rid)
        record = {"request": payload, "result_id": rid, "response": data}
        cached.write_text(json.dumps(record))
        return record

    @staticmethod
    def tables(record: dict) -> list[IEDBTable]:
        out: list[IEDBTable] = []
        for res in record.get("response", {}).get("data", {}).get("results", []):
            cols = [c["name"] for c in res.get("table_columns", [])]
            out.append(
                IEDBTable(
                    table_type=res.get("type", "unknown"),
                    columns=cols,
                    rows=res.get("table_data", []),
                )
            )
        return out

    # ---------------------------------------------------------------- cleavage

    def netchop_scores(
        self,
        sequences: dict[str, str],
        *,
        network_method: str = "c_term",
        threshold: float = 0.5,
        force: bool = False,
    ) -> dict[str, list[float]]:
        """Per-residue proteasome cleavage probability, no binding predictor.

        ``c_term`` is the network trained on C-terminal cleavage, which is the
        relevant one: the proteasome defines the MHC class I peptide C-terminus
        while ERAP1 trims the N-terminus afterwards.
        """
        names = list(sequences)
        fasta = "".join(f">{n}\n{sequences[n]}\n" for n in names)
        payload = {
            "pipeline_title": "netchop-cleavage",
            "run_stage_range": [1, 1],
            "stages": [
                {
                    "stage_number": 1,
                    "tool_group": "mhci",
                    "input_sequence_text": fasta,
                    "input_parameters": {
                        "alleles": "HLA-A*02:01",
                        "peptide_length_range": [9, 9],
                        "predictors": [
                            {
                                "type": "processing",
                                "method": "netchop",
                                "network_method": network_method,
                                "threshold": threshold,
                            }
                        ],
                    },
                }
            ],
        }
        record = self.run(payload, tag="netchop", force=force)

        scores: dict[str, list[float]] = {n: [] for n in names}
        for table in self.tables(record):
            if table.table_type != "residue_table":
                continue
            for row in table.rows:
                if len(row) < 4:
                    continue
                seq_num, _pos, _aa, score = row[0], row[1], row[2], row[3]
                idx = int(seq_num) - 1
                if 0 <= idx < len(names):
                    scores[names[idx]].append(float(score))
        for n, s in scores.items():
            if s and len(s) != len(sequences[n]):
                log.warning("netchop returned %d scores for %s of length %d", len(s), n, len(sequences[n]))
        return scores

    # -------------------------------------------------------------- presentation

    def score_peptides(
        self,
        peptides: list[str],
        *,
        tool_group: str,
        alleles: str | None = None,
        predictors: list[dict] | None = None,
        batch_size: int = 400,
        force: bool = False,
    ) -> list[dict]:
        """Score an explicit peptide list.

        ``peptide_length_range: null`` tells IEDB to take each line as a whole
        peptide rather than sliding a window over it, which is what keeps the
        digest in control of which peptides exist.
        """
        if not peptides:
            return []
        if tool_group == "mhci":
            alleles = alleles or DEFAULT_MHCI_ALLELES
            # Binding only by default. Bundling the immunogenicity predictor into
            # the same stage makes the job fail, so it is requested separately
            # when wanted rather than silently taking the binding call down.
            predictors = predictors or [{"type": "binding", "method": "netmhcpan_el"}]
        elif tool_group == "mhcii":
            alleles = alleles or DEFAULT_MHCII_ALLELES
            predictors = predictors or [{"type": "binding", "method": "netmhciipan_el"}]
        else:
            raise ValueError(f"unsupported tool_group {tool_group!r}")

        unique = sorted(set(peptides))
        rows: list[dict] = []
        for i in range(0, len(unique), batch_size):
            batch = unique[i : i + batch_size]
            payload = {
                "pipeline_title": f"{tool_group}-peptides",
                "run_stage_range": [1, 1],
                "stages": [
                    {
                        "stage_number": 1,
                        "stage_type": "prediction",
                        "tool_group": tool_group,
                        "input_sequence_text": "\n".join(batch),
                        "input_parameters": {
                            "alleles": alleles,
                            "peptide_length_range": None,
                            "predictors": predictors,
                        },
                    }
                ],
            }
            record = self.run(payload, tag=f"{tool_group}-pep", force=force)
            for table in self.tables(record):
                if table.table_type == "residue_table":
                    continue
                rows.extend(table.as_dicts())
            log.info("%s: scored batch %d-%d of %d", tool_group, i, i + len(batch), len(unique))
        return rows
