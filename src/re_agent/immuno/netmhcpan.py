"""NetMHCpan 4.1 teacher labeling through the documented IEDB Tools API.

The EL and BA channels are retained separately. EL approximates presentation,
while BA is the binding-only channel that can replace the original MHCflurry
binding head without silently double-counting cleavage and TAP evidence.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import httpx
import pandas as pd

IEDB_MHCI_API = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
NETMHCPAN_VERSION = "4.1"
DEFAULT_ALLELE = "HLA-A*02:01"
DEFAULT_PEPTIDE_LENGTH = 9
CHANNELS = ("el", "ba")


def _normalized_score(row: dict[str, str]) -> float:
    """Return a bounded higher-is-better score for either API channel."""

    if row.get("score"):
        return float(row["score"])
    if row.get("ic50"):
        affinity = max(float(row["ic50"]), 1.0)
        return max(0.0, min(1.0, 1.0 - math.log(affinity) / math.log(50_000.0)))
    raise RuntimeError(f"NetMHCpan row has neither score nor ic50: columns={sorted(row)}")


@dataclass(frozen=True)
class NetMHCpanTeacherConfig:
    allele: str = DEFAULT_ALLELE
    peptide_length: int = DEFAULT_PEPTIDE_LENGTH
    version: str = NETMHCPAN_VERSION
    batch_size: int = 500
    timeout_seconds: float = 300.0
    max_retries: int = 3


class NetMHCpanTeacher:
    """Label fixed-length peptides with separate NetMHCpan EL and BA outputs."""

    def __init__(
        self,
        cache_dir: Path,
        config: NetMHCpanTeacherConfig | None = None,
        endpoint: str = IEDB_MHCI_API,
    ) -> None:
        self.cache_dir = Path(cache_dir)
        self.config = config or NetMHCpanTeacherConfig()
        self.endpoint = endpoint

    def label(self, frame: pd.DataFrame, peptide_column: str = "peptide") -> pd.DataFrame:
        """Return ``frame`` with NetMHCpan columns attached.

        Duplicate peptides are submitted only once and mapped back to every
        source row. Cache keys include the ordered peptide batch and all model
        parameters, so interrupted jobs can resume without relabeling completed
        batches.
        """

        if peptide_column not in frame:
            raise ValueError(f"missing peptide column {peptide_column!r}")
        peptides = frame[peptide_column].astype(str).str.upper()
        invalid = peptides.str.len().ne(self.config.peptide_length)
        if invalid.any():
            sample = peptides[invalid].head(3).tolist()
            raise ValueError(
                f"NetMHCpan cohort must contain {self.config.peptide_length}-mers; "
                f"found {int(invalid.sum())} invalid rows such as {sample}"
            )
        unique = list(dict.fromkeys(peptides))
        tables = []
        for channel in CHANNELS:
            channel_tables = []
            for offset in range(0, len(unique), self.config.batch_size):
                batch = unique[offset : offset + self.config.batch_size]
                channel_tables.append(self._label_batch(batch, channel))
            channel_frame = pd.concat(channel_tables, ignore_index=True)
            channel_frame = channel_frame.rename(
                columns={
                    "score": f"netmhcpan_{channel}_score",
                    "percentile_rank": f"netmhcpan_{channel}_rank",
                    "ic50": f"netmhcpan_{channel}_ic50_nm",
                }
            )
            keep = [
                column
                for column in (
                    "peptide",
                    f"netmhcpan_{channel}_score",
                    f"netmhcpan_{channel}_rank",
                    f"netmhcpan_{channel}_ic50_nm",
                )
                if column in channel_frame
            ]
            tables.append(channel_frame[keep])

        labels = tables[0].merge(tables[1], on="peptide", how="outer", validate="one_to_one")
        out = frame.copy()
        out[peptide_column] = peptides
        out = out.merge(
            labels,
            left_on=peptide_column,
            right_on="peptide",
            how="left",
            validate="many_to_one",
            suffixes=("", "_teacher"),
        )
        if peptide_column != "peptide":
            out = out.drop(columns=["peptide_teacher"])
        out["netmhcpan_allele"] = self.config.allele
        out["netmhcpan_version"] = self.config.version
        out["netmhcpan_source"] = self.endpoint
        return out

    def label_parent_sequences(
        self,
        parents: pd.DataFrame,
        *,
        id_column: str = "parent_sequence_id",
        sequence_column: str = "sequence",
        parents_per_batch: int = 20,
    ) -> pd.DataFrame:
        """Tile and label whole proteins in NetMHCpan itself.

        Sending whole parents avoids one expensive API startup per peptide and
        guarantees that every possible 9-mer is considered before any sampling.
        IEDB coordinates are converted from 1-indexed inclusive bounds to this
        project's 0-indexed half-open convention.
        """

        missing = {id_column, sequence_column} - set(parents)
        if missing:
            raise ValueError(f"parent table is missing columns: {sorted(missing)}")
        if parents[id_column].duplicated().any():
            raise ValueError(f"{id_column} must be unique")

        parent_rows = parents.to_dict(orient="records")
        channel_tables: list[pd.DataFrame] = []
        for channel in CHANNELS:
            batches = []
            for offset in range(0, len(parent_rows), parents_per_batch):
                batch = parent_rows[offset : offset + parents_per_batch]
                batches.append(
                    self._label_parent_batch(
                        batch,
                        channel,
                        id_column=id_column,
                        sequence_column=sequence_column,
                    )
                )
            channel_frame = pd.concat(batches, ignore_index=True)
            channel_frame = channel_frame.rename(
                columns={
                    "score": f"netmhcpan_{channel}_score",
                    "percentile_rank": f"netmhcpan_{channel}_rank",
                    "ic50": f"netmhcpan_{channel}_ic50_nm",
                }
            )
            channel_tables.append(channel_frame)

        keys = [id_column, "start", "end", "peptide"]
        labels = channel_tables[0].merge(
            channel_tables[1],
            on=keys,
            how="outer",
            validate="one_to_one",
            suffixes=("_el", "_ba"),
        )
        parent_metadata = parents.drop(columns=[sequence_column])
        out = labels.merge(parent_metadata, on=id_column, how="left", validate="many_to_one")
        sequence_lookup = parents.set_index(id_column)[sequence_column].astype(str).to_dict()
        flank = 4
        out["n_flank"] = [
            sequence_lookup[parent_id][max(0, start - flank) : start]
            for parent_id, start in zip(out[id_column], out["start"], strict=True)
        ]
        out["c_flank"] = [
            sequence_lookup[parent_id][end : end + flank]
            for parent_id, end in zip(out[id_column], out["end"], strict=True)
        ]
        out["length"] = out["peptide"].str.len()
        out["hla_allele"] = self.config.allele
        out["netmhcpan_allele"] = self.config.allele
        out["netmhcpan_version"] = self.config.version
        out["netmhcpan_source"] = self.endpoint
        return out

    def _label_batch(self, peptides: list[str], channel: str) -> pd.DataFrame:
        cache_path = self._cache_path(peptides, channel)
        if cache_path.exists():
            rows = json.loads(cache_path.read_text())["rows"]
            return pd.DataFrame(rows)

        fasta = "".join(f">peptide_{index}\n{peptide}\n" for index, peptide in enumerate(peptides))
        method = f"netmhcpan_{channel}-{self.config.version}"
        response_text = self._post_with_retries(
            {
                "method": method,
                "sequence_text": fasta,
                "allele": self.config.allele,
                "length": str(self.config.peptide_length),
            }
        )
        if response_text.lstrip().startswith("<"):
            raise RuntimeError("IEDB returned HTML instead of a NetMHCpan prediction table")
        rows = list(csv.DictReader(io.StringIO(response_text), delimiter="\t"))
        if len(rows) != len(peptides):
            raise RuntimeError(
                f"{method} returned {len(rows)} rows for {len(peptides)} fixed-length peptides"
            )

        normalized = []
        for row in rows:
            seq_index = int(row["seq_num"]) - 1
            if not 0 <= seq_index < len(peptides):
                raise RuntimeError(f"{method} returned invalid seq_num={row['seq_num']}")
            record: dict[str, str | float] = {
                "peptide": peptides[seq_index],
                "score": _normalized_score(row),
                "percentile_rank": float(row["percentile_rank"]),
            }
            if row.get("ic50"):
                record["ic50"] = float(row["ic50"])
            normalized.append(record)

        payload = {
            "schema_version": "1.0.0",
            "provider": "netmhcpan",
            "version": self.config.version,
            "channel": channel,
            "allele": self.config.allele,
            "peptide_length": self.config.peptide_length,
            "endpoint": self.endpoint,
            "rows": normalized,
        }
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(payload, indent=2) + "\n")
        return pd.DataFrame(normalized)

    def _label_parent_batch(
        self,
        parents: list[dict],
        channel: str,
        *,
        id_column: str,
        sequence_column: str,
    ) -> pd.DataFrame:
        cache_path = self._parent_cache_path(
            parents,
            channel,
            id_column=id_column,
            sequence_column=sequence_column,
        )
        if cache_path.exists():
            return pd.DataFrame(json.loads(cache_path.read_text())["rows"])

        fasta = "".join(
            f">parent_{index}\n{str(parent[sequence_column]).upper()}\n"
            for index, parent in enumerate(parents)
        )
        method = f"netmhcpan_{channel}-{self.config.version}"
        response_text = self._post_with_retries(
            {
                "method": method,
                "sequence_text": fasta,
                "allele": self.config.allele,
                "length": str(self.config.peptide_length),
            }
        )
        if response_text.lstrip().startswith("<"):
            raise RuntimeError("IEDB returned HTML instead of a NetMHCpan prediction table")
        rows = list(csv.DictReader(io.StringIO(response_text), delimiter="\t"))
        if not rows:
            raise RuntimeError(f"{method} returned no rows")

        normalized = []
        for row in rows:
            seq_index = int(row["seq_num"]) - 1
            if not 0 <= seq_index < len(parents):
                raise RuntimeError(f"{method} returned invalid seq_num={row['seq_num']}")
            record: dict[str, str | float | int] = {
                id_column: str(parents[seq_index][id_column]),
                "start": int(row["start"]) - 1,
                "end": int(row["end"]),
                "peptide": row["peptide"],
                "score": _normalized_score(row),
                "percentile_rank": float(row["percentile_rank"]),
            }
            if row.get("ic50"):
                record["ic50"] = float(row["ic50"])
            normalized.append(record)

        payload = {
            "schema_version": "1.0.0",
            "provider": "netmhcpan",
            "version": self.config.version,
            "channel": channel,
            "allele": self.config.allele,
            "peptide_length": self.config.peptide_length,
            "endpoint": self.endpoint,
            "rows": normalized,
        }
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(payload, indent=2) + "\n")
        return pd.DataFrame(normalized)

    def _post_with_retries(self, data: dict[str, str]) -> str:
        last_error: Exception | None = None
        for attempt in range(self.config.max_retries):
            try:
                response = httpx.post(
                    self.endpoint,
                    data=data,
                    timeout=self.config.timeout_seconds,
                )
                response.raise_for_status()
                return response.text
            except (httpx.HTTPError, TimeoutError) as exc:
                last_error = exc
                if attempt + 1 < self.config.max_retries:
                    time.sleep(2**attempt)
        raise RuntimeError(f"NetMHCpan request failed after retries: {last_error}") from last_error

    def _cache_path(self, peptides: list[str], channel: str) -> Path:
        digest = hashlib.sha256(
            json.dumps(
                {
                    "provider": "netmhcpan",
                    "version": self.config.version,
                    "channel": channel,
                    "allele": self.config.allele,
                    "peptide_length": self.config.peptide_length,
                    "peptides": peptides,
                },
                sort_keys=True,
            ).encode()
        ).hexdigest()
        return self.cache_dir / channel / f"{digest}.json"

    def _parent_cache_path(
        self,
        parents: list[dict],
        channel: str,
        *,
        id_column: str,
        sequence_column: str,
    ) -> Path:
        digest = hashlib.sha256(
            json.dumps(
                {
                    "provider": "netmhcpan",
                    "version": self.config.version,
                    "channel": channel,
                    "allele": self.config.allele,
                    "peptide_length": self.config.peptide_length,
                    "parents": [
                        {
                            "id": str(parent[id_column]),
                            "sequence": str(parent[sequence_column]).upper(),
                        }
                        for parent in parents
                    ],
                },
                sort_keys=True,
            ).encode()
        ).hexdigest()
        return self.cache_dir / "parents" / channel / f"{digest}.json"
