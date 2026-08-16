#!/usr/bin/env python3
"""Run the real adapter-based immunogenicity screening pipeline."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from re_agent.immuno.adapters import (  # noqa: E402
    ArtifactResponseAdapter,
    ChallengerArtifactProvider,
    NetMHCIIpanArtifactProvider,
    UnavailableResponseAdapter,
    sequence_sha256,
)
from re_agent.immuno.calibration import FrozenCalibratedAdapter  # noqa: E402
from re_agent.immuno.iedb_mhc import IEDBNetMHCIIpanProvider  # noqa: E402
from re_agent.immuno.pipeline import (  # noqa: E402
    ImmunogenicityScreeningAgent,
    write_assessments,
)
from re_agent.immuno.registry import MHCProviderRegistry, ResponseModelRegistry  # noqa: E402


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].split()[0]
            chunks = []
        elif name is None:
            raise ValueError("FASTA sequence appeared before its header")
        else:
            chunks.append(line.upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        raise ValueError("FASTA contains no records")
    return records


def parse_challenger(value: str) -> ChallengerArtifactProvider:
    try:
        provider_id, version, path = value.split("=", 2)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "challenger must be PROVIDER=VERSION=PATH"
        ) from exc
    return ChallengerArtifactProvider(provider_id, version, Path(path))


def parse_mapping(value: str) -> tuple[str, Path]:
    try:
        key, path = value.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be ADAPTER_ID=PATH") from exc
    return key, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    candidate_source = parser.add_mutually_exclusive_group(required=True)
    candidate_source.add_argument("--fasta", type=Path)
    candidate_source.add_argument("--sequence", help="Single raw amino-acid sequence")
    parser.add_argument("--sequence-id", default="query")
    parser.add_argument("--response-artifact", action="append", type=Path, default=[])
    parser.add_argument(
        "--response-unavailable",
        action="append",
        default=[],
        metavar="ADAPTER_ID",
        help="Declare an expected response adapter unavailable without fabricating scores",
    )
    parser.add_argument("--default-response-adapter", required=True)
    parser.add_argument("--calibration", action="append", default=[], type=parse_mapping)
    mhc_source = parser.add_mutually_exclusive_group(required=True)
    mhc_source.add_argument("--netmhciipan-artifact", type=Path)
    mhc_source.add_argument(
        "--iedb-live",
        action="store_true",
        help="Call documented IEDB NetMHCIIpan 4.3 EL and BA endpoints with disk cache",
    )
    parser.add_argument("--challenger", action="append", default=[], type=parse_challenger)
    parser.add_argument("--accessibility-artifact", type=Path)
    parser.add_argument("--cleavage-artifact", type=Path)
    parser.add_argument("--shared-hla-artifact", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "reports" / "candidate_screen.json",
    )
    args = parser.parse_args()

    calibrations = dict(args.calibration)
    artifact_adapters = [ArtifactResponseAdapter(path) for path in args.response_artifact]
    response_adapters = [
        *(
            FrozenCalibratedAdapter(adapter, calibrations[adapter.adapter_id])
            if adapter.adapter_id in calibrations
            else adapter
            for adapter in artifact_adapters
        ),
        *(
            UnavailableResponseAdapter(
                adapter_id,
                "checkpoint or immutable prediction artifact has not been handed off",
            )
            for adapter_id in args.response_unavailable
        ),
    ]
    if not response_adapters:
        parser.error("provide --response-artifact or --response-unavailable")
    response_registry = ResponseModelRegistry(response_adapters)
    primary_mhc = (
        NetMHCIIpanArtifactProvider(args.netmhciipan_artifact)
        if args.netmhciipan_artifact
        else IEDBNetMHCIIpanProvider(ROOT / "data" / "processed" / "mhc_cache" / "iedb")
    )
    mhc_registry = MHCProviderRegistry([primary_mhc, *args.challenger])
    agent = ImmunogenicityScreeningAgent(
        response_registry=response_registry,
        mhc_registry=mhc_registry,
        default_response_adapter=args.default_response_adapter,
        hla_panel_path=ROOT / "docs" / "hla_class_ii_panel.v1.json",
        fusion_rule_path=ROOT / "docs" / "fusion_rule.v1.json",
        self_proteome_path=ROOT / "data" / "processed" / "self_proteome.parquet",
    )
    records = (
        read_fasta(args.fasta)
        if args.fasta
        else [(args.sequence_id, str(args.sequence).upper())]
    )
    accessibility = (
        json.loads(args.accessibility_artifact.read_text())["values_by_sha256"]
        if args.accessibility_artifact
        else {}
    )
    cleavage = (
        json.loads(args.cleavage_artifact.read_text())["values_by_sha256"]
        if args.cleavage_artifact
        else {}
    )
    shared_hla = (
        json.loads(args.shared_hla_artifact.read_text())["values_by_sha256"]
        if args.shared_hla_artifact
        else {}
    )
    assessments = []
    for candidate_id, sequence in records:
        key = sequence_sha256(sequence)
        assessments.append(
            agent.assess(
                candidate_id,
                sequence,
                accessibility_by_residue=accessibility.get(key),
                cleavage_by_window=cleavage.get(key),
                shared_hla_by_window=shared_hla.get(key),
            )
        )
    output = write_assessments(assessments, args.output)
    print(output)


if __name__ == "__main__":
    main()
