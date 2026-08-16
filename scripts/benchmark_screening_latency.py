#!/usr/bin/env python3
"""Measure the one advantage a distilled student has over its own teacher.

The student cannot beat NetMHCpan on accuracy; it is trained to imitate it. The
argument for it is that scanning every 9-mer of a design is affordable enough to
run inside a generation loop. That claim was asserted in the writeup and never
measured, so this measures it.

The comparison is deliberately end-to-end and unflattering to us: the student
timing includes the ESM-2 forward pass, which dominates, not just the small head.
Both paths answer the same question - score every 9-mer window in this protein
against HLA-A*02:01 - so the numbers are directly comparable as an operational
choice, even though one is local compute and the other is a network round trip.
"""

from __future__ import annotations

import argparse
import json
import tempfile
import time
from pathlib import Path

import pandas as pd

from re_agent.e2e_pls import data
from re_agent.immuno.e2e_pls_pickle import TeamE2EPLSAdapter, _ESM2SegmentEncoder
from re_agent.immuno.netmhcpan import NetMHCpanTeacher

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CHECKPOINT = REPO_ROOT / "models/chao1/cv5_heads.pkl 2"
DEFAULT_STUDENT = REPO_ROOT / "models/a0201-netmhcpan-pda-cv5-v4/checkpoint"
DEFAULT_OUTPUT = REPO_ROOT / "results/benchmarks/screening_latency"
# A representative single-domain design: long enough that per-call overhead does
# not dominate, short enough to stay inside one IEDB request.
PROTEIN = (
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"
    "QTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG"
)
WARMUP_PROTEIN = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEFELRRQACGRTHVDLNQHMKSAWDEVL"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--student-dir", type=Path, default=DEFAULT_STUDENT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--skip-teacher", action="store_true")
    args = parser.parse_args()

    windows = data.tile_protein(
        parent_id="latency-probe", sequence=PROTEIN, source_domain="de_novo"
    )
    peptides = [str(window["peptide"]) for window in windows]
    print(f"scoring {len(peptides)} 9-mer windows from a {len(PROTEIN)}-residue protein\n")

    adapter = TeamE2EPLSAdapter(args.checkpoint, netmhcpan_checkpoint_dir=args.student_dir)
    # One untimed pass on a decoy so lazy weight loading does not land in the
    # measurement. It must not be PROTEIN: the encoder memoizes its last input,
    # so warming up on the timed sequence would measure a cache hit instead.
    adapter.predict(WARMUP_PROTEIN)

    student_times = []
    for index in range(args.repeats):
        _ESM2SegmentEncoder._last_sequences = None
        _ESM2SegmentEncoder._last_embeddings = None
        start = time.perf_counter()
        result = adapter.predict(PROTEIN)
        elapsed = time.perf_counter() - start
        student_times.append(elapsed)
        print(f"  student run {index + 1}: {elapsed:.2f}s ({len(result.predictions)} windows)")
    student_best = min(student_times)

    report = {
        "n_windows": len(peptides),
        "protein_length": len(PROTEIN),
        "student": {
            "runs_seconds": student_times,
            "best_seconds": student_best,
            "peptides_per_second": len(peptides) / student_best,
            "includes": "ESM-2 embedding, chao1 cleavage and TAP heads, student ensemble",
            "device": "cpu",
        },
    }

    if not args.skip_teacher:
        # A cold cache directory, so this times a real request rather than a
        # cache hit from an earlier corpus build.
        with tempfile.TemporaryDirectory() as cache_dir:
            teacher = NetMHCpanTeacher(Path(cache_dir))
            start = time.perf_counter()
            teacher.label(pd.DataFrame({"peptide": peptides}))
            teacher_seconds = time.perf_counter() - start
        report["teacher"] = {
            "seconds": teacher_seconds,
            "peptides_per_second": len(peptides) / teacher_seconds,
            "transport": "IEDB Tools API, single batched request",
        }
        report["speedup"] = teacher_seconds / student_best

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# Scanning cost: student vs NetMHCpan teacher",
        "",
        f"Every one of the {report['n_windows']} 9-mer windows in a "
        f"{report['protein_length']}-residue design, scored against HLA-A*02:01.",
        "",
        "| Path | Wall clock | Windows/second |",
        "| --- | ---: | ---: |",
        f"| Student (CPU, end to end) | {student_best:.2f}s "
        f"| {report['student']['peptides_per_second']:.0f} |",
    ]
    if "teacher" in report:
        teacher = report["teacher"]
        lines.append(
            f"| NetMHCpan via IEDB API | {teacher['seconds']:.2f}s "
            f"| {teacher['peptides_per_second']:.0f} |"
        )
        lines += [
            "",
            f"**{report['speedup']:.1f}x faster end to end.** The student figure includes the "
            "ESM-2 forward pass, which dominates it; the head itself is negligible.",
        ]
    lines += [
        "",
        "## Boundaries",
        "",
        "- This compares local compute against a network round trip, which is the real "
        "operational choice but is not a like-for-like model benchmark. A local "
        "NetMHCpan binary would close most of the gap.",
        "- Measured on CPU. A GPU widens the student's margin, and the embedding cost "
        "amortizes further when scoring many designs in one batch.",
        "- The teacher is rate-limited in a way a single timed call does not capture, so "
        "this understates the difference over a full campaign.",
        "",
    ]
    (args.output_dir / "REPORT.md").write_text("\n".join(lines))
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
