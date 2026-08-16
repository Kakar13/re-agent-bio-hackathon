"""Run artifacts and the pilot gates.

Everything a reader needs to disagree with the conclusion is written out: the
site table with its evidence columns, the peptides with their proteolytic
provenance, the per-arm geometry, and the gate verdicts. A gate that fails is
reported as failed; that is the point of declaring them before the run.
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

log = logging.getLogger(__name__)


@dataclass
class Gate:
    """One declared-in-advance pass/fail check."""

    gate_id: str
    question: str
    passed: bool | None  # None when the run could not evaluate it
    detail: str
    value: Any = None

    @property
    def verdict(self) -> str:
        if self.passed is None:
            return "NOT EVALUATED"
        return "PASS" if self.passed else "FAIL"


def write_csv(path: Path, rows: Sequence[dict], columns: list[str] | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return path
    # Union the keys so a row with extra fields cannot silently drop columns.
    cols = columns or list({k: None for r in rows for k in r})
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})
    return path


def _mannwhitney(a: list[float], b: list[float]) -> tuple[float | None, float | None]:
    """One-sided test that ``a`` is greater than ``b``."""
    if len(a) < 3 or len(b) < 3:
        return None, None
    try:
        from scipy.stats import mannwhitneyu

        stat, p = mannwhitneyu(a, b, alternative="greater")
        return float(stat), float(p)
    except Exception:  # noqa: BLE001
        return None, None


def _mean(xs: list[float]) -> float | None:
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def evaluate_gates(
    *,
    geometries: list[dict],
    peptide_summary: dict,
    positive_control: list[dict],
    clean_designs: dict,
) -> list[Gate]:
    """The five pilot gates, evaluated on whatever the run produced."""
    gates: list[Gate] = []

    def dists(arm: str) -> list[float]:
        return [
            g["nucleophile_distance"]
            for g in geometries
            if g.get("arm") == arm and g.get("nucleophile_distance") is not None
        ]

    # 1. Positive control -------------------------------------------------
    if positive_control:
        in_clip = [p for p in positive_control if p.get("in_clip_window")]
        engaged = [p for p in positive_control if p.get("engaged")]
        passed = bool(in_clip) and bool(engaged)
        gates.append(
            Gate(
                "cd74_recovery",
                "Does cathepsin S engage a CD74 site inside the CLIP window?",
                passed,
                f"{len(in_clip)}/{len(positive_control)} positive-control sites fall in the "
                f"CLIP window (97-125); {len(engaged)} are engaged in the co-fold.",
                {"in_clip": len(in_clip), "engaged": len(engaged)},
            )
        )
    else:
        gates.append(
            Gate("cd74_recovery", "Does cathepsin S engage a CD74 site inside the CLIP window?",
                 None, "No positive-control co-folds completed.")
        )

    # 2. Scramble separation ---------------------------------------------
    real, scram = dists("real"), dists("scramble")
    if real and scram:
        # Shorter distance is better, so test that scrambles sit further away.
        _, p = _mannwhitney(scram, real)
        mr, ms = _mean(real), _mean(scram)
        passed = bool(mr is not None and ms is not None and mr < ms)
        gates.append(
            Gate(
                "scramble_separation",
                "Do real sites sit closer to the nucleophile than composition-matched scrambles?",
                passed,
                f"mean nucleophile distance {mr:.2f} A real vs {ms:.2f} A scrambled"
                + (f", Mann-Whitney p={p:.3g}" if p is not None else ""),
                {"real": mr, "scramble": ms, "p_value": p},
            )
        )
    else:
        gates.append(Gate("scramble_separation", "Real versus scrambled separation", None,
                          f"insufficient data (real n={len(real)}, scramble n={len(scram)})"))

    # 3. Protease swap ----------------------------------------------------
    swap = dists("protease_swap")
    if real and swap:
        _, p = _mannwhitney(swap, real)
        mr, msw = _mean(real), _mean(swap)
        passed = bool(mr is not None and msw is not None and mr < msw)
        gates.append(
            Gate(
                "protease_swap_collapse",
                "Does the signal weaken when the segment is folded against the wrong protease?",
                passed,
                f"mean nucleophile distance {mr:.2f} A correct enzyme vs {msw:.2f} A swapped"
                + (f", Mann-Whitney p={p:.3g}" if p is not None else ""),
                {"real": mr, "swap": msw, "p_value": p},
            )
        )
    else:
        gates.append(Gate("protease_swap_collapse", "Protease-swap negative collapse", None,
                          f"insufficient data (real n={len(real)}, swap n={len(swap)})"))

    # 4. Clean designs stay clean ----------------------------------------
    if clean_designs.get("evaluated"):
        gates.append(
            Gate(
                "clean_designs_clean",
                "Do designs with no confident cleavage sites also yield few presented peptides?",
                clean_designs.get("passed"),
                clean_designs.get("detail", ""),
                clean_designs.get("value"),
            )
        )
    else:
        gates.append(Gate("clean_designs_clean",
                          "Do predicted-clean designs stay clean?", None,
                          clean_designs.get("detail", "not evaluated")))

    # 5. Constraint actually matters --------------------------------------
    c = peptide_summary.get("constrained_binders")
    u = peptide_summary.get("unconstrained_binders")
    if c is not None and u is not None and u > 0:
        frac = c / u
        passed = frac < 0.95  # the digest must remove a non-trivial share
        gates.append(
            Gate(
                "constraint_matters",
                "Does cleavage-first digestion change the epitope set versus a naive scan?",
                passed,
                f"{c} strong binders survive the digest against {u} from the unconstrained "
                f"scan ({100 * (1 - frac):.1f}% removed)",
                {"constrained": c, "unconstrained": u, "retained_fraction": frac},
            )
        )
    else:
        gates.append(Gate("constraint_matters", "Constrained versus unconstrained difference",
                          None, "peptide scoring did not complete for both arms"))

    return gates


def write_report(
    run_dir: Path,
    *,
    manifest: dict,
    sites: list[dict],
    peptides: list[dict],
    geometries: list[dict],
    gates: list[Gate],
    summary: dict,
) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)

    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    write_csv(run_dir / "sites.csv", sites)
    write_csv(run_dir / "peptides.csv", peptides)
    write_csv(run_dir / "concordance.csv", geometries)
    (run_dir / "gates.json").write_text(
        json.dumps([asdict(g) | {"verdict": g.verdict} for g in gates], indent=2, default=str)
    )

    md = _render_markdown(manifest, sites, peptides, geometries, gates, summary)
    path = run_dir / "report.md"
    path.write_text(md)
    log.info("wrote report to %s", path)
    return path


def _render_markdown(
    manifest: dict,
    sites: list[dict],
    peptides: list[dict],
    geometries: list[dict],
    gates: list[Gate],
    summary: dict,
) -> str:
    lines: list[str] = []
    a = lines.append

    a("# Cleavage-first protease attribution for de novo designs")
    a("")
    a(f"Run `{manifest.get('run_id')}` - {manifest.get('started_at')}")
    a("")
    a("Cleavage sites are resolved first, from MEROPS specificity matrices and the")
    a("design's own crystal structure. Only peptides that proteolysis could produce")
    a("are then scored for MHC presentation, and Boltz-2 co-folding is run afterwards")
    a("as an independent geometric check that never sees the matrix scores.")
    a("")

    a("## Pilot gates")
    a("")
    a("| Gate | Question | Verdict | Detail |")
    a("| --- | --- | --- | --- |")
    for g in gates:
        a(f"| `{g.gate_id}` | {g.question} | **{g.verdict}** | {g.detail} |")
    a("")

    a("## What was run")
    a("")
    for k, v in summary.items():
        a(f"- **{k.replace('_', ' ')}**: {v}")
    a("")

    if geometries:
        a("## Co-fold geometry by arm")
        a("")
        a("Distance is from the catalytic nucleophile to the carbonyl carbon of the")
        a("bond the sequence model predicted would be cut. Lower is better.")
        a("")
        a("| Arm | n | mean distance (A) | engaged | mean ipTM |")
        a("| --- | --- | --- | --- | --- |")
        arms: dict[str, list[dict]] = {}
        for g in geometries:
            arms.setdefault(g.get("arm", "?"), []).append(g)
        for arm, rows in sorted(arms.items()):
            d = _mean([r.get("nucleophile_distance") for r in rows])
            it = _mean([r.get("iptm") for r in rows])
            eng = sum(1 for r in rows if r.get("engaged"))
            a(
                f"| {arm} | {len(rows)} | "
                f"{'-' if d is None else f'{d:.2f}'} | {eng}/{len(rows)} | "
                f"{'-' if it is None else f'{it:.3f}'} |"
            )
        a("")

    if sites:
        a("## Top cleavage sites")
        a("")
        a("| Design | Bond | Protease | PWM score | p | Accessibility |")
        a("| --- | --- | --- | --- | --- | --- |")
        top = sorted(sites, key=lambda s: s.get("p_value", 1.0))[:20]
        for s in top:
            a(
                f"| {s.get('design_id')} | {s.get('scissile_bond')} | {s.get('protease')} | "
                f"{s.get('pwm_score', 0):.2f} | {s.get('p_value', 1):.4f} | "
                f"{s.get('accessibility')} |"
            )
        a("")

    a("## Limitations")
    a("")
    a("- MEROPS matrices are pooled across substrates and organisms, so they capture")
    a("  subsite preference rather than the kinetics of any one substrate.")
    a("- Accessibility is read from the crystal structure of the folded design. The")
    a("  endolysosome reduces disulfides and partially unfolds antigens, which is why")
    a("  buried sites are labelled `unfolding-required` rather than discarded.")
    a("- A Boltz-2 co-fold is a model, not a Michaelis complex. Distances are used")
    a("  comparatively against controls, never as absolute evidence of catalysis.")
    a("- Class II peptides are scored against a six-allele DRB1 panel; DP and DQ are")
    a("  not covered.")
    a("- Nothing here is experimental validation. The output is a ranked, falsifiable")
    a("  hypothesis set.")
    a("")
    return "\n".join(lines)


def new_run_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
