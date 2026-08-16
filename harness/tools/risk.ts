import type {
  AggregationReport,
  ConfidenceReport,
  MhcHit,
  RiskBreakdown,
  StructureFeatures,
  ToleranceHit,
} from "./types.ts";

/**
 * Prefer Python-computed risk from the bridge. This fallback only combines
 * already-evidence-backed MHC + tolerance hits (no stub_hash ranks).
 */
export function scoreImmunoRisk(input: {
  sequenceId: string;
  features?: StructureFeatures;
  mhc: MhcHit[];
  tolerance: ToleranceHit[];
  aggregation?: AggregationReport;
  confidence?: ConfidenceReport;
  pythonRisk?: RiskBreakdown;
}): RiskBreakdown {
  if (input.pythonRisk) return input.pythonRisk;

  const factors: RiskBreakdown["factors"] = [];
  let score = 0;

  const strong = input.mhc.filter((h) => h.binder || h.binderStub);
  const foreignStrong = strong.filter((h) =>
    input.tolerance.some((t) => t.peptide === h.peptide && t.status === "foreign_like"),
  );

  const mhcContrib = Math.min(45, strong.filter((h) => h.mhcClass === "I").length * 6);
  factors.push({
    name: "mhc_i_binders",
    contribution: mhcContrib,
    note: `${strong.filter((h) => h.mhcClass === "I").length} MHC-I binders; ${foreignStrong.length} foreign_like`,
  });
  score += mhcContrib;

  const foreign = input.tolerance.filter((t) => t.status === "foreign_like").length;
  const foreignContrib = Math.min(25, foreign * 4);
  factors.push({
    name: "foreign_like_peptides",
    contribution: foreignContrib,
    note: `${foreign} foreign_like vs Atlas reference`,
  });
  score += foreignContrib;

  const ii = strong.filter((h) => h.mhcClass === "II").length;
  if (ii) {
    const c = Math.min(15, ii * 2);
    factors.push({
      name: "mhc_ii_presentation_thin",
      contribution: c,
      note: `${ii} MHC-II presentation hits (not ADA)`,
    });
    score += c;
  }

  if (input.aggregation && input.aggregation.score0to100 >= 55) {
    factors.push({
      name: "high_aggregation_context",
      contribution: 5,
      note: "High aggregation report — contextual only",
    });
    score += 5;
  }

  score = Math.min(100, score);
  const overall: RiskBreakdown["overall"] = score >= 55 ? "high" : score >= 30 ? "moderate" : "low";

  return {
    sequenceId: input.sequenceId,
    overall,
    score0to100: score,
    factors,
    peptidesFlagged: [
      ...new Set([
        ...foreignStrong.map((h) => h.peptide),
        ...input.tolerance.filter((t) => t.status === "foreign_like").map((t) => t.peptide),
      ]),
    ].slice(0, 40),
    method: "ts_evidence_combine_v1",
    caveat:
      "Uses Python MHC/Atlas evidence when available. Screening score — not clinical probability.",
  };
}
