import type { MhcHit, RiskBreakdown, StructureFeatures, ToleranceHit } from "./types.ts";

export function scoreImmunoRisk(input: {
  sequenceId: string;
  features: StructureFeatures;
  mhc: MhcHit[];
  tolerance: ToleranceHit[];
}): RiskBreakdown {
  const factors: RiskBreakdown["factors"] = [];
  let score = 0;

  const strong = input.mhc.filter((h) => h.binderStub);
  const foreignStrong = strong.filter((h) =>
    input.tolerance.some((t) => t.peptide === h.peptide && t.status === "foreign_like"),
  );

  const mhcContrib = Math.min(40, strong.length * 8);
  factors.push({
    name: "mhc_stub_binders",
    contribution: mhcContrib,
    note: `${strong.length} stub MHC binders (rank≤2%); ${foreignStrong.length} also foreign_like`,
  });
  score += mhcContrib;

  const foreign = input.tolerance.filter((t) => t.status === "foreign_like").length;
  const foreignContrib = Math.min(30, foreign * 3);
  factors.push({
    name: "foreign_like_peptides",
    contribution: foreignContrib,
    note: `${foreign} peptides marked foreign_like vs stub self set`,
  });
  score += foreignContrib;

  // Low accessibility / high stability → clearance concern (team discussion)
  if (input.features.meanRsaProxy < 0.35 && input.features.unfoldingDgProxy > 6) {
    factors.push({
      name: "hyperstable_low_rsa",
      contribution: 15,
      note: "Low RSA proxy + high ΔG proxy — may resist protease clearance",
    });
    score += 15;
  }

  if (input.features.disorderFractionProxy > 0.5) {
    factors.push({
      name: "high_disorder",
      contribution: 10,
      note: "High disorder proxy — more cleavage accessibility (can raise peptide load)",
    });
    score += 10;
  }

  score = Math.min(100, score);
  const overall: RiskBreakdown["overall"] = score >= 55 ? "high" : score >= 30 ? "moderate" : "low";

  const peptidesFlagged = [
    ...new Set([
      ...foreignStrong.map((h) => h.peptide),
      ...input.tolerance.filter((t) => t.status === "foreign_like").map((t) => t.peptide),
    ]),
  ].slice(0, 40);

  return {
    sequenceId: input.sequenceId,
    overall,
    score0to100: score,
    factors,
    peptidesFlagged,
    method: "weighted_stub_v0",
    caveat: "Demo risk score only — recalibrate when real MHC + Atlas tolerance replace stubs.",
  };
}
