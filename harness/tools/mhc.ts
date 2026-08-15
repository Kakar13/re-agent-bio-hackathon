import type { MhcHit } from "./types.ts";

const DEFAULT_MHC1 = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"];
const DEFAULT_MHC2 = ["HLA-DRB1*01:01", "HLA-DRB1*15:01"];

function stubRank(peptide: string, allele: string): number {
  // Deterministic pseudo-rank from char codes — NOT a real binder predictor.
  let h = 0;
  const s = peptide + "|" + allele;
  for (let i = 0; i < s.length; i++) h = (h * 33 + s.charCodeAt(i)) >>> 0;
  return ((h % 5000) / 100); // 0.00–49.99
}

/**
 * Stub MHC presentation scores. Swap for NetMHCpan / MHCflurry / MixMHC2pred later.
 */
export function scoreMhc(
  peptides: string[],
  opts?: { mhcClass?: "I" | "II" | "both"; alleles?: string[]; topN?: number },
): MhcHit[] {
  const which = opts?.mhcClass ?? "both";
  const topN = opts?.topN ?? 25;
  const hits: MhcHit[] = [];

  const run = (mhcClass: "I" | "II", alleles: string[], minL: number, maxL: number) => {
    for (const peptide of peptides) {
      if (peptide.length < minL || peptide.length > maxL) continue;
      for (const allele of alleles) {
        const rankPctStub = Math.round(stubRank(peptide, allele) * 100) / 100;
        hits.push({
          peptide,
          allele,
          mhcClass,
          length: peptide.length,
          rankPctStub,
          binderStub: rankPctStub <= 2.0,
          method: "stub_hash_rank_v0",
          caveat:
            "NOT NetMHCpan/MHCflurry — deterministic placeholder so the agent loop is wired. Replace before trusting binder calls.",
        });
      }
    }
  };

  if (which === "I" || which === "both") {
    run("I", opts?.alleles?.filter((a) => a.includes("HLA-A") || a.includes("HLA-B") || a.includes("HLA-C")) ?? DEFAULT_MHC1, 8, 11);
  }
  if (which === "II" || which === "both") {
    const a2 =
      opts?.alleles?.filter((a) => a.includes("DR") || a.includes("DQ") || a.includes("DP")) ?? DEFAULT_MHC2;
    run("II", a2, 12, 18);
  }

  hits.sort((a, b) => a.rankPctStub - b.rankPctStub);
  return hits.slice(0, topN);
}
