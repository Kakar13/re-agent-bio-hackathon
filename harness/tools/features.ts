import type { StructureFeatures } from "./types.ts";

const HYDROPHOBIC = new Set("AILMFVWY".split(""));
const DISORDER_PRONE = new Set("GSTNPQEDKR".split(""));

/**
 * Sequence-only proxies for RSA / disorder / ΔG until real structure or MPNN scorefiles land.
 * Labelled clearly so judges do not confuse them with DSSP / Rosetta numbers.
 */
export function estimateStructureFeatures(sequenceId: string, sequence: string): StructureFeatures {
  const n = sequence.length || 1;
  let hydro = 0;
  let disorder = 0;
  for (const aa of sequence) {
    if (HYDROPHOBIC.has(aa)) hydro++;
    if (DISORDER_PRONE.has(aa)) disorder++;
  }
  const hydroFrac = hydro / n;
  const disorderFractionProxy = disorder / n;
  // Crude RSA: ends more exposed; hydrophobic core buried
  const meanRsaProxy = Math.min(1, Math.max(0, 0.55 - 0.35 * hydroFrac + 0.15 * disorderFractionProxy));
  // Long "loops": runs of disorder-prone ≥ 8
  let longLoopCountProxy = 0;
  let run = 0;
  for (const aa of sequence) {
    if (DISORDER_PRONE.has(aa)) {
      run++;
      if (run === 8) longLoopCountProxy++;
    } else run = 0;
  }
  // Unfolding ΔG proxy: higher hydro + length → more stable (arbitrary kcal scale)
  const unfoldingDgProxy = Math.round((2.5 + 8 * hydroFrac + 0.01 * n - 4 * disorderFractionProxy) * 10) / 10;
  const secondaryGuess =
    hydroFrac > 0.4 && disorderFractionProxy < 0.35
      ? "helix/sheet-leaning (hydrophobic)"
      : disorderFractionProxy > 0.45
        ? "disorder-leaning"
        : "mixed";

  return {
    sequenceId,
    length: sequence.length,
    meanRsaProxy: Math.round(meanRsaProxy * 1000) / 1000,
    disorderFractionProxy: Math.round(disorderFractionProxy * 1000) / 1000,
    longLoopCountProxy,
    unfoldingDgProxy,
    secondaryGuess,
    method: "sequence_heuristic_v0",
    caveat:
      "PLACEHOLDER — replace with structure RSA / MPNN scorefile / Proto metrics when upstream artifacts exist.",
  };
}
