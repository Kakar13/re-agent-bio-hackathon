import type { ToleranceHit } from "./types.ts";

/** Tiny self reference until HLA Ligand Atlas / proteome peptides are loaded. */
export const SELF_PEPTIDE_STUBS = [
  "GILGFVFTL", // flu matrix A2 epitope — for demo contrast, treat as "known ligand-like"
  "NLVPMVATV",
  "GLCTLVAML",
  "STAPPVHNV",
  "LLFGYPVYV",
  "KVFGSLAFV",
  "ILKEPVHGV",
  "SLLMWITQV",
];

function identity(a: string, b: string): number {
  if (!a.length || !b.length) return 0;
  const n = Math.min(a.length, b.length);
  let m = 0;
  for (let i = 0; i < n; i++) if (a[i] === b[i]) m++;
  return m / Math.max(a.length, b.length);
}

/**
 * Tolerance / self-likeness check. Real path: HLA Ligand Atlas healthy refs.
 */
export function checkTolerance(
  peptides: string[],
  selfRef: string[] = SELF_PEPTIDE_STUBS,
): ToleranceHit[] {
  return peptides.map((peptide) => {
    let best: string | null = null;
    let bestId = 0;
    for (const s of selfRef) {
      const id = identity(peptide, s);
      if (id > bestId) {
        bestId = id;
        best = s;
      }
    }
    const status: ToleranceHit["status"] =
      bestId >= 0.85 ? "self_like" : bestId <= 0.35 ? "foreign_like" : "unknown";
    return {
      peptide,
      nearestSelf: best,
      identity: Math.round(bestId * 1000) / 1000,
      status,
      method: "identity_vs_stub_self_set_v0",
      caveat:
        "PLACEHOLDER self set — wire HLA Ligand Atlas / healthy immunopeptidome before claiming clinical tolerance.",
    };
  });
}
