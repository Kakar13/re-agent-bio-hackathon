/**
 * Tolerance check via Python Atlas join (no fake viral self set).
 */
import { scoreSequenceViaPython } from "./mhc.ts";
import type { ToleranceHit } from "./types.ts";

export function checkTolerance(
  peptides: string[],
  opts?: { sequence?: string },
): ToleranceHit[] {
  if (!peptides.length) return [];
  const seq = opts?.sequence ?? peptides.join("");
  if (seq.length < 8) {
    return peptides.map((peptide) => ({
      peptide,
      nearestSelf: null,
      identity: 0,
      status: "unknown" as const,
      method: "insufficient_sequence",
      caveat: "Need source sequence context for Atlas join via Python backend.",
    }));
  }
  try {
    const result = scoreSequenceViaPython(seq, { write: false });
    const ev = (result.peptides as Array<Record<string, unknown>>) ?? [];
    const byPep = new Map<string, ToleranceHit>();
    for (const e of ev) {
      const tol = e.tolerance as Record<string, unknown> | null | undefined;
      if (!tol) continue;
      const peptide = String(e.peptide ?? tol.peptide);
      byPep.set(peptide, {
        peptide,
        nearestSelf: (tol.nearest_self as string | null) ?? null,
        identity: Number(tol.identity ?? 0),
        status: (tol.status as ToleranceHit["status"]) ?? "unknown",
        atlasHit: Boolean(tol.atlas_hit),
        method: String(tol.method ?? "atlas"),
        caveat: String(tol.caveat ?? ""),
      });
    }
    return peptides.map(
      (p) =>
        byPep.get(p) ??
        byPep.get(p.toUpperCase()) ?? {
          peptide: p,
          nearestSelf: null,
          identity: 0,
          status: "unknown" as const,
          method: "not_in_top_hits",
          caveat: "Peptide not among scored top hits for Atlas join.",
        },
    );
  } catch (e) {
    throw new Error(
      `check_tolerance requires Python immuno backend: ${e instanceof Error ? e.message : String(e)}`,
    );
  }
}
