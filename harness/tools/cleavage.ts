import { getCatalyticSites } from "./catalytic-sites.ts";
import type { CatalyticSite, CleavageEvent } from "./types.ts";

function pushCut(
  out: CleavageEvent[],
  site: CatalyticSite,
  seq: string,
  p1Index: number,
): void {
  if (p1Index < 0 || p1Index >= seq.length - 1) return;
  const p1 = seq[p1Index]!;
  const p1Prime = seq[p1Index + 1]!;
  if (site.blockedP1Prime?.includes(p1Prime)) return;
  out.push({
    siteId: site.id,
    siteName: site.name,
    position: p1Index,
    p1,
    p1Prime,
    nTerminalProduct: seq.slice(0, p1Index + 1),
    cTerminalProduct: seq.slice(p1Index + 1),
  });
}

function matchFurin(seq: string, site: CatalyticSite, out: CleavageEvent[]): void {
  // R-X-[KR]-R↓  (P1 = final R)
  for (let i = 0; i <= seq.length - 4; i++) {
    if (seq[i] === "R" && (seq[i + 2] === "K" || seq[i + 2] === "R") && seq[i + 3] === "R") {
      pushCut(out, site, seq, i + 3);
    }
  }
}

function matchMmp(seq: string, site: CatalyticSite, out: CleavageEvent[]): void {
  for (let i = 0; i < seq.length - 1; i++) {
    if (seq[i] === "G" && seq[i + 1] === "P") pushCut(out, site, seq, i);
  }
}

/** Predict cleavage events against the curated catalytic-site catalog. */
export function predictCleavage(sequence: string, siteIds?: string[]): CleavageEvent[] {
  const sites = getCatalyticSites(siteIds);
  const out: CleavageEvent[] = [];

  for (const site of sites) {
    if (site.id === "furin_rxkr") {
      matchFurin(sequence, site, out);
      continue;
    }
    if (site.id === "mmp_gp") {
      matchMmp(sequence, site, out);
      continue;
    }
    if (!site.p1.length) continue;
    const p1set = new Set(site.p1);
    for (let i = 0; i < sequence.length - 1; i++) {
      if (p1set.has(sequence[i]!)) pushCut(out, site, sequence, i);
    }
  }

  out.sort((a, b) => a.position - b.position || a.siteId.localeCompare(b.siteId));
  return out;
}

/** Unique peptide fragments from cleavage (including full chain if no cuts). */
export function peptidePool(sequence: string, cleavages: CleavageEvent[], maxLen = 25): string[] {
  const cuts = [...new Set(cleavages.map((c) => c.position + 1))].sort((a, b) => a - b);
  const frags: string[] = [];
  let start = 0;
  for (const end of cuts) {
    if (end > start) frags.push(sequence.slice(start, end));
    start = end;
  }
  if (start < sequence.length) frags.push(sequence.slice(start));
  if (!frags.length) frags.push(sequence);

  // Always include the intact chain so dense cut maps still yield MHC-length windows.
  const sources = [...new Set([...frags, sequence])];

  const windowed = new Set<string>();
  for (const f of sources) {
    if (f.length >= 8 && f.length <= maxLen) windowed.add(f);
    for (let L = 8; L <= 11; L++) {
      for (let i = 0; i + L <= f.length; i++) windowed.add(f.slice(i, i + L));
    }
    for (let L = 12; L <= 18; L++) {
      for (let i = 0; i + L <= f.length; i++) windowed.add(f.slice(i, i + L));
    }
  }
  return [...windowed].sort((a, b) => a.length - b.length || a.localeCompare(b));
}
