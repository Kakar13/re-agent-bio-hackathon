import type { CatalyticSite } from "./types.ts";

/**
 * Starter catalog (~10 characterized sites). Motifs are simplified for a weekend demo;
 * replace with Mark's curated set + literature citations when ready.
 */
export const DEFAULT_CATALYTIC_SITES: CatalyticSite[] = [
  {
    id: "trypsin_kr",
    name: "Trypsin-like (K/R)",
    proteaseClass: "serine",
    motif: "P1 = K/R; blocked if P1' = P",
    p1: ["K", "R"],
    blockedP1Prime: ["P"],
    notes: "Classic tryptic cut; intracellular processing proxy.",
  },
  {
    id: "chymotrypsin_fwy",
    name: "Chymotrypsin-like (F/Y/W)",
    proteaseClass: "serine",
    motif: "P1 = F/Y/W",
    p1: ["F", "Y", "W"],
    notes: "Aromatic P1 preference.",
  },
  {
    id: "caspase_asp",
    name: "Caspase-like (D)",
    proteaseClass: "cysteine",
    motif: "P1 = D (simplified; real caspases need DXXD-like context)",
    p1: ["D"],
    notes: "Stub — tighten to DXXD when evaluating apoptosis-linked paths.",
  },
  {
    id: "furin_rxkr",
    name: "Furin-like (R-X-K/R-R)",
    proteaseClass: "serine",
    motif: "R-X-[KR]-R↓",
    p1: [],
    notes: "Handled by dedicated pattern matcher, not plain P1 list.",
  },
  {
    id: "thrombin_r",
    name: "Thrombin-like (R)",
    proteaseClass: "serine",
    motif: "P1 = R (simplified)",
    p1: ["R"],
    notes: "Overlaps trypsin; kept as separate labeled site for demos.",
  },
  {
    id: "pepsin_fl",
    name: "Pepsin-like (F/L)",
    proteaseClass: "aspartic",
    motif: "P1 = F/L (simplified)",
    p1: ["F", "L"],
  },
  {
    id: "elastase_agv",
    name: "Elastase-like (A/G/V)",
    proteaseClass: "serine",
    motif: "P1 = A/G/V",
    p1: ["A", "G", "V"],
  },
  {
    id: "proteasome_hydrophobic",
    name: "Proteasome-like (hydrophobic C-term preference)",
    proteaseClass: "threonine",
    motif: "P1 = L/I/V/F/Y (cytosolic MHC I antigen processing proxy)",
    p1: ["L", "I", "V", "F", "Y"],
    notes: "Coarse stand-in for immunoproteasome cuts feeding MHC I.",
  },
  {
    id: "mmp_gp",
    name: "MMP-like (G-P soft site)",
    proteaseClass: "metallo",
    motif: "P1–P1' ≈ G-P / P-X (simplified GP motif scan)",
    p1: [],
    notes: "Pattern: GP or PX soft sites via custom rule.",
  },
  {
    id: "legumain_n",
    name: "Legumain-like (N)",
    proteaseClass: "cysteine",
    motif: "P1 = N",
    p1: ["N"],
    notes: "Asn-specific endopeptidase stub.",
  },
];

export function getCatalyticSites(ids?: string[]): CatalyticSite[] {
  if (!ids?.length) return DEFAULT_CATALYTIC_SITES;
  const want = new Set(ids);
  return DEFAULT_CATALYTIC_SITES.filter((s) => want.has(s.id));
}
