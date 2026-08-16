import type { CatalyticSite } from "./types.ts";

/**
 * Catalytic-site catalog — keep in sync with
 * ``src/re_agent/immuno_risk/cleavage.py`` (Python is source of truth for cohort scans).
 * Motifs are simplified for weekend demo; cathepsin rules are P2-driven.
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
    pattern: "furin",
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
    pattern: "mmp",
    notes: "Pattern: GP soft sites via custom rule.",
  },
  {
    id: "legumain_n",
    name: "Legumain / AEP (N)",
    proteaseClass: "cysteine",
    motif: "P1 = N; prefer non-Pro P1'",
    p1: ["N"],
    blockedP1Prime: ["P"],
    notes: "Asn-specific endopeptidase (AEP/legumain) — MHC-II pathway relevant.",
  },
  {
    id: "cathepsin_s",
    name: "Cathepsin S (P2 hydrophobic)",
    proteaseClass: "cysteine",
    motif: "P2 = V/L/I/M/F; P1 broad (not P)",
    p1: ["A", "G", "S", "T", "N", "Q", "K", "R", "H", "L", "I", "V", "M", "F", "Y", "W", "E", "D"],
    p2: ["V", "L", "I", "M", "F"],
    blockedP1Prime: ["P"],
    notes:
      "Endosomal MHC-II processing. CatS is P2-driven (bulky hydrophobic); P1 is relatively permissive.",
  },
  {
    id: "cathepsin_l",
    name: "Cathepsin L (P2 aromatic/hydrophobic)",
    proteaseClass: "cysteine",
    motif: "P2 = F/Y/W/L; P1 broad",
    p1: ["A", "G", "S", "T", "N", "Q", "K", "R", "H", "L", "I", "V", "M", "F", "Y", "W", "E", "D"],
    p2: ["F", "Y", "W", "L"],
    blockedP1Prime: ["P"],
    notes: "Lysosomal endopeptidase; aromatic/Leu P2 preference.",
  },
  {
    id: "cathepsin_b",
    name: "Cathepsin B (Arg P1 + hydrophobic P2)",
    proteaseClass: "cysteine",
    motif: "P1 = R/K; P2 = hydrophobic",
    p1: ["R", "K"],
    p2: ["V", "L", "I", "M", "F", "A"],
    blockedP1Prime: ["P"],
    notes: "Endopeptidase mode: Arg/Lys P1 with hydrophobic P2.",
  },
  {
    id: "cathepsin_b_cpx",
    name: "Cathepsin B dipeptidyl-CPX",
    proteaseClass: "cysteine",
    motif: "Removes C-terminal dipeptides",
    p1: [],
    pattern: "ctsb_dipeptidyl",
    notes: "Dipeptidyl carboxypeptidase mode; cuts two residues from the C-terminus.",
  },
  {
    id: "immunoproteasome_b5i",
    name: "Immunoproteasome β5i / LMP7 (chymotrypsin-like)",
    proteaseClass: "threonine",
    motif: "P1 = F/Y/W/L (hydrophobic/aromatic)",
    p1: ["F", "Y", "W", "L"],
    notes: "IFN-γ-induced immunoproteasome subunit; MHC-I epitope C-termini.",
  },
  {
    id: "immunoproteasome_b1i",
    name: "Immunoproteasome β1i / LMP2 (caspase-like)",
    proteaseClass: "threonine",
    motif: "P1 = D/E (acidic)",
    p1: ["D", "E"],
    notes: "Immunoproteasome caspase-like activity.",
  },
  {
    id: "immunoproteasome_b2i",
    name: "Immunoproteasome β2i / MECL-1 (trypsin-like)",
    proteaseClass: "threonine",
    motif: "P1 = K/R (basic)",
    p1: ["K", "R"],
    blockedP1Prime: ["P"],
    notes: "Immunoproteasome trypsin-like activity.",
  },
];

export function getCatalyticSites(ids?: string[]): CatalyticSite[] {
  if (!ids?.length) return DEFAULT_CATALYTIC_SITES;
  const want = new Set(ids);
  return DEFAULT_CATALYTIC_SITES.filter((s) => want.has(s.id));
}
