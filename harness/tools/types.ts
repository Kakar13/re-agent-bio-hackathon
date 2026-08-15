/** Shared types for immuno-risk Pi tools (late-stage de novo pipeline). */

export interface CatalyticSite {
  id: string;
  name: string;
  proteaseClass: string;
  /** One-letter motif description for humans / judges. */
  motif: string;
  /** Residues after which cleavage is preferred (P1). Empty = custom rule. */
  p1: string[];
  /** Optional residues blocked at P1'. */
  blockedP1Prime?: string[];
  notes?: string;
}

export interface CleavageEvent {
  siteId: string;
  siteName: string;
  position: number; // 0-based index of P1 residue
  p1: string;
  p1Prime: string;
  nTerminalProduct: string;
  cTerminalProduct: string;
}

export interface StructureFeatures {
  sequenceId: string;
  length: number;
  /** Mean relative solvent accessibility proxy in [0,1] (heuristic until real RSA). */
  meanRsaProxy: number;
  disorderFractionProxy: number;
  longLoopCountProxy: number;
  unfoldingDgProxy: number;
  secondaryGuess: string;
  method: string;
  caveat: string;
}

export interface MhcHit {
  peptide: string;
  allele: string;
  mhcClass: "I" | "II";
  length: number;
  /** Stub rank percentile (lower = stronger binder). Not NetMHCpan. */
  rankPctStub: number;
  binderStub: boolean;
  method: string;
  caveat: string;
}

export interface ToleranceHit {
  peptide: string;
  nearestSelf: string | null;
  identity: number;
  status: "self_like" | "foreign_like" | "unknown";
  method: string;
  caveat: string;
}

export interface RiskBreakdown {
  sequenceId: string;
  overall: "low" | "moderate" | "high";
  score0to100: number;
  factors: Array<{ name: string; contribution: number; note: string }>;
  peptidesFlagged: string[];
  method: string;
  caveat: string;
}

export interface PipelineResult {
  sequenceId: string;
  sequence: string;
  features: StructureFeatures;
  cleavages: CleavageEvent[];
  peptides: string[];
  mhc: MhcHit[];
  tolerance: ToleranceHit[];
  risk: RiskBreakdown;
  writtenPath?: string;
}
