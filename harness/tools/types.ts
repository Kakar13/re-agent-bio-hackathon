/** Shared types for immuno-risk Pi tools (evidence-backed dual-arm pipeline). */

export interface CatalyticSite {
  id: string;
  name: string;
  proteaseClass: string;
  motif: string;
  p1: string[];
  /** Preferred residues at P2 (N-terminal of P1). Empty = unrestricted. */
  p2?: string[];
  /** Preferred residues at P3. Empty = unrestricted. */
  p3?: string[];
  /** Preferred residues at P1'. Empty = unrestricted. */
  p1Prime?: string[];
  blockedP1Prime?: string[];
  notes?: string;
  /** Special pattern: furin | mmp | ctsb_dipeptidyl */
  pattern?: string;
}

export interface CleavageEvent {
  siteId: string;
  siteName: string;
  proteaseClass?: string;
  position: number;
  p1: string;
  p1Prime: string;
  p2?: string;
  p3?: string;
  score?: number;
  nTerminalProduct: string;
  cTerminalProduct: string;
}

export interface StructureFeatures {
  sequenceId: string;
  length: number;
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
  start?: number | null;
  end?: number | null;
  affinityNm?: number | null;
  presentationScore?: number | null;
  processingScore?: number | null;
  percentileRank?: number | null;
  binder: boolean;
  method: string;
  version?: string | null;
  caveat: string;
  /** @deprecated stub field — prefer percentileRank */
  rankPctStub?: number;
  binderStub?: boolean;
}

export interface ToleranceHit {
  peptide: string;
  nearestSelf: string | null;
  identity: number;
  status: "self_like" | "foreign_like" | "unknown";
  atlasHit?: boolean;
  method: string;
  caveat: string;
}

export interface AggregationReport {
  sequenceId: string;
  overall: "low" | "moderate" | "high";
  score0to100: number;
  factors: Array<{ name: string; contribution: number; note: string }>;
  method: string;
  caveat: string;
}

export interface ConfidenceReport {
  score0to1: number;
  factors: Array<{ name: string; contribution: number }>;
  method: string;
}

export interface ResidueRisk {
  position: number;
  residue: string;
  risk: number;
  peptideCount: number;
  peptides: string[];
}

export interface RiskBreakdown {
  sequenceId: string;
  overall: "low" | "moderate" | "high";
  score0to100: number;
  factors: Array<{ name: string; contribution: number; note: string }>;
  peptidesFlagged: string[];
  method: string;
  caveat: string;
  mhcI?: { binderCount: number; baselineScore?: number | null; iedbRiskScore?: number | null };
  mhcII?: { binderCount: number };
}

export interface PipelineResult {
  runId?: string;
  sequenceId: string;
  sequence: string;
  deliveryMode?: string;
  features?: StructureFeatures;
  cleavages?: CleavageEvent[];
  peptides: string[] | unknown[];
  mhc: MhcHit[];
  tolerance: ToleranceHit[];
  risk: RiskBreakdown;
  confidence?: ConfidenceReport;
  aggregation?: AggregationReport;
  residueRisk?: ResidueRisk[];
  predictorVersions?: Record<string, string>;
  caveats?: string[];
  artifactDir?: string;
  writtenPath?: string;
}
