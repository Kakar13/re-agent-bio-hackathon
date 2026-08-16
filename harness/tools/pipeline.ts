import { mkdirSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { predictCleavage, peptidePool } from "./cleavage.ts";
import { estimateStructureFeatures } from "./features.ts";
import { scoreSequenceViaPython } from "./mhc.ts";
import { scoreImmunoRisk } from "./risk.ts";
import { parseFastaOrSequence } from "./text.ts";
import type {
  AggregationReport,
  ConfidenceReport,
  MhcHit,
  PipelineResult,
  ResidueRisk,
  RiskBreakdown,
  ToleranceHit,
} from "./types.ts";

function mapPythonResult(
  raw: Record<string, unknown>,
  features: ReturnType<typeof estimateStructureFeatures>,
  cleavages: ReturnType<typeof predictCleavage>,
  peptides: string[],
): PipelineResult {
  const riskRaw = raw.risk as Record<string, unknown>;
  const mhcI = riskRaw?.mhc_i as Record<string, unknown> | undefined;
  const mhcIi = riskRaw?.mhc_ii as Record<string, unknown> | undefined;

  const risk: RiskBreakdown = {
    sequenceId: String(raw.sequence_id),
    overall: (riskRaw?.overall as RiskBreakdown["overall"]) ?? "moderate",
    score0to100: Number(riskRaw?.score0to100 ?? 0),
    factors: (riskRaw?.factors as RiskBreakdown["factors"]) ?? [],
    peptidesFlagged: (riskRaw?.peptides_flagged as string[]) ?? [],
    method: String(riskRaw?.method ?? "python_dual_arm_v1"),
    caveat: String(riskRaw?.caveat ?? ""),
    mhcI: mhcI
      ? {
          binderCount: Number(mhcI.binder_count ?? 0),
          baselineScore: (mhcI.baseline_score as number | null) ?? null,
          iedbRiskScore: (mhcI.iedb_risk_score as number | null) ?? null,
        }
      : undefined,
    mhcII: mhcIi ? { binderCount: Number(mhcIi.binder_count ?? 0) } : undefined,
  };

  const peptidesEv = (raw.peptides as Array<Record<string, unknown>>) ?? [];
  const mhc: MhcHit[] = peptidesEv.map((e) => {
    const h = (e.mhc as Record<string, unknown>) ?? e;
    const binder = Boolean(h.binder);
    const pct = (h.percentile_rank as number | null) ?? null;
    return {
      peptide: String(e.peptide),
      allele: String(e.allele),
      mhcClass: e.mhc_class === "II" ? "II" : "I",
      length: Number(h.length ?? String(e.peptide).length),
      percentileRank: pct,
      presentationScore: (h.presentation_score as number | null) ?? null,
      affinityNm: (h.affinity_nm as number | null) ?? null,
      processingScore: (h.processing_score as number | null) ?? null,
      binder,
      method: String(h.method ?? ""),
      version: (h.version as string | null) ?? null,
      caveat: String(h.caveat ?? ""),
      binderStub: binder,
      rankPctStub: pct ?? undefined,
    };
  });

  const tolMap = new Map<string, ToleranceHit>();
  for (const e of peptidesEv) {
    const t = e.tolerance as Record<string, unknown> | null | undefined;
    if (!t) continue;
    tolMap.set(String(e.peptide), {
      peptide: String(e.peptide),
      nearestSelf: (t.nearest_self as string | null) ?? null,
      identity: Number(t.identity ?? 0),
      status: (t.status as ToleranceHit["status"]) ?? "unknown",
      atlasHit: Boolean(t.atlas_hit),
      method: String(t.method ?? ""),
      caveat: String(t.caveat ?? ""),
    });
  }

  const aggRaw = raw.aggregation as Record<string, unknown> | undefined;
  const aggregation: AggregationReport | undefined = aggRaw
    ? {
        sequenceId: String(aggRaw.sequence_id ?? raw.sequence_id),
        overall: aggRaw.overall as AggregationReport["overall"],
        score0to100: Number(aggRaw.score0to100 ?? 0),
        factors: (aggRaw.factors as AggregationReport["factors"]) ?? [],
        method: String(aggRaw.method ?? ""),
        caveat: String(aggRaw.caveat ?? ""),
      }
    : undefined;

  const confRaw = raw.confidence as Record<string, unknown> | undefined;
  const confidence: ConfidenceReport | undefined = confRaw
    ? {
        score0to1: Number(confRaw.score0to1 ?? 0),
        factors: (confRaw.factors as ConfidenceReport["factors"]) ?? [],
        method: String(confRaw.method ?? ""),
      }
    : undefined;

  return {
    runId: String(raw.run_id ?? ""),
    sequenceId: String(raw.sequence_id),
    sequence: String(raw.sequence),
    deliveryMode: String(raw.delivery_mode ?? "intracellular_plasmid"),
    features,
    cleavages,
    peptides: peptides.slice(0, 200),
    mhc,
    tolerance: [...tolMap.values()],
    risk,
    confidence,
    aggregation,
    residueRisk: (raw.residue_risk as ResidueRisk[]) ?? [],
    predictorVersions: (raw.predictor_versions as Record<string, string>) ?? {},
    caveats: (raw.caveats as string[]) ?? [],
    artifactDir: (raw.artifact_dir as string | null) ?? undefined,
  };
}

export function runImmunoPipeline(
  rawSequence: string,
  opts?: {
    sequenceId?: string;
    siteIds?: string[];
    mhcClass?: "I" | "II" | "both";
    writeDir?: string;
    repoRoot?: string;
  },
): PipelineResult {
  const parsed = parseFastaOrSequence(rawSequence);
  const sequenceId = opts?.sequenceId ?? parsed.id;
  const sequence = parsed.sequence;
  if (sequence.length < 8) {
    throw new Error(`sequence too short (${sequence.length}); need ≥8 AA`);
  }

  const features = estimateStructureFeatures(sequenceId, sequence);
  const cleavages = predictCleavage(sequence, opts?.siteIds);
  const peptides = peptidePool(sequence, cleavages);

  const py = scoreSequenceViaPython(sequence, {
    sequenceId,
    mhcClass: opts?.mhcClass ?? "both",
    write: Boolean(opts?.writeDir),
    repoRoot: opts?.repoRoot,
  });

  const result = mapPythonResult(py, features, cleavages, peptides);
  // Ensure risk object always present
  result.risk = scoreImmunoRisk({
    sequenceId,
    features,
    mhc: result.mhc,
    tolerance: result.tolerance,
    aggregation: result.aggregation,
    pythonRisk: result.risk,
  });

  if (opts?.writeDir && !result.artifactDir) {
    mkdirSync(opts.writeDir, { recursive: true });
    const path = join(opts.writeDir, `${sequenceId.replace(/[^\w.-]+/g, "_")}.json`);
    writeFileSync(path, JSON.stringify(result, null, 2));
    result.writtenPath = path;
  }

  return result;
}

export function defaultResultsDir(repoRoot: string): string {
  return join(repoRoot, "results", "immuno_risk");
}

export function ensureParent(filePath: string): void {
  mkdirSync(dirname(filePath), { recursive: true });
}
