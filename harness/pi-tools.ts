/**
 * Immuno-risk custom tools for Pi — evidence-backed MHC-I + thin MHC-II + aggregation + Benchling.
 */
import { Type, type Static } from "typebox";
import { resolve } from "node:path";
import { existsSync } from "node:fs";
import { spawnSync } from "node:child_process";
import { getCatalyticSites } from "./tools/catalytic-sites.ts";
import { peptidePool, predictCleavage } from "./tools/cleavage.ts";
import { estimateStructureFeatures } from "./tools/features.ts";
import { scoreMhc } from "./tools/mhc.ts";
import { defaultResultsDir, runImmunoPipeline } from "./tools/pipeline.ts";
import { scoreImmunoRisk } from "./tools/risk.ts";
import { checkTolerance } from "./tools/tolerance.ts";
import { jsonResult, parseFastaOrSequence, text } from "./tools/text.ts";
import type { MhcHit, StructureFeatures, ToleranceHit } from "./tools/types.ts";

export interface ImmunoToolsOptions {
  repoRoot?: string;
}

export interface ImmunoToolDefinition {
  name: string;
  label: string;
  description: string;
  parameters: unknown;
  execute: (
    toolCallId: string,
    params: never,
    signal?: AbortSignal,
  ) => Promise<{ content: Array<{ type: string; text?: string }>; details?: unknown }>;
}

function repoRootFromCwd(): string {
  return resolve(process.cwd(), "..");
}

function runPythonCli(root: string, args: string[]): Record<string, unknown> {
  const env = {
    ...process.env,
    PYTHONPATH: `${root}/src${process.env.PYTHONPATH ? `:${process.env.PYTHONPATH}` : ""}`,
    IMMUNO_ALLOW_HEURISTIC_MHC: process.env.IMMUNO_ALLOW_HEURISTIC_MHC ?? "1",
  };
  const candidates = [
    resolve(root, ".venv/bin/python"),
    "/tmp/immuno-venv/bin/python",
    process.env.IMMUNO_PYTHON ?? "",
    "python3",
  ].filter(Boolean);
  let proc: ReturnType<typeof spawnSync> | null = null;
  for (const cmd of candidates) {
    if (cmd !== "python3" && !existsSync(cmd)) continue;
    proc = spawnSync(cmd, ["-m", "re_agent.immuno_risk.cli", ...args], {
      cwd: root,
      encoding: "utf-8",
      env,
      maxBuffer: 32 * 1024 * 1024,
    });
    if (proc.status === 0) break;
  }
  if (!proc || proc.status !== 0) {
    throw new Error((proc?.stderr || proc?.stdout || "python cli failed").slice(0, 800));
  }
  try {
    return JSON.parse(proc.stdout) as Record<string, unknown>;
  } catch {
    return { raw: proc.stdout };
  }
}

export function createImmunoRiskTools(opts?: ImmunoToolsOptions): ImmunoToolDefinition[] {
  const root = opts?.repoRoot ?? repoRootFromCwd();

  const listSitesParams = Type.Object({});
  const structureParams = Type.Object({
    sequence: Type.String({ description: "FASTA or raw amino-acid sequence" }),
    sequence_id: Type.Optional(Type.String()),
  });
  const cleavageParams = Type.Object({
    sequence: Type.String({ description: "FASTA or raw amino-acid sequence" }),
    site_ids: Type.Optional(Type.Array(Type.String())),
  });
  const mhcParams = Type.Object({
    peptides: Type.Optional(Type.Array(Type.String())),
    sequence: Type.Optional(Type.String({ description: "Preferred: full protein sequence" })),
    mhc_class: Type.Optional(
      Type.Union([Type.Literal("I"), Type.Literal("II"), Type.Literal("both")]),
    ),
    top_n: Type.Optional(Type.Number()),
  });
  const toleranceParams = Type.Object({
    peptides: Type.Array(Type.String()),
    sequence: Type.Optional(Type.String()),
  });
  const riskParams = Type.Object({
    sequence_id: Type.String(),
    features: Type.Optional(Type.Record(Type.String(), Type.Unknown())),
    mhc: Type.Array(Type.Record(Type.String(), Type.Unknown())),
    tolerance: Type.Array(Type.Record(Type.String(), Type.Unknown())),
  });
  const pipelineParams = Type.Object({
    sequence: Type.String({ description: "FASTA or raw amino-acid sequence" }),
    sequence_id: Type.Optional(Type.String()),
    mhc_class: Type.Optional(
      Type.Union([Type.Literal("I"), Type.Literal("II"), Type.Literal("both")]),
    ),
    site_ids: Type.Optional(Type.Array(Type.String())),
    write: Type.Optional(Type.Boolean()),
  });
  const benchPullParams = Type.Object({
    ids: Type.Optional(Type.String({ description: "Comma-separated Benchling AA sequence IDs" })),
    name_includes: Type.Optional(Type.String()),
    dry_run: Type.Optional(Type.Boolean()),
  });
  const benchPubParams = Type.Object({
    run_dir: Type.String({ description: "Path to results/immuno_risk/<run-id>/" }),
    dry_run: Type.Optional(Type.Boolean()),
  });

  const tools: ImmunoToolDefinition[] = [
    {
      name: "list_catalytic_sites",
      label: "List catalytic sites",
      description:
        "List curated protease / catalytic sites used as a diagnostic cleavage layer (~10 starters).",
      parameters: listSitesParams,
      async execute() {
        const sites = getCatalyticSites();
        return jsonResult({ sites, count: sites.length }, "Curated catalytic sites:");
      },
    },
    {
      name: "structure_features",
      label: "Structure features",
      description:
        "Estimate RSA / disorder / loop / unfolding-ΔG proxies (sequence heuristics until real structure provided).",
      parameters: structureParams,
      async execute(_id, params: Static<typeof structureParams>) {
        const parsed = parseFastaOrSequence(params.sequence);
        const id = params.sequence_id ?? parsed.id;
        if (parsed.sequence.length < 8) return text("error: sequence too short (need ≥8 AA)");
        return jsonResult(estimateStructureFeatures(id, parsed.sequence));
      },
    },
    {
      name: "predict_cleavage",
      label: "Predict cleavage",
      description:
        "Diagnostic cleavage vs catalytic sites + peptide pool. MHC-I funnel is driven by MHCflurry processing/presentation, not this motif count.",
      parameters: cleavageParams,
      async execute(_id, params: Static<typeof cleavageParams>) {
        const parsed = parseFastaOrSequence(params.sequence);
        if (parsed.sequence.length < 8) return text("error: sequence too short (need ≥8 AA)");
        const cleavages = predictCleavage(parsed.sequence, params.site_ids);
        const peptides = peptidePool(parsed.sequence, cleavages);
        return jsonResult({
          sequenceId: parsed.id,
          length: parsed.sequence.length,
          cleavages,
          peptideCount: peptides.length,
          peptides: peptides.slice(0, 150),
        });
      },
    },
    {
      name: "score_mhc",
      label: "Score MHC",
      description:
        "Score MHC-I (MHCflurry ± optional NetMHCpan) and thin MHC-II (HLAIIPred heuristic fallback). Prefer passing full sequence.",
      parameters: mhcParams,
      async execute(_id, params: Static<typeof mhcParams>) {
        try {
          const hits = scoreMhc(params.peptides ?? [], {
            mhcClass: params.mhc_class ?? "both",
            topN: params.top_n ?? 40,
            sequence: params.sequence,
          });
          return jsonResult({
            hits,
            binderCount: hits.filter((h) => h.binder).length,
          });
        } catch (e) {
          return text(`error: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
    {
      name: "check_tolerance",
      label: "Check tolerance",
      description:
        "Compare peptides to HLA Ligand Atlas-style benign ligands (exact + identity proxy). Not clinical tolerance.",
      parameters: toleranceParams,
      async execute(_id, params: Static<typeof toleranceParams>) {
        if (!params.peptides?.length) return text("error: peptides array required");
        try {
          return jsonResult({
            hits: checkTolerance(params.peptides, { sequence: params.sequence }),
          });
        } catch (e) {
          return text(`error: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
    {
      name: "score_immuno_risk",
      label: "Score immuno risk",
      description:
        "Combine MHC + Atlas tolerance into screening risk (prefer run_immuno_pipeline for full dual-arm + aggregation).",
      parameters: riskParams,
      async execute(_id, params: Static<typeof riskParams>) {
        const risk = scoreImmunoRisk({
          sequenceId: params.sequence_id,
          features: params.features as unknown as StructureFeatures | undefined,
          mhc: params.mhc as unknown as MhcHit[],
          tolerance: params.tolerance as unknown as ToleranceHit[],
        });
        return jsonResult(risk);
      },
    },
    {
      name: "run_immuno_pipeline",
      label: "Run immuno pipeline",
      description:
        "End-to-end: MHC-I (MHCflurry + IEDB head) + thin MHC-II + Atlas + aggregation + residue map. Writes results/immuno_risk/<run-id>/ when write=true.",
      parameters: pipelineParams,
      async execute(_id, params: Static<typeof pipelineParams>) {
        try {
          const write = params.write !== false;
          const result = runImmunoPipeline(params.sequence, {
            sequenceId: params.sequence_id,
            siteIds: params.site_ids,
            mhcClass: params.mhc_class ?? "both",
            writeDir: write ? defaultResultsDir(root) : undefined,
            repoRoot: root,
          });
          return jsonResult(
            result,
            `Pipeline complete. overall=${result.risk.overall} score=${result.risk.score0to100} conf=${result.confidence?.score0to1 ?? "n/a"}`,
          );
        } catch (e) {
          return text(`error: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
    {
      name: "benchling_pull_candidates",
      label: "Benchling pull",
      description:
        "Pull AA Sequence candidates from Benchling (requires BENCHLING_* env). Use dry_run=true to preview.",
      parameters: benchPullParams,
      async execute(_id, params: Static<typeof benchPullParams>) {
        try {
          const args = ["benchling-pull"];
          if (params.ids) args.push("--ids", params.ids);
          if (params.name_includes) args.push("--name-includes", params.name_includes);
          if (params.dry_run !== false) args.push("--dry-run");
          return jsonResult(runPythonCli(root, args));
        } catch (e) {
          return text(`error: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
    {
      name: "benchling_publish_run",
      label: "Benchling publish",
      description:
        "Publish an immuno-risk run summary to Benchling (idempotent by run_id). Explicit external action — prefer dry_run first.",
      parameters: benchPubParams,
      async execute(_id, params: Static<typeof benchPubParams>) {
        try {
          const args = ["benchling-publish", "--run-dir", params.run_dir];
          if (params.dry_run !== false) args.push("--dry-run");
          return jsonResult(runPythonCli(root, args));
        } catch (e) {
          return text(`error: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
  ];

  return tools;
}
