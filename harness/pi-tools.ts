/**
 * Immuno-risk custom tools for Pi — pattern inspired by
 * AutopsyAI qm/src/harness/pi-tools.ts (createPiTools factory + typed tools),
 * adapted for re:AGENT late-stage cleavage → MHC → tolerance → risk.
 *
 * Plain tool objects (no defineTool import) so they load in Pi extensions and
 * in local smoke tests without hitting pi-coding-agent package exports.
 */
import { Type, type Static } from "typebox";
import { resolve } from "node:path";
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
  /** Absolute path to repo root (parent of harness/). */
  repoRoot?: string;
}

/** Minimal tool shape accepted by pi.registerTool / createAgentSession customTools. */
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
    peptides: Type.Array(Type.String()),
    mhc_class: Type.Optional(
      Type.Union([Type.Literal("I"), Type.Literal("II"), Type.Literal("both")]),
    ),
    top_n: Type.Optional(Type.Number()),
  });
  const toleranceParams = Type.Object({
    peptides: Type.Array(Type.String()),
  });
  const riskParams = Type.Object({
    sequence_id: Type.String(),
    features: Type.Record(Type.String(), Type.Unknown()),
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

  const tools: ImmunoToolDefinition[] = [
    {
      name: "list_catalytic_sites",
      label: "List catalytic sites",
      description:
        "List the curated protease / catalytic sites used for cleavage prediction (~10 starters). Call before predict_cleavage.",
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
        "Estimate RSA / disorder / loop / unfolding-ΔG proxies for a protein sequence. Sequence-only heuristics until real structure or MPNN scorefiles are provided.",
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
        "Compare a sequence against curated catalytic sites and return cleavage events plus a peptide pool (including MHC-length windows).",
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
        "Score peptides for MHC class I and/or II. CURRENTLY A STUB — replace with NetMHCpan/MHCflurry before trusting binder calls. Prefer class I for intracellular delivery.",
      parameters: mhcParams,
      async execute(_id, params: Static<typeof mhcParams>) {
        if (!params.peptides?.length) return text("error: peptides array required");
        const hits = scoreMhc(params.peptides, {
          mhcClass: params.mhc_class ?? "both",
          topN: params.top_n ?? 40,
        });
        return jsonResult({ hits, binderCount: hits.filter((h) => h.binderStub).length });
      },
    },
    {
      name: "check_tolerance",
      label: "Check tolerance",
      description:
        "Compare peptides to a self / healthy-ligand reference. CURRENTLY a tiny stub set — wire HLA Ligand Atlas before claiming clinical tolerance.",
      parameters: toleranceParams,
      async execute(_id, params: Static<typeof toleranceParams>) {
        if (!params.peptides?.length) return text("error: peptides array required");
        return jsonResult({ hits: checkTolerance(params.peptides) });
      },
    },
    {
      name: "score_immuno_risk",
      label: "Score immuno risk",
      description:
        "Combine structure features + MHC stub hits + tolerance into a polarized risk score (low/moderate/high). Demo weights only.",
      parameters: riskParams,
      async execute(_id, params: Static<typeof riskParams>) {
        const risk = scoreImmunoRisk({
          sequenceId: params.sequence_id,
          features: params.features as unknown as StructureFeatures,
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
        "End-to-end late-stage pipeline: features → cleavage → peptide pool → MHC stub → tolerance → risk. Writes JSON+MD under results/immuno_risk/ when write=true.",
      parameters: pipelineParams,
      async execute(_id, params: Static<typeof pipelineParams>) {
        try {
          const write = params.write !== false;
          const result = runImmunoPipeline(params.sequence, {
            sequenceId: params.sequence_id,
            siteIds: params.site_ids,
            mhcClass: params.mhc_class ?? "both",
            writeDir: write ? defaultResultsDir(root) : undefined,
          });
          return jsonResult(
            result,
            `Pipeline complete. overall=${result.risk.overall} score=${result.risk.score0to100}`,
          );
        } catch (e) {
          return text(`error: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
  ];

  return tools;
}
