/** Pi tools backed by the real Python adapter pipeline, never placeholder scores. */
import { existsSync, readFileSync } from "node:fs";
import { resolve } from "node:path";
import { Type, type Static } from "typebox";
import { runRealImmunoPipeline } from "./tools/real-pipeline.ts";
import { jsonResult, text } from "./tools/text.ts";

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

  const statusParams = Type.Object({});
  const pipelineParams = Type.Object({
    sequence: Type.String({ description: "FASTA or raw amino-acid sequence" }),
    sequence_id: Type.Optional(Type.String()),
    response_artifacts: Type.Array(Type.String(), {
      minItems: 1,
      description: "Versioned JSON artifacts produced by response-model adapters",
    }),
    default_response_adapter: Type.String(),
    calibrations: Type.Optional(Type.Array(Type.String(), {
      description: "Optional ADAPTER_ID=PATH frozen calibration artifacts",
    })),
    netmhciipan_artifact: Type.Optional(Type.String({
      description: "Versioned NetMHCIIpan cache retaining separate EL and BA columns",
    })),
    iedb_live: Type.Optional(Type.Boolean({
      description: "Call documented IEDB NetMHCIIpan 4.3 EL and BA APIs with local cache",
    })),
    challengers: Type.Optional(
      Type.Array(Type.String(), {
        description: "Optional PROVIDER=VERSION=PATH challenger artifacts",
      }),
    ),
    accessibility_artifact: Type.Optional(Type.String()),
    cleavage_artifact: Type.Optional(Type.String()),
    shared_hla_artifact: Type.Optional(Type.String({
      description: "Matched-self/query shared-HLA evidence keyed by sequence hash",
    })),
  });

  const tools: ImmunoToolDefinition[] = [
    {
      name: "immuno_architecture_status",
      label: "Immuno architecture status",
      description: "Report the frozen HLA panel, fusion rule, benchmark, and real provider artifacts.",
      parameters: statusParams,
      async execute() {
        const paths = {
          hlaPanel: resolve(root, "docs", "hla_class_ii_panel.v1.json"),
          fusionRule: resolve(root, "docs", "fusion_rule.v1.json"),
          benchmarkManifest: resolve(
            root,
            "data",
            "processed",
            "benchmarks",
            "iedb-class-ii-response-v1.manifest.json",
          ),
          selfProteome: resolve(root, "data", "processed", "self_proteome.parquet"),
        };
        const available = Object.fromEntries(
          Object.entries(paths).map(([key, path]) => [key, { path, exists: existsSync(path) }]),
        );
        const benchmark = existsSync(paths.benchmarkManifest)
          ? JSON.parse(readFileSync(paths.benchmarkManifest, "utf8"))
          : null;
        return jsonResult({
          available,
          benchmark,
          registeredTools: ["immuno_architecture_status", "run_immuno_pipeline"],
          excludedPlaceholders: [
            "hash-derived MHC ranks",
            "tiny hard-coded self peptide set",
            "weighted stub risk score",
          ],
        });
      },
    },
    {
      name: "run_immuno_pipeline",
      label: "Run immuno pipeline",
      description:
        "Run response-model adapters, NetMHCIIpan EL+BA, challenger providers, self-tolerance evidence, and transparent late fusion. Requires real versioned provider artifacts.",
      parameters: pipelineParams,
      async execute(_id, params: Static<typeof pipelineParams>, signal?: AbortSignal) {
        try {
          if (Boolean(params.netmhciipan_artifact) === Boolean(params.iedb_live)) {
            return text("error: provide exactly one of netmhciipan_artifact or iedb_live=true");
          }
          const result = await runRealImmunoPipeline(params.sequence, {
            repoRoot: root,
            sequenceId: params.sequence_id,
            responseArtifacts: params.response_artifacts,
            defaultResponseAdapter: params.default_response_adapter,
            calibrations: params.calibrations,
            netmhciipanArtifact: params.netmhciipan_artifact,
            iedbLive: params.iedb_live,
            challengers: params.challengers,
            accessibilityArtifact: params.accessibility_artifact,
            cleavageArtifact: params.cleavage_artifact,
            sharedHlaArtifact: params.shared_hla_artifact,
          }, signal);
          if (result.exitCode !== 0) {
            return text(`pipeline failed (${result.exitCode}): ${result.stderr || result.stdout}`);
          }
          return jsonResult(
            result,
            `Pipeline complete. Inspectable artifact: ${result.outputPath}`,
          );
        } catch (e) {
          return text(`error: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
  ];

  return tools;
}
