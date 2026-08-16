import { mkdirSync, writeFileSync } from "node:fs";
import { spawn } from "node:child_process";
import { join } from "node:path";
import { parseFastaOrSequence } from "./text.ts";

export interface RealPipelineOptions {
  repoRoot: string;
  sequenceId?: string;
  responseArtifacts: string[];
  defaultResponseAdapter: string;
  calibrations?: string[];
  netmhciipanArtifact?: string;
  iedbLive?: boolean;
  challengers?: string[];
  accessibilityArtifact?: string;
  cleavageArtifact?: string;
  sharedHlaArtifact?: string;
}

export interface CommandResult {
  exitCode: number;
  stdout: string;
  stderr: string;
  outputPath: string;
}

export async function runRealImmunoPipeline(
  rawSequence: string,
  options: RealPipelineOptions,
  signal?: AbortSignal,
): Promise<CommandResult> {
  const parsed = parseFastaOrSequence(rawSequence);
  const sequenceId = options.sequenceId ?? parsed.id;
  const runId = `${Date.now()}-${sequenceId.replace(/[^\w.-]+/g, "_")}`;
  const inputDir = join(options.repoRoot, "results", "immuno_risk", "inputs");
  const outputDir = join(options.repoRoot, "results", "immuno_risk");
  mkdirSync(inputDir, { recursive: true });
  mkdirSync(outputDir, { recursive: true });
  const fastaPath = join(inputDir, `${runId}.fasta`);
  const outputPath = join(outputDir, `${runId}.json`);
  writeFileSync(fastaPath, `>${sequenceId}\n${parsed.sequence}\n`);

  const args = [
    "run",
    "python",
    join(options.repoRoot, "scripts", "screen_candidates.py"),
    "--fasta",
    fastaPath,
    "--default-response-adapter",
    options.defaultResponseAdapter,
    "--output",
    outputPath,
  ];
  if (options.netmhciipanArtifact) {
    args.push("--netmhciipan-artifact", options.netmhciipanArtifact);
  } else if (options.iedbLive) {
    args.push("--iedb-live");
  } else {
    throw new Error("a NetMHCIIpan artifact or iedbLive=true is required");
  }
  for (const path of options.responseArtifacts) {
    args.push("--response-artifact", path);
  }
  for (const calibration of options.calibrations ?? []) {
    args.push("--calibration", calibration);
  }
  for (const challenger of options.challengers ?? []) {
    args.push("--challenger", challenger);
  }
  if (options.accessibilityArtifact) {
    args.push("--accessibility-artifact", options.accessibilityArtifact);
  }
  if (options.cleavageArtifact) {
    args.push("--cleavage-artifact", options.cleavageArtifact);
  }
  if (options.sharedHlaArtifact) {
    args.push("--shared-hla-artifact", options.sharedHlaArtifact);
  }

  return new Promise((resolve, reject) => {
    const child = spawn("uv", args, {
      cwd: options.repoRoot,
      env: { ...process.env, PYTHONPATH: join(options.repoRoot, "src") },
      signal,
    });
    let stdout = "";
    let stderr = "";
    child.stdout.on("data", (chunk) => {
      stdout += chunk.toString();
    });
    child.stderr.on("data", (chunk) => {
      stderr += chunk.toString();
    });
    child.on("error", reject);
    child.on("close", (code) => {
      resolve({ exitCode: code ?? 1, stdout, stderr, outputPath });
    });
  });
}
