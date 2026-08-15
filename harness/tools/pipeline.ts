import { mkdirSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { predictCleavage, peptidePool } from "./cleavage.ts";
import { estimateStructureFeatures } from "./features.ts";
import { scoreMhc } from "./mhc.ts";
import { scoreImmunoRisk } from "./risk.ts";
import { checkTolerance } from "./tolerance.ts";
import { parseFastaOrSequence } from "./text.ts";
import type { PipelineResult } from "./types.ts";

export function runImmunoPipeline(rawSequence: string, opts?: {
  sequenceId?: string;
  siteIds?: string[];
  mhcClass?: "I" | "II" | "both";
  writeDir?: string;
}): PipelineResult {
  const parsed = parseFastaOrSequence(rawSequence);
  const sequenceId = opts?.sequenceId ?? parsed.id;
  const sequence = parsed.sequence;
  if (sequence.length < 8) {
    throw new Error(`sequence too short (${sequence.length}); need ≥8 AA`);
  }

  const features = estimateStructureFeatures(sequenceId, sequence);
  const cleavages = predictCleavage(sequence, opts?.siteIds);
  const peptides = peptidePool(sequence, cleavages);
  const mhc = scoreMhc(peptides, { mhcClass: opts?.mhcClass ?? "both", topN: 40 });
  const topPeptides = [...new Set(mhc.slice(0, 30).map((h) => h.peptide))];
  const tolerance = checkTolerance(topPeptides);
  const risk = scoreImmunoRisk({ sequenceId, features, mhc, tolerance });

  const result: PipelineResult = {
    sequenceId,
    sequence,
    features,
    cleavages,
    peptides: peptides.slice(0, 200),
    mhc,
    tolerance,
    risk,
  };

  if (opts?.writeDir) {
    mkdirSync(opts.writeDir, { recursive: true });
    const path = join(opts.writeDir, `${sequenceId.replace(/[^\w.-]+/g, "_")}.json`);
    writeFileSync(path, JSON.stringify(result, null, 2));
    result.writtenPath = path;

    const md = join(opts.writeDir, `${sequenceId.replace(/[^\w.-]+/g, "_")}.md`);
    writeFileSync(
      md,
      [
        `# Immuno-risk: ${sequenceId}`,
        "",
        `- Overall: **${risk.overall}** (score ${risk.score0to100}/100)`,
        `- Length: ${sequence.length}`,
        `- Cleavage events: ${cleavages.length}`,
        `- MHC stub binders: ${mhc.filter((h) => h.binderStub).length}`,
        `- Method caveats: stubs for MHC + tolerance — see JSON`,
        "",
        "## Factors",
        ...risk.factors.map((f) => `- ${f.name} (+${f.contribution}): ${f.note}`),
        "",
        `Full JSON: \`${path}\``,
        "",
      ].join("\n"),
    );
  }

  return result;
}

export function defaultResultsDir(repoRoot: string): string {
  return join(repoRoot, "results", "immuno_risk");
}

export function ensureParent(filePath: string): void {
  mkdirSync(dirname(filePath), { recursive: true });
}
