/**
 * Smoke-test immuno-risk tools without launching the Pi TUI.
 * Usage: cd harness && npx tsx scripts/smoke-immuno-tools.ts
 */
import { createImmunoRiskTools } from "../pi-tools.ts";
import { runImmunoPipeline, defaultResultsDir } from "../tools/pipeline.ts";
import { resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));

const DEMO =
  ">demo_natural\n" +
  "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGDPDVLT";

async function main() {
  const root = resolve(here, "../..");
  const tools = createImmunoRiskTools({ repoRoot: root });
  console.log(
    "tools:",
    tools.map((t) => t.name).join(", "),
  );

  process.env.IMMUNO_ALLOW_HEURISTIC_MHC ??= "1";
  const result = runImmunoPipeline(DEMO, {
    writeDir: defaultResultsDir(root),
    mhcClass: "I",
    repoRoot: root,
  });
  console.log(
    JSON.stringify(
      {
        sequenceId: result.sequenceId,
        runId: result.runId,
        cleavages: result.cleavages?.length ?? 0,
        peptides: Array.isArray(result.peptides) ? result.peptides.length : 0,
        mhcBinders: result.mhc.filter((h) => h.binder || h.binderStub).length,
        risk: result.risk.overall,
        score: result.risk.score0to100,
        confidence: result.confidence?.score0to1,
        aggregation: result.aggregation?.overall,
        artifactDir: result.artifactDir,
        writtenPath: result.writtenPath,
      },
      null,
      2,
    ),
  );
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
