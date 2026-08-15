/**
 * Smoke-test immuno-risk tools without launching the Pi TUI.
 * Usage: cd harness && npx tsx scripts/smoke-immuno-tools.ts
 */
import { createImmunoRiskTools } from "../pi-tools.ts";
import { runImmunoPipeline, defaultResultsDir } from "../tools/pipeline.ts";
import { resolve } from "node:path";

const DEMO =
  ">demo_natural\n" +
  "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGDPDVLT";

async function main() {
  const tools = createImmunoRiskTools({ repoRoot: resolve(import.meta.dirname, "../..") });
  console.log(
    "tools:",
    tools.map((t) => t.name).join(", "),
  );

  const root = resolve(import.meta.dirname, "../..");
  const result = runImmunoPipeline(DEMO, {
    writeDir: defaultResultsDir(root),
    mhcClass: "I",
  });
  console.log(
    JSON.stringify(
      {
        sequenceId: result.sequenceId,
        cleavages: result.cleavages.length,
        peptides: result.peptides.length,
        mhcBinders: result.mhc.filter((h) => h.binderStub).length,
        risk: result.risk,
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
