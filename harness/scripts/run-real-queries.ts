/**
 * Run immuno-risk tools on real UniProt / known sequences.
 * Usage: cd harness && npx tsx scripts/run-real-queries.ts
 */
import { readFileSync, readdirSync } from "node:fs";
import { join, resolve } from "node:path";
import { createImmunoRiskTools } from "../pi-tools.ts";
import { defaultResultsDir, runImmunoPipeline } from "../tools/pipeline.ts";

const ROOT = resolve(import.meta.dirname, "../..");
const FIXTURES = join(ROOT, "data/raw/immuno_tests");

async function callTool(
  tools: ReturnType<typeof createImmunoRiskTools>,
  name: string,
  params: Record<string, unknown>,
) {
  const tool = tools.find((t) => t.name === name);
  if (!tool) throw new Error(`missing tool ${name}`);
  const result = await tool.execute("test", params as never);
  const text = result.content
    .filter((c) => c.type === "text")
    .map((c) => c.text ?? "")
    .join("\n");
  return text;
}

async function main() {
  const tools = createImmunoRiskTools({ repoRoot: ROOT });
  const files = readdirSync(FIXTURES).filter((f) => f.endsWith(".fasta") && !f.includes("spike_sars2.fasta"));
  // prefer n200 spike over full spike
  const fastaFiles = files.filter((f) => f !== "spike_sars2.fasta");

  console.log("=== list_catalytic_sites ===");
  const sitesRaw = await callTool(tools, "list_catalytic_sites", {});
  const sitesJson = JSON.parse(sitesRaw.split("\n\n").slice(1).join("\n\n") || sitesRaw);
  console.log(`sites: ${sitesJson.count}`);

  const summary: Array<Record<string, unknown>> = [];

  for (const file of fastaFiles.sort()) {
    const raw = readFileSync(join(FIXTURES, file), "utf8").trim();
    if (!raw || raw.length < 20) {
      console.log(`skip empty: ${file}`);
      continue;
    }
    const id = file.replace(/\.fasta$/, "");
    console.log(`\n=== ${id} ===`);

    // Staged query: cleavage only (real motif scan)
    const cleavageText = await callTool(tools, "predict_cleavage", { sequence: raw });
    const jsonStart = cleavageText.indexOf("{");
    const cleavage =
      jsonStart >= 0
        ? (JSON.parse(cleavageText.slice(jsonStart)) as { cleavages?: unknown[]; peptideCount?: number })
        : {};
    console.log(
      `cleavage tool: events=${cleavage.cleavages?.length ?? "?"} peptides=${cleavage.peptideCount ?? "?"}`,
    );

    // End-to-end with MHC I focus (intracellular delivery path)
    const result = runImmunoPipeline(raw, {
      sequenceId: id,
      mhcClass: "I",
      writeDir: defaultResultsDir(ROOT),
    });

    const row = {
      id,
      length: result.sequence.length,
      cleavageEvents: result.cleavages.length,
      peptides: result.peptides.length,
      mhcStubBinders: result.mhc.filter((h) => h.binderStub).length,
      foreignLike: result.tolerance.filter((t) => t.status === "foreign_like").length,
      selfLike: result.tolerance.filter((t) => t.status === "self_like").length,
      risk: result.risk.overall,
      score: result.risk.score0to100,
      topMhc: result.mhc.slice(0, 3).map((h) => `${h.peptide}@${h.allele} rank=${h.rankPctStub}`),
      written: result.writtenPath,
    };
    summary.push(row);
    console.log(JSON.stringify(row, null, 2));
  }

  console.log("\n=== SUMMARY ===");
  console.table(
    summary.map((r) => ({
      id: r.id,
      len: r.length,
      cuts: r.cleavageEvents,
      peptides: r.peptides,
      mhcI_stub: r.mhcStubBinders,
      risk: r.risk,
      score: r.score,
    })),
  );
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
