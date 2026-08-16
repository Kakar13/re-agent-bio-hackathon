/** Run real provider artifacts on the checked-in FASTA controls. */
import { readFileSync, readdirSync } from "node:fs";
import { join, resolve } from "node:path";
import { createImmunoRiskTools } from "../pi-tools.ts";

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
  const responseArtifacts = process.env.RESPONSE_ARTIFACTS?.split(",").filter(Boolean) ?? [];
  const defaultResponseAdapter = process.env.DEFAULT_RESPONSE_ADAPTER;
  const netmhciipanArtifact = process.env.NETMHCIIPAN_ARTIFACT;
  if (!responseArtifacts.length || !defaultResponseAdapter || !netmhciipanArtifact) {
    throw new Error(
      "set RESPONSE_ARTIFACTS, DEFAULT_RESPONSE_ADAPTER, and NETMHCIIPAN_ARTIFACT; " +
        "the real-query runner will not fall back to placeholder scores",
    );
  }
  const files = readdirSync(FIXTURES).filter((f) => f.endsWith(".fasta") && !f.includes("spike_sars2.fasta"));
  // prefer n200 spike over full spike
  const fastaFiles = files.filter((f) => f !== "spike_sars2.fasta");

  for (const file of fastaFiles.sort()) {
    const raw = readFileSync(join(FIXTURES, file), "utf8").trim();
    if (!raw || raw.length < 20) {
      console.log(`skip empty: ${file}`);
      continue;
    }
    const id = file.replace(/\.fasta$/, "");
    console.log(`\n=== ${id} ===`);
    const result = await callTool(tools, "run_immuno_pipeline", {
      sequence: raw,
      sequence_id: id,
      response_artifacts: responseArtifacts,
      default_response_adapter: defaultResponseAdapter,
      netmhciipan_artifact: netmhciipanArtifact,
    });
    console.log(result);
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
