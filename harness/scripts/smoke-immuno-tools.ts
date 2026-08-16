/** Smoke-test real tool registration without requiring provider artifacts. */
import { createImmunoRiskTools } from "../pi-tools.ts";
import { resolve } from "node:path";

async function main() {
  const tools = createImmunoRiskTools({ repoRoot: resolve(import.meta.dirname, "../..") });
  const names = tools.map((tool) => tool.name);
  if (names.join(",") !== "immuno_architecture_status,run_immuno_pipeline") {
    throw new Error(`unexpected registered tools: ${names.join(",")}`);
  }
  const status = tools[0];
  const result = await status.execute("smoke", {} as never);
  const rendered = result.content.map((item) => item.text ?? "").join("\n");
  if (!rendered.includes("excludedPlaceholders")) {
    throw new Error("status output did not declare excluded placeholders");
  }
  console.log(`tools: ${names.join(", ")}`);
  console.log("real architecture status: ok");
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
