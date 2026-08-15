/**
 * Pi extension: register immuno-risk tools for the late-stage de novo pipeline.
 *
 * Pattern: AutopsyAI createPiTools → customTools; here we register the same
 * tool definitions into the interactive Pi session via pi.registerTool.
 */
import type { ExtensionAPI } from "@earendil-works/pi-coding-agent";
import { createImmunoRiskTools } from "../../pi-tools.ts";

export default function immunoRiskToolsExtension(pi: ExtensionAPI) {
  const tools = createImmunoRiskTools();
  for (const tool of tools) {
    // registerTool accepts the same shape as AutopsyAI customTools / defineTool output
    pi.registerTool(tool as Parameters<ExtensionAPI["registerTool"]>[0]);
  }

  pi.on("session_start", async (_event, ctx) => {
    if (ctx.mode === "tui") {
      ctx.ui.notify(
        `Immuno-risk tools ready (${tools.map((t) => t.name).join(", ")}). Try /denovo then run_immuno_pipeline.`,
        "info",
      );
    }
  });

  pi.registerCommand("immuno-tools", {
    description: "List immuno-risk custom tools",
    handler: async (_args, ctx) => {
      const lines = tools.map((t) => `- ${t.name}: ${t.description}`);
      ctx.ui.notify(["Immuno-risk tools:", ...lines].join("\n"), "info");
    },
  });
}
