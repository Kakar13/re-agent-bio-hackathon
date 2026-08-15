/**
 * Interactive immuno-risk REPL with chat context + tool calling.
 *
 * Usage:
 *   cd harness && ./repl.sh
 *   cd harness && npm run repl
 *
 * Slash commands: /help /clear /history /tools /quit
 */
import { createInterface } from "node:readline/promises";
import { stdin as input, stdout as output } from "node:process";
import { existsSync, readFileSync } from "node:fs";
import { resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import Anthropic from "@anthropic-ai/sdk";
import { createImmunoRiskTools, type ImmunoToolDefinition } from "./pi-tools.ts";

const __dirname = dirname(fileURLToPath(import.meta.url));
const HARNESS_ROOT = __dirname;
const REPO_ROOT = resolve(HARNESS_ROOT, "..");

const MODEL = process.env.IMMUNO_REPL_MODEL || "claude-sonnet-4-20250514";
const MAX_TOKENS = 4096;
const MAX_TOOL_ROUNDS = 12;

function loadEnvFiles(): void {
  for (const path of [resolve(REPO_ROOT, ".env"), resolve(HARNESS_ROOT, ".env")]) {
    if (!existsSync(path)) continue;
    for (const line of readFileSync(path, "utf8").split(/\r?\n/)) {
      const m = line.match(/^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$/);
      if (!m) continue;
      const key = m[1]!;
      let val = m[2]!.trim();
      if (
        (val.startsWith('"') && val.endsWith('"')) ||
        (val.startsWith("'") && val.endsWith("'"))
      ) {
        val = val.slice(1, -1);
      }
      if (process.env[key] === undefined) process.env[key] = val;
    }
  }
}

function toInputSchema(parameters: unknown): Anthropic.Tool["input_schema"] {
  const raw = JSON.parse(JSON.stringify(parameters ?? { type: "object", properties: {} })) as Record<
    string,
    unknown
  >;
  if (raw.type !== "object") {
    return { type: "object", properties: {} };
  }
  return raw as Anthropic.Tool["input_schema"];
}

function toolsForAnthropic(tools: ImmunoToolDefinition[]): Anthropic.Tool[] {
  return tools.map((t) => ({
    name: t.name,
    description: t.description,
    input_schema: toInputSchema(t.parameters),
  }));
}

function toolResultText(
  result: Awaited<ReturnType<ImmunoToolDefinition["execute"]>>,
): string {
  return result.content
    .filter((c) => c.type === "text")
    .map((c) => c.text ?? "")
    .join("\n");
}

const SYSTEM = `You are the re:AGENT immuno-risk REPL assistant (late-stage de novo binder pipeline).

You help the user run and interpret:
sequence (+ structure) → accessibility features → cleavage vs ~10 catalytic sites → peptide pool → MHC I/II → tolerance → risk score.

You have custom tools for that pipeline. Prefer tools over guessing sequences or scores.
- MHC and tolerance tools are STUBS until NetMHCpan / HLA Ligand Atlas are wired — always say so when reporting binder/tolerance claims.
- Cleavage motifs are a curated starter set (~10 sites), not a full protease proteome.
- Write durable artifacts under results/immuno_risk/ via run_immuno_pipeline (write=true).
- Keep answers concise. Cite tool outputs. Maintain continuity with earlier turns in this chat.
- If the user pastes a FASTA or raw sequence, run the pipeline (or the specific stage they ask for).`;

function printHelp(): void {
  console.log(`
Commands:
  /help      Show this help
  /clear     Clear chat context (start fresh)
  /history   Show turn count + last user/assistant snippets
  /tools     List immuno-risk tools
  /quit      Exit

Anything else is sent to Claude with full chat context and tool access.
Paste a FASTA or ask e.g. "run immuno pipeline on insulin" .
`);
}

async function main(): Promise<void> {
  loadEnvFiles();
  const apiKey = process.env.ANTHROPIC_API_KEY?.trim();
  if (!apiKey) {
    console.error("Missing ANTHROPIC_API_KEY. Set it in ../.env or harness/.env");
    process.exit(1);
  }

  const immunoTools = createImmunoRiskTools({ repoRoot: REPO_ROOT });
  const toolByName = new Map(immunoTools.map((t) => [t.name, t]));
  const anthropicTools = toolsForAnthropic(immunoTools);
  const client = new Anthropic({ apiKey });

  type Msg = Anthropic.MessageParam;
  const messages: Msg[] = [];

  const rl = createInterface({ input, output, terminal: true });
  console.log(`Immuno-risk REPL  model=${MODEL}`);
  console.log(`Repo root: ${REPO_ROOT}`);
  console.log(`Tools: ${immunoTools.map((t) => t.name).join(", ")}`);
  console.log(`Type /help for commands. Chat context is kept until /clear.\n`);

  const handleSlash = async (line: string): Promise<boolean> => {
    const cmd = line.trim().toLowerCase();
    if (cmd === "/help" || cmd === "/h" || cmd === "?") {
      printHelp();
      return true;
    }
    if (cmd === "/quit" || cmd === "/exit" || cmd === "/q") {
      rl.close();
      process.exit(0);
    }
    if (cmd === "/clear") {
      messages.length = 0;
      console.log("(context cleared)\n");
      return true;
    }
    if (cmd === "/tools") {
      for (const t of immunoTools) console.log(`- ${t.name}: ${t.description.slice(0, 100)}…`);
      console.log();
      return true;
    }
    if (cmd === "/history") {
      console.log(`messages in context: ${messages.length}`);
      for (const m of messages.slice(-6)) {
        const role = m.role;
        const preview =
          typeof m.content === "string"
            ? m.content.slice(0, 120)
            : Array.isArray(m.content)
              ? m.content
                  .map((b) => {
                    if (b.type === "text") return b.text.slice(0, 80);
                    if (b.type === "tool_use") return `[tool_use ${b.name}]`;
                    if (b.type === "tool_result") return `[tool_result]`;
                    return `[${(b as { type: string }).type}]`;
                  })
                  .join(" ")
              : "";
        console.log(`  ${role}: ${preview}`);
      }
      console.log();
      return true;
    }
    return false;
  };

  const runTurn = async (userText: string): Promise<void> => {
    messages.push({ role: "user", content: userText });

    for (let round = 0; round < MAX_TOOL_ROUNDS; round++) {
      process.stdout.write(round === 0 ? "… thinking\n" : "… tools\n");
      const response = await client.messages.create({
        model: MODEL,
        max_tokens: MAX_TOKENS,
        system: SYSTEM,
        tools: anthropicTools,
        messages,
      });

      messages.push({ role: "assistant", content: response.content });

      const toolUses = response.content.filter(
        (b): b is Anthropic.ToolUseBlock => b.type === "tool_use",
      );
      const texts = response.content.filter((b): b is Anthropic.TextBlock => b.type === "text");

      for (const t of texts) {
        if (t.text.trim()) console.log(`\nassistant> ${t.text.trim()}\n`);
      }

      if (response.stop_reason !== "tool_use" || toolUses.length === 0) {
        return;
      }

      const toolResults: Anthropic.ToolResultBlockParam[] = [];
      for (const use of toolUses) {
        const tool = toolByName.get(use.name);
        console.log(`  ↳ ${use.name}(${JSON.stringify(use.input).slice(0, 100)}…)`);
        if (!tool) {
          toolResults.push({
            type: "tool_result",
            tool_use_id: use.id,
            content: `error: unknown tool ${use.name}`,
            is_error: true,
          });
          continue;
        }
        try {
          const out = await tool.execute(use.id, use.input as never);
          const body = toolResultText(out);
          const capped = body.length > 60_000 ? `${body.slice(0, 60_000)}\n…[truncated]` : body;
          toolResults.push({
            type: "tool_result",
            tool_use_id: use.id,
            content: capped,
          });
        } catch (e) {
          toolResults.push({
            type: "tool_result",
            tool_use_id: use.id,
            content: `error: ${e instanceof Error ? e.message : String(e)}`,
            is_error: true,
          });
        }
      }
      messages.push({ role: "user", content: toolResults });
    }
    console.log("(stopped: too many tool rounds — say /clear or continue)\n");
  };

  while (true) {
    let line: string;
    try {
      line = await rl.question("you> ");
    } catch {
      break;
    }
    const trimmed = line.trim();
    if (!trimmed) continue;
    if (await handleSlash(trimmed)) continue;
    try {
      await runTurn(trimmed);
    } catch (e) {
      console.error(`error: ${e instanceof Error ? e.message : String(e)}\n`);
    }
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
