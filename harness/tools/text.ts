/** Result helpers inspired by AutopsyAI qm/src/harness/pi-tools.ts */

const MAX_TOOL_RESULT_CHARS = 80_000;

export function text(s: string) {
  return { content: [{ type: "text" as const, text: s }], details: {} };
}

export function jsonResult(value: unknown, preamble?: string) {
  const body = JSON.stringify(value, null, 2);
  const full = preamble ? `${preamble}\n\n${body}` : body;
  const capped =
    full.length > MAX_TOOL_RESULT_CHARS
      ? `${full.slice(0, MAX_TOOL_RESULT_CHARS)}\n…[truncated ${full.length} chars]`
      : full;
  return {
    content: [{ type: "text" as const, text: capped }],
    details: value as Record<string, unknown>,
  };
}

export function cleanSequence(raw: string): string {
  return raw
    .replace(/^>.*$/gm, "")
    .replace(/[^A-Za-z]/g, "")
    .toUpperCase()
    .replace(/[^ACDEFGHIKLMNPQRSTVWY]/g, "");
}

export function parseFastaOrSequence(raw: string): { id: string; sequence: string } {
  const trimmed = raw.trim();
  if (trimmed.startsWith(">")) {
    const lines = trimmed.split(/\r?\n/);
    const id = (lines[0] ?? ">seq").replace(/^>/, "").trim() || "seq";
    const sequence = cleanSequence(lines.slice(1).join(""));
    return { id, sequence };
  }
  return { id: "seq", sequence: cleanSequence(trimmed) };
}
