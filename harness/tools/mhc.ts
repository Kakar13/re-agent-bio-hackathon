/**
 * MHC scoring — calls Python immuno_risk backend (MHCflurry / HLAIIPred).
 * No stub_hash_rank_v0 in user-facing outputs.
 */
import { spawnSync } from "node:child_process";
import { existsSync } from "node:fs";
import { resolve } from "node:path";
import type { MhcHit } from "./types.ts";

const DEFAULT_MHC1 = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"];
const DEFAULT_MHC2 = ["DRB1*01:01", "DRB1*15:01"];

function repoRootGuess(): string {
  // harness/tools -> repo root
  const here = new URL(".", import.meta.url).pathname;
  return resolve(here, "..", "..");
}

function pythonCmd(root: string): { cmd: string; argsPrefix: string[] } {
  const candidates = [
    resolve(root, ".venv/bin/python"),
    "/tmp/immuno-venv/bin/python",
    process.env.IMMUNO_PYTHON ?? "",
  ].filter(Boolean);
  for (const c of candidates) {
    if (existsSync(c)) return { cmd: c, argsPrefix: [] };
  }
  // Avoid `uv run` when it would sync optional proto git deps in restricted envs
  return { cmd: "python3", argsPrefix: [] };
}

/** Run full Python pipeline bridge for a sequence (preferred). */
export function scoreSequenceViaPython(
  sequence: string,
  opts?: {
    sequenceId?: string;
    mhcClass?: "I" | "II" | "both";
    allelesI?: string[];
    allelesIi?: string[];
    write?: boolean;
    repoRoot?: string;
  },
): Record<string, unknown> {
  const root = opts?.repoRoot ?? repoRootGuess();
  const payload = {
    action: "run",
    sequence,
    sequence_id: opts?.sequenceId ?? "query",
    mhc_class: opts?.mhcClass ?? "both",
    alleles_i: opts?.allelesI,
    alleles_ii: opts?.allelesIi,
    write: opts?.write ?? false,
  };
  const env = {
    ...process.env,
    PYTHONPATH: `${root}/src${process.env.PYTHONPATH ? `:${process.env.PYTHONPATH}` : ""}`,
    IMMUNO_ALLOW_HEURISTIC_MHC: process.env.IMMUNO_ALLOW_HEURISTIC_MHC ?? "1",
  };
  const py = pythonCmd(root);
  const proc = spawnSync(
    py.cmd,
    [...py.argsPrefix, "-m", "re_agent.immuno_risk.bridge"],
    {
      cwd: root,
      input: JSON.stringify(payload),
      encoding: "utf-8",
      env,
      maxBuffer: 32 * 1024 * 1024,
    },
  );
  if (proc.error) {
    throw new Error(`Python bridge failed to start: ${proc.error.message}`);
  }
  if (proc.status !== 0) {
    throw new Error(
      `Python bridge exited ${proc.status}: ${(proc.stderr || proc.stdout || "").slice(0, 800)}`,
    );
  }
  return JSON.parse(proc.stdout) as Record<string, unknown>;
}

function mapHit(raw: Record<string, unknown>): MhcHit {
  const pct = (raw.percentile_rank as number | null | undefined) ?? null;
  const binder = Boolean(raw.binder);
  return {
    peptide: String(raw.peptide),
    allele: String(raw.allele),
    mhcClass: raw.mhc_class === "II" ? "II" : "I",
    length: Number(raw.length ?? String(raw.peptide).length),
    start: (raw.start as number | null | undefined) ?? null,
    end: (raw.end as number | null | undefined) ?? null,
    affinityNm: (raw.affinity_nm as number | null | undefined) ?? null,
    presentationScore: (raw.presentation_score as number | null | undefined) ?? null,
    processingScore: (raw.processing_score as number | null | undefined) ?? null,
    percentileRank: pct,
    binder,
    method: String(raw.method ?? "python_backend"),
    version: (raw.version as string | null | undefined) ?? null,
    caveat: String(raw.caveat ?? ""),
    // compat for older risk code paths
    rankPctStub: pct ?? undefined,
    binderStub: binder,
  };
}

/**
 * Score peptides for MHC. Prefer scoreSequenceViaPython for full evidence.
 * Peptide-only path runs a short sequence join (peptides concatenated with GGGS)
 * only when a full sequence is not available — prefer pipeline.
 */
export function scoreMhc(
  peptides: string[],
  opts?: { mhcClass?: "I" | "II" | "both"; alleles?: string[]; topN?: number; sequence?: string },
): MhcHit[] {
  const seq =
    opts?.sequence ??
    peptides
      .filter((p) => p.length >= 8)
      .slice(0, 30)
      .join("");
  if (seq.length < 8) return [];

  const allelesI =
    opts?.alleles?.filter((a) => /HLA-[ABC]\*/.test(a) || /^[ABC]\*/.test(a)) ?? DEFAULT_MHC1;
  const allelesIi =
    opts?.alleles?.filter((a) => /DR|DQ|DP/.test(a)) ?? DEFAULT_MHC2;

  try {
    const result = scoreSequenceViaPython(seq, {
      mhcClass: opts?.mhcClass ?? "both",
      allelesI,
      allelesIi,
      write: false,
    });
    const peptidesEv = (result.peptides as Array<Record<string, unknown>>) ?? [];
    const hits = peptidesEv.map((e) => mapHit((e.mhc as Record<string, unknown>) ?? e));
    const wanted = new Set(peptides.map((p) => p.toUpperCase()));
    const filtered = wanted.size
      ? hits.filter((h) => wanted.has(h.peptide.toUpperCase()))
      : hits;
    const topN = opts?.topN ?? 40;
    return filtered.slice(0, topN);
  } catch (e) {
    throw new Error(
      `score_mhc requires Python immuno backend: ${e instanceof Error ? e.message : String(e)}`,
    );
  }
}
