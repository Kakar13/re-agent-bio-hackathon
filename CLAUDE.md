# re:AGENT project

Hackathon repo for scientific agents. Weekend of Aug 15–16, 2026. Submission Sunday 10:45 AM.

Primary agent harness: **[Pi](https://pi.dev/)** from **[`harness/`](harness/)** (de novo task). Project instructions also live in `AGENTS.md` (Pi loads both). **Keep `AGENTS.md` and this file aligned.**

## First thing for humans (and for you)

**Onboarding entry point is always [`START_HERE.md`](START_HERE.md).**

If the user asks how to set up, clone, init, install tools, get keys, run checks, or “what do I do first”:

1. Point them to **`START_HERE.md`** (read it; summarize the next 1–2 steps if helpful).
2. For deeper troubleshooting, use **`docs/SETUP.md`**.
3. For partners / credits / booths, use **`SPONSORS.md`**.
4. Do **not** invent a parallel setup path. Do **not** invent API keys.

Each teammate uses their own `.env` and Paperclip login. Never copy keys between laptops.

## Skills (read before using a tool)

Sundial re:AGENT skills live in [`skills/`](skills/). Read **`skills/reagent/SKILL.md` first** when scoping the weekend, then the matching tool skill before you run anything:

| Skill | When |
| --- | --- |
| [`skills/reagent/`](skills/reagent/SKILL.md) | Track choice, judging bar, weekend plan |
| [`skills/paperclip/`](skills/paperclip/SKILL.md) | Literature, FDA, trials, patents |
| [`skills/cellxgene-census/`](skills/cellxgene-census/SKILL.md) | Single-cell / Census queries |
| [`skills/proto/`](skills/proto/SKILL.md) | Structure, design, docking, bioinformatics tools |
| [`skills/boltz/`](skills/boltz/SKILL.md) | Fold / affinity with Boltz-2 |

Project skills are also linked at `.claude/skills` → `skills/` for Claude Code auto-discovery. Follow each skill’s rules exactly.

## Tools to use

- **Paperclip** — literature, FDA, trials, UniProt/PDB/ChEMBL. Prefer MCP (Pi: [pi-mcp-adapter](https://pi.dev/packages/pi-mcp-adapter) + `.mcp.json`; Cursor: `.cursor/mcp.json`) or the CLI (`paperclip search|map|grep|sql|cat`). Load `/paperclip` or `skills/paperclip` before Paperclip work. Do not download entire papers by hand.
- **Proto** — two paths: (1) hosted MCP `proto-bio` at `https://mcp.evodesign.org/mcp` (`PROTO_API_KEY` + `proto_client`); (2) local `proto-tools` / `proto-language` via `uv sync --extra proto`, run on Modal with `device="modal"` after `uv run modal setup` and deploy. Do not invent Modal tokens.
- **CELLxGENE Census / Boltz** — no API key; follow the skills above for install and gotchas.
- **Claude API** — `ANTHROPIC_API_KEY` in `.env` (also works for Pi auth).

If a check fails, tell the user which booth / Discord channel to get the credential from (`SPONSORS.md`).

### Pi + MCP

1. `cd harness && npm install && npx pi install -l --approve npm:pi-mcp-adapter` ([package](https://pi.dev/packages/pi-mcp-adapter)).
2. Shared servers live in **[`harness/.mcp.json`](harness/.mcp.json)** (and root `.mcp.json` for other tools). Export `PROTO_API_KEY` before starting Pi.
3. Inside Pi: `/denovo`, `/mcp`, `/mcp-auth <server>` if needed. Lazy proxy keeps context small.

Task brief for this team: [`harness/TASK.md`](harness/TASK.md).

## Working rules

- Keep secrets in `.env` only. Never commit `.env`.
- Write durable code under `src/re_agent/`. Scratch work goes in `notebooks/` or `scripts/`.
- Large artifacts (FASTA, PDB, paper dumps) go in `data/` or `results/` — those dirs are gitignored except `.gitkeep`.
- Prefer inspectable structured outputs (JSON/CSV/Markdown tables) over chatty prose. Judges need to see the reasoning.
- Paperclip workflow that works: `search` → `map` over the result id → `reduce` / synthesize. Do not download entire papers by hand.
- After setup changes, run `uv run python scripts/check_setup.py`.

## Tracks (pick one and stay on it)

- A: end-to-end scientific agent with tools + inspectable trace
- B: dataset / meta-analysis with a finding no single paper shows
- C: biological design to a spec, plus evidence it could hold up
