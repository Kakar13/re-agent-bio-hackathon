# re:AGENT — Pi harness (late immuno-risk stage)

You were started from `harness/`. The **project root is the parent directory** (`..`).

## Mission

This harness is for the **last part** of the de novo binder pipeline (Track C):

**sequence (+ structure) → cleavage peptides → MHC I/II → tolerance → risk score**

Upstream design (diffusion / MPNN) happens elsewhere. Do **not** spend the weekend reinventing backbone generation here unless `TASK.md` explicitly expands scope.

1. Read [`TASK.md`](TASK.md) first and follow that pipeline.
2. Load [`../skills/reagent/SKILL.md`](../skills/reagent/SKILL.md), then Paperclip / Proto / Boltz / ESM as needed.
3. Start or resume with `/denovo` in Pi, or use the text REPL: `./repl.sh` (keeps chat context + same immuno tools). Custom tools are registered by `.pi/extensions/immuno-risk-tools.ts` (see `pi-tools.ts`). Prefer `run_immuno_pipeline` for an end-to-end pass; use the staged tools when debugging one step.

## Paths

| What | Where |
| --- | --- |
| Task brief | `TASK.md` (this folder) |
| Code | `../src/re_agent/` (e.g. `immuno_risk/`) |
| Skills | `../skills/` (also via `.pi/settings.json`) |
| Secrets | `../.env` (never commit; source before launch) |
| Artifacts | `../data/`, `../results/immuno_risk/` |
| Setup check | `cd .. && uv run python scripts/check_setup.py` |

Prefer editing and running commands against the **repo root**, not this `harness/` folder, except when changing Pi settings under `.pi/` or updating `TASK.md`.

## Rules

Follow parent [`../AGENTS.md`](../AGENTS.md). Cite or measure every claim. Keep reasoning inspectable for judges. De novo candidates have **no ground-truth labels** — say so and use natural holdouts + structure checks as proxies.

MCP (Paperclip + Proto) comes from [`.mcp.json`](.mcp.json) via `pi-mcp-adapter`. Use `/mcp` after auth; `/mcp-auth paperclip` / `/mcp-auth proto-bio` if prompted.
