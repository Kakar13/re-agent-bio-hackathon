# re:AGENT — Pi harness (late immuno-risk stage)

You were started from `harness/`. The **project root is the parent directory** (`..`).

## Mission

This harness is for the **last part** of the de novo binder pipeline (Track C):

**sequence → MHC-I (MHCflurry + IEDB head) + thin MHC-II → Atlas evidence → aggregation → risk + Benchling**

Upstream design (diffusion / MPNN) happens elsewhere. Do **not** reinvent backbone generation here unless `TASK.md` expands scope.

1. Read [`TASK.md`](TASK.md) and [`../docs/IMMUNO_RISK_DESIGN.md`](../docs/IMMUNO_RISK_DESIGN.md).
2. Load [`../skills/reagent/SKILL.md`](../skills/reagent/SKILL.md), then Paperclip / Proto / Boltz as needed.
3. Start with `/denovo` in Pi, or `./repl.sh`. Prefer `run_immuno_pipeline`. Use `benchling_pull_candidates` / `benchling_publish_run` when credentials exist (dry-run first).

## Paths

| What | Where |
| --- | --- |
| Task brief | `TASK.md` |
| Design doc | `../docs/IMMUNO_RISK_DESIGN.md` |
| Python backend | `../src/re_agent/immuno_risk/` |
| Skills | `../skills/` |
| Secrets | `../.env` |
| Artifacts | `../results/immuno_risk/<run-id>/` |
| Setup check | `cd .. && uv run python scripts/check_setup.py` |

## Rules

Follow parent [`../AGENTS.md`](../AGENTS.md). Screening scores only — not clinical probability. Keep MHC-I, MHC-II, Atlas, and aggregation separate. Never reintroduce `stub_hash_rank_v0`. Ask before commits.
