# Pi harness (re:AGENT — late immuno-risk stage)

Project-local [Pi](https://pi.dev/) workspace for the **last part** of the de novo binder pipeline: cleavage → MHC → tolerance → risk score. Upstream design (diffusion / MPNN) is out of scope here unless `TASK.md` says otherwise.

## One-time setup

**Node ≥ 22.19** (Pi 0.84+). Check with `node -v`.

From the **repo root**:

```bash
# Node 20+
cd harness
npm install

# Install project packages into .pi/npm (needs project trust)
npx pi install -l --approve npm:pi-mcp-adapter
npx pi install -l --approve npm:@langchain/langsmith-pi-extension
# or: npm run install-packages

# Keys live in the parent .env — do not commit it
set -a && source ../.env && set +a
# Add LANGSMITH_API_KEY + TRACE_TO_LANGSMITH=true (see ../.env.example)
```

## LangSmith traces

Pi sessions are traced to LangSmith via [`@langchain/langsmith-pi-extension`](https://docs.langchain.com/langsmith/trace-with-pi).

1. Create an API key: [LangSmith API keys](https://docs.langchain.com/langsmith/create-account-api-key)
2. In `../.env` (your values are fine):
   ```bash
   LANGSMITH_TRACING=true
   LANGSMITH_ENDPOINT=https://api.smith.langchain.com
   LANGSMITH_API_KEY=lsv2_pt_...
   LANGSMITH_PROJECT=reAgent-hackathon
   ```
   `./run.sh` maps `LANGSMITH_TRACING` → `TRACE_TO_LANGSMITH` (required by the Pi extension).
3. Project metadata file: [`.pi/langsmith.json`](.pi/langsmith.json) (env `LANGSMITH_PROJECT` wins)
4. `./run.sh` then inside Pi: `/langsmith-tracing`

UI: [smith.langchain.com](https://smith.langchain.com) → project **reAgent-hackathon**.

## REPL (chat + tools, keeps context)

Text chat with the immuno-risk tools and full conversation memory (no Pi TUI):

```bash
cd harness
./repl.sh
# or: npm run repl
```

Needs `ANTHROPIC_API_KEY` in `../.env` or `harness/.env`.

| Slash | Action |
| --- | --- |
| `/help` | Commands |
| `/clear` | Reset chat context |
| `/history` | Show recent turns |
| `/tools` | List immuno tools |
| `/quit` | Exit |

Example: paste a FASTA, or ask `run MHC I immuno pipeline on the insulin sequence we discussed`.

## Launch (full Pi TUI)

```bash
cd harness
./run.sh
# or:
set -a && source ../.env && set +a && npm run pi
```

Inside Pi (trust the project when prompted, or `/trust`):

- `/langsmith-tracing` — confirm LangSmith export is on
- `/denovo` — start/resume late immuno-risk workflow (reads `TASK.md`)
- `/immuno-tools` — list custom immuno-risk tools
- Custom tools: `list_catalytic_sites`, `structure_features`, `predict_cleavage`, `score_mhc`, `check_tolerance`, `score_immuno_risk`, `run_immuno_pipeline`
- `/skill:reagent` then Proto / Boltz / Paperclip skills as needed
- `/mcp` — Paperclip + Proto status
- `/mcp-auth paperclip` / `/mcp-auth proto-bio` if prompted
- `/reload` after editing settings, skills, or `TASK.md`

Auth: `ANTHROPIC_API_KEY` in `../.env`, or `/login`.

## Before the agent builds

Confirm [`TASK.md`](TASK.md) matches the team brief (late immuno-risk stage). Design doc: [`../docs/IMMUNO_RISK_DESIGN.md`](../docs/IMMUNO_RISK_DESIGN.md).

## Immuno-risk Python backend

From **repo root**:

```bash
uv sync --extra immuno          # MHCflurry + benchling-sdk (+ torch)
# Optional real models:
# mhcflurry-downloads fetch models_class1_presentation

# Offline / first-boot demo without MHCflurry weights:
export IMMUNO_ALLOW_HEURISTIC_MHC=1

uv run python -m re_agent.immuno_risk.cli ensure-fixtures
uv run python -m re_agent.immuno_risk.cli train
uv run python -m re_agent.immuno_risk.cli run --sequence GILGFVFTLAAAAKKKLLLGGGG --sequence-id demo
uv run python scripts/benchmark_immuno_risk.py
uv run pytest tests/test_immuno_risk.py -q

# Benchling (credentials from booth / Discord)
uv run python -m re_agent.immuno_risk.cli benchling-pull --dry-run
# uv run python -m re_agent.immuno_risk.cli benchling-publish --run-dir results/immuno_risk/<run-id> --dry-run
```

Optional licensed comparators: set `NETMHCPAN_BIN` / `NETMHCIIPAN_BIN` to local DTU academic binaries (pin 4.2e / 4.3k). Do not automate the public web form.

## Layout

| Path | Purpose |
| --- | --- |
| [`pi-tools.ts`](pi-tools.ts) | Custom immuno-risk + Benchling tools |
| [`tools/`](tools/) | Cleavage / MHC bridge / tolerance / risk / pipeline |
| [`../src/re_agent/immuno_risk/`](../src/re_agent/immuno_risk/) | MHCflurry, IEDB head, Atlas, aggregation, Benchling |
| [`.pi/extensions/immuno-risk-tools.ts`](.pi/extensions/immuno-risk-tools.ts) | Registers tools into Pi |
| [`.pi/langsmith.json`](.pi/langsmith.json) | LangSmith project + metadata (no API key) |
| [`.pi/settings.json`](.pi/settings.json) | Provider, packages, skills, sessions |
| [`.pi/prompts/denovo.md`](.pi/prompts/denovo.md) | `/denovo` prompt template |
| [`.mcp.json`](.mcp.json) | Paperclip + Proto MCP |
| [`AGENTS.md`](AGENTS.md) | Harness instructions (parent `AGENTS.md` also loads) |
| `../skills/` | Sundial skills (reagent, paperclip, census, proto, boltz) |

Custom tools: `list_catalytic_sites`, `structure_features`, `predict_cleavage`, `score_mhc`, `check_tolerance`, `score_immuno_risk`, `run_immuno_pipeline`, `benchling_pull_candidates`, `benchling_publish_run`

`node_modules/`, `.pi/npm/`, and `.pi/sessions/` stay local (gitignored).

## Docs

- Design: [`../docs/IMMUNO_RISK_DESIGN.md`](../docs/IMMUNO_RISK_DESIGN.md)
- [pi.dev quickstart](https://pi.dev/docs/latest/quickstart)
- [pi-mcp-adapter](https://pi.dev/packages/pi-mcp-adapter)
- Repo onboarding: [`../START_HERE.md`](../START_HERE.md)
