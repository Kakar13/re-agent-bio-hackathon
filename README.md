# re:AGENT

Two-day build weekend for scientific agents: better datasets, sharper tools, and evaluation you can trust.

**When:** Saturday Aug 15 – Sunday Aug 16, 2026  
**Where:** 2 Marina Boulevard, Building C, 3rd floor (left off the elevator)  
**Event:** [luma.com/g6org075](https://luma.com/g6org075)  
**Deadline:** Sunday 10:45 AM — final submission

You keep what you build. Open-source it, keep going, or spin it out.

**Teammates:** open **[START_HERE.md](START_HERE.md)** first. Then [docs/SETUP.md](docs/SETUP.md) for the full walkthrough. Do not reuse someone else's `.env`.

**Agent skills (Sundial):** [`skills/`](skills/) — start with [`skills/reagent/SKILL.md`](skills/reagent/SKILL.md), then Paperclip / Census / Proto / Boltz.

**Sponsors & co-hosts (use them):** [SPONSORS.md](SPONSORS.md)

## Tracks

| Track | What you demo |
| --- | --- |
| **A — Build an AI Scientist** | An agent that runs a scientific or drug-dev workflow end to end: gather evidence, use tools/DBs, generate and test hypotheses, produce inspectable structured output. |
| **B — Build a Dataset or Meta-Analysis** | Queries across thousands of papers, then a pattern no single paper could show. |
| **C — Build the Biological Design** | A pipeline that designs biology to a spec — sequence or system that did not exist before you built it. |

Or bring your own project.

## Weekend tools

Credits and access come from lightning talks + Discord. Do not commit keys.

| Tool | What it is | Get started |
| --- | --- | --- |
| **Pi** ([pi.dev](https://pi.dev/)) | Primary coding-agent harness (TUI, skills, AGENTS.md) | `curl -fsSL https://pi.dev/install.sh \| sh` |
| **Paperclip** (GXL) | Agent-native search over 11M+ papers, FDA docs, trials, UniProt / PDB / ChEMBL | `curl -fsSL https://paperclip.gxl.ai/install.sh \| bash` · MCP: `https://paperclip.gxl.ai/mcp` |
| **Proto** (Arc Institute) | proto-tools + proto-language (AlphaFold, Boltz, Evo2, ESMC, …) | `uv sync --extra proto` · MCP: `https://mcp.evodesign.org/mcp` |
| **Claude** (Anthropic) | API credits for Pi / Claude Code | `ANTHROPIC_API_KEY` in `.env` or Pi `/login` |
| **Modal** | Remote GPU for Proto tools ($100 re:AGENT credits) | `uv run modal setup` · [docs](https://proto.evodesign.org/docs/tools/modal-integration) |
| **Benchling, BenchFlow, LatchBio, Boltz, Strand** | Lab, eval, compute, structure, sequence partners | on-site + Discord |

## Quick start

Full teammate walkthrough (clone → keys → Paperclip → Cursor MCP → verify): **[docs/SETUP.md](docs/SETUP.md)**

```bash
git clone https://github.com/Kakar13/re-agent-bio-hackathon.git
cd re-agent-bio-hackathon
./scripts/setup.sh
```

Then, in **your** terminal (these open a browser):

```bash
# Paperclip CLI + agent skill — each person logs in as themselves
curl -fsSL https://paperclip.gxl.ai/install.sh | bash
paperclip login
paperclip install

# Pi — primary coding agent (https://pi.dev/)
curl -fsSL https://pi.dev/install.sh | sh
pi install npm:pi-mcp-adapter   # Paperclip + Proto via .mcp.json

# Keys from Discord / lightning talks — your keys, not a teammate's
# setup.sh already copied .env.example → .env
# edit .env: ANTHROPIC_API_KEY, PROTO_API_KEY

uv run python scripts/check_setup.py
set -a && source .env && set +a && pi
# /mcp
```

In Pi, project rules load from [`AGENTS.md`](AGENTS.md). MCP: [pi-mcp-adapter](https://pi.dev/packages/pi-mcp-adapter) + [`.mcp.json`](.mcp.json). Cursor: optional **Cmd+Shift+P → Tools & MCPs**.

Need write access? Ask Vikas ([@Kakar13](https://github.com/Kakar13)) to add you as a collaborator.

### Smoke the tools

```bash
paperclip search "CRISPR base editing efficiency" -n 5
paperclip config

# Proto stack (proto-tools + proto-language + modal) — opt-in
uv sync --extra proto
uv run modal setup                    # browser login; $100 re:AGENT credits day-of
# Hosted SDK (needs PROTO_API_KEY in .env):
uv run python -c "from proto_client import ProtoClient; print(ProtoClient())"
```

Full Proto + Modal steps: [docs/SETUP.md](docs/SETUP.md#7-proto--modal-when-you-design--fold--optimize-sequences).

## Repo layout

```
src/re_agent/     shared config + setup checks
scripts/          setup.sh, check_setup.py
notebooks/        exploration
data/raw          inputs you collect (gitignored)
data/processed    cleaned tables (gitignored)
results/          demo artifacts (gitignored)
AGENTS.md         Pi project instructions
.mcp.json         Shared MCP (Paperclip + Proto) for Pi / tools
.cursor/mcp.json  Same servers for Cursor
```

## Useful Paperclip commands

```bash
paperclip search "GLP-1 receptor agonist efficacy" -n 20
paperclip search -s fda "pembrolizumab"
paperclip search -s trials/us "phase 3 NSCLC"
paperclip search -s uniprot "p53"
paperclip map --from s_RESULT_ID "What were the primary endpoints?"
paperclip sql "SELECT title, doi FROM documents WHERE authors ILIKE '%Doudna%' LIMIT 5"
```

Every paper is a directory: `/papers/<id>/meta.json`, `content.lines`, `sections/`, `figures/`.

## Schedule (local)

**Sat:** 8:30 check-in · 9:15 welcome · 9:35 lightning talks · 10:25 ideation / teams · 12:10 lunch · 1:00–9:45 build · 9:45 overnight checkpoint  
**Sun:** 9:00 final build · **10:45 submit** · 11:30 panel · 12:30 demos · 2:10 awards

## Team

- Vikas Kakar ([@Kakar13](https://github.com/Kakar13))

Add names here after formation. New teammates: [docs/SETUP.md](docs/SETUP.md).

## Links

- Sponsors & co-hosts: [SPONSORS.md](SPONSORS.md)
- Paperclip docs: [paperclip.gxl.ai/docs](https://paperclip.gxl.ai/docs)
- Pi: [pi.dev](https://pi.dev/) · [quickstart](https://pi.dev/docs/latest/quickstart)
- Proto: [proto.evodesign.org](https://proto.evodesign.org)
- Proto tools: [github.com/evo-design/proto-tools](https://github.com/evo-design/proto-tools)
- Proto language: [github.com/evo-design/proto-language](https://github.com/evo-design/proto-language)
- Proto Modal: [docs/tools/modal-integration](https://proto.evodesign.org/docs/tools/modal-integration)
- Claude Code: [code.claude.com/docs](https://code.claude.com/docs/en/quickstart)
