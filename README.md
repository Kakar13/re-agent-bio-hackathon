# re:AGENT

Two-day build weekend for scientific agents: better datasets, sharper tools, and evaluation you can trust.

**When:** Saturday Aug 15 – Sunday Aug 16, 2026  
**Where:** 2 Marina Boulevard, Building C, 3rd floor (left off the elevator)  
**Event:** [luma.com/g6org075](https://luma.com/g6org075)  
**Deadline:** Sunday 10:45 AM — final submission

You keep what you build. Open-source it, keep going, or spin it out.

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
| **Paperclip** (GXL) | Agent-native search over 11M+ papers, FDA docs, trials, UniProt / PDB / ChEMBL | `curl -fsSL https://paperclip.gxl.ai/install.sh \| bash` · MCP: `https://paperclip.gxl.ai/mcp` |
| **Proto** (Arc Institute) | Generative biology language + 80+ tools (fold, design, optimize DNA/RNA/protein) | `pip install proto-client` · MCP: `https://mcp.evodesign.org/mcp` |
| **Claude** (Anthropic) | Claude Code + API credits | already on this machine (`claude --version`) |
| **Modal** | Remote GPU for Proto / Boltz | token from sponsor booth |
| **Benchling, BenchFlow, LatchBio, Boltz, Strand** | Lab, eval, compute, structure, sequence partners | on-site + Discord |

## Quick start

```bash
cd re-agent-bio-hackathon
./scripts/setup.sh
```

Then, in **your** terminal (these open a browser):

```bash
# Paperclip CLI + agent skill
curl -fsSL https://paperclip.gxl.ai/install.sh | bash
paperclip login
paperclip install

# Keys from Discord / lightning talks
cp .env.example .env   # if setup.sh did not already
# edit .env: ANTHROPIC_API_KEY, PROTO_API_KEY

uv run python scripts/check_setup.py
```

In Cursor: **Cmd+Shift+P → Tools & MCPs** → enable `paperclip`, authenticate. Proto uses `PROTO_API_KEY` from the environment via `.cursor/mcp.json`.

### Smoke the tools

```bash
paperclip search "CRISPR base editing efficiency" -n 5
paperclip config

# After PROTO_API_KEY is set
uv pip install proto-client
uv run python -c "from proto_client import ProtoClient; print(ProtoClient().whoami())"
```

## Repo layout

```
src/re_agent/     shared config + setup checks
scripts/          setup.sh, check_setup.py
notebooks/        exploration
data/raw          inputs you collect (gitignored)
data/processed    cleaned tables (gitignored)
results/          demo artifacts (gitignored)
.cursor/mcp.json  Paperclip + Proto MCP for this project
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

Add teammates here after formation.

## Links

- Paperclip docs: [paperclip.gxl.ai/docs](https://paperclip.gxl.ai/docs)
- Proto: [proto.evodesign.org](https://proto.evodesign.org)
- Proto client: [github.com/evo-design/proto-client](https://github.com/evo-design/proto-client)
- Claude Code: [code.claude.com/docs](https://code.claude.com/docs/en/quickstart)
