# re:AGENT project

Hackathon repo for scientific agents. Weekend of Aug 15–16, 2026. Submission Sunday 10:45 AM.

## Tools to use

- **Paperclip** — literature, FDA, trials, UniProt/PDB/ChEMBL. Prefer the CLI (`paperclip search|map|grep|sql|cat`) or the project MCP server `paperclip`. Mention `/paperclip` so the skill loads.
- **Proto** — sequence design and 80+ bio tools. Hosted MCP `proto-bio` at `https://mcp.evodesign.org/mcp` (needs `PROTO_API_KEY`). Python: `from proto_client import ProtoClient`.
- **Claude API** — `ANTHROPIC_API_KEY` in `.env`.

Do not invent API keys. If a check fails, tell the user which booth / Discord channel to get the credential from.

New teammates follow `docs/SETUP.md`. Each person uses their own `.env` and Paperclip login. Never copy keys between laptops.

## Working rules

- Keep secrets in `.env` only. Never commit `.env`.
- Write durable code under `src/re_agent/`. Scratch work goes in `notebooks/` or `scripts/`.
- Large artifacts (FASTA, PDB, paper dumps) go in `data/` or `results/` — those dirs are gitignored except `.gitkeep`.
- Prefer inspectable structured outputs (JSON/CSV/Markdown tables) over chatty prose. Judges need to see the reasoning.
- Paperclip workflow that works: `search` → `map` over the result id → `reduce` / synthesize. Do not download entire papers by hand.

## Tracks (pick one and stay on it)

- A: end-to-end scientific agent with tools + inspectable trace
- B: dataset / meta-analysis with a finding no single paper shows
- C: biological design to a spec, plus evidence it could hold up
