# re:AGENT — Pi project instructions

Hackathon repo for scientific agents. Weekend of Aug 15–16, 2026. Submission Sunday 10:45 AM.

We use the [Pi](https://pi.dev/) coding agent harness (`pi` in this directory). See [docs/SETUP.md](docs/SETUP.md).

## Tools to use

- **Paperclip** — literature, FDA, trials, UniProt/PDB/ChEMBL. Prefer the CLI (`paperclip search|map|grep|sql|cat`). Pi has no built-in MCP; drive Paperclip via `bash`. Do not download entire papers by hand.
- **Proto** — two paths: (1) hosted SDK `proto_client` with `PROTO_API_KEY`; (2) local `proto-tools` / `proto-language` via `uv sync --extra proto`, run on Modal with `device="modal"` after `uv run modal setup` and deploy. Do not invent Modal tokens.
- **Claude / Anthropic** — `ANTHROPIC_API_KEY` in `.env` (Pi can also `/login` for subscription providers).

Do not invent API keys. If a check fails, tell the user which booth / Discord channel to get the credential from. Full partner list: `SPONSORS.md`.

New teammates: open `START_HERE.md`, then follow `docs/SETUP.md`. Each person uses their own `.env` and Paperclip login. Never copy keys between laptops.

## Working rules

- Keep secrets in `.env` only. Never commit `.env`.
- Write durable code under `src/re_agent/`. Scratch work goes in `notebooks/` or `scripts/`.
- Large artifacts (FASTA, PDB, paper dumps) go in `data/` or `results/` — those dirs are gitignored except `.gitkeep`.
- Prefer inspectable structured outputs (JSON/CSV/Markdown tables) over chatty prose. Judges need to see the reasoning.
- Paperclip workflow that works: `search` → `map` over the result id → `reduce` / synthesize.
- After setup changes, run `uv run python scripts/check_setup.py`.
- Use `!command` in Pi for shell; keep commits intentional (ask before committing unless the user asked).

## Tracks (pick one and stay on it)

- A: end-to-end scientific agent with tools + inspectable trace
- B: dataset / meta-analysis with a finding no single paper shows
- C: biological design to a spec, plus evidence it could hold up
