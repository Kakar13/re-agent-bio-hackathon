# Start here

Welcome to **re:AGENT**. This repo is the shared project for the weekend.

**Repo:** https://github.com/Kakar13/re-agent-bio-hackathon  
**Deadline:** Sunday 10:45 AM · **Event:** [luma.com/g6org075](https://luma.com/g6org075)

Each person sets up on their own laptop. **Do not copy someone else's `.env`.**

---

## 1. Clone

```bash
git clone https://github.com/Kakar13/re-agent-bio-hackathon.git
cd re-agent-bio-hackathon
```

Need write access? Ask Vikas ([@Kakar13](https://github.com/Kakar13)) to add you as a collaborator.

## 2. Base Python env

```bash
# Need uv? https://docs.astral.sh/uv/
curl -LsSf https://astral.sh/uv/install.sh | sh

chmod +x scripts/setup.sh scripts/check_setup.py
./scripts/setup.sh
```

This creates `.env` from `.env.example` and installs base deps (Python 3.12).

## 3. Keys (your own)

Edit `.env`:

| Key | From |
| --- | --- |
| `ANTHROPIC_API_KEY` | Anthropic booth / Discord / lightning talks |
| `PROTO_API_KEY` | Proto workspace (Arc Institute) |

Optional later: `HF_TOKEN`, Modal via `uv run modal setup` (not required for first boot).

Full partner list: [SPONSORS.md](SPONSORS.md)

## 4. Paperclip

```bash
curl -fsSL https://paperclip.gxl.ai/install.sh | bash
paperclip login
paperclip install
paperclip config
```

## 5. Pi (primary coding agent)

```bash
curl -fsSL https://pi.dev/install.sh | sh

set -a && source .env && set +a
pi
```

Pi loads [AGENTS.md](AGENTS.md). Auth: API key in `.env`, or `/login` inside Pi.  
Docs: [pi.dev](https://pi.dev/) · [quickstart](https://pi.dev/docs/latest/quickstart)

**Note:** Pi has no built-in MCP — use `!paperclip …` and `uv run …` from the shell.

## 6. Verify

```bash
uv run python scripts/check_setup.py
```

`pi` + `paperclip` should be OK before you start building. Proto/Modal are opt-in (§ below).

## 7. Proto + Modal (when you design / fold)

```bash
uv sync --extra proto
uv run modal setup                          # $100 re:AGENT credits, day-of
uv run modal environment create proto-env
```

Details: [docs/SETUP.md](docs/SETUP.md#7-proto--modal-when-you-design--fold--optimize-sequences)

## 8. Branch and build

```bash
git checkout main && git pull
git checkout -b yourname/short-topic
# ... work ...
git add -p
git commit -m "Why this change exists."
git push -u origin HEAD
```

---

## Read next

| Doc | Purpose |
| --- | --- |
| [docs/SETUP.md](docs/SETUP.md) | Full teammate setup (troubleshooting, Cursor MCP, etc.) |
| [SPONSORS.md](SPONSORS.md) | Co-hosts & sponsors — use them in the demo |
| [AGENTS.md](AGENTS.md) | Rules for Pi / agents |
| [README.md](README.md) | Tracks, schedule, Paperclip cheatsheet |

## Tracks (pick one)

- **A** — AI scientist (end-to-end agent + inspectable trace)
- **B** — Dataset / meta-analysis (finding no single paper shows)
- **C** — Biological design (sequence/system to a spec)

Stuck? Discord + partner booths all weekend.
