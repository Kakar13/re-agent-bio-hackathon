# Teammate setup

You should be ready in about 10 minutes. Each person does this on their own laptop. **Do not copy someone else's `.env`.** Keys are per-person from Discord / lightning talks.

Repo: https://github.com/Kakar13/re-agent-bio-hackathon

## 0. What you need

| Need | Notes |
| --- | --- |
| macOS, Linux, or WSL | Proto's local stack is macOS/Linux. Windows: use WSL. |
| GitHub account | Repo is public — clone is enough. Ask Vikas (`@Kakar13`) for write access. |
| Cursor or Claude Code | Cursor is enough. Claude Code is optional but useful. |
| Event Discord | Paperclip login, Proto key, Claude credits, Modal tokens. |

Install these if you do not already have them:

```bash
# GitHub CLI (optional, for PRs)
brew install gh          # macOS
# or: https://cli.github.com

# uv — Python package manager we use
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Python 3.12 is pinned in the repo. `uv` will download it. You do not need conda.

## 1. Clone

```bash
git clone https://github.com/Kakar13/re-agent-bio-hackathon.git
cd re-agent-bio-hackathon
```

If you already have write access:

```bash
gh repo clone Kakar13/re-agent-bio-hackathon
cd re-agent-bio-hackathon
```

## 2. Python env

```bash
./scripts/setup.sh
```

That copies `.env.example` → `.env` (if missing) and runs `uv sync` on Python 3.12.

If `setup.sh` is not executable:

```bash
chmod +x scripts/setup.sh scripts/check_setup.py
./scripts/setup.sh
```

## 3. Fill your `.env`

Open `.env` and add **your** keys. Leave a line blank if you do not have that tool yet.

| Variable | Where to get it | Required for |
| --- | --- | --- |
| `ANTHROPIC_API_KEY` | Anthropic booth / Discord / lightning talks | Claude API calls in code |
| `PROTO_API_KEY` | Proto workspace (Arc Institute) | Proto SDK + Cursor MCP `proto-bio` |
| `HF_TOKEN` | huggingface.co → Settings → Access Tokens | Gated models only (ESM3, AlphaFold3, AlphaGenome) |
| `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET` | Modal booth | Remote GPU for Proto / Boltz |

Rules:

- Never commit `.env`. It is gitignored.
- Never paste keys into Slack/Discord/chat with the agent.
- `PYTHONPATH=src` stays in `.env` so `uv run` can import `re_agent`.

## 4. Paperclip (everyone)

Run this in **your** terminal, not through an agent. The installer opens a browser.

```bash
curl -fsSL https://paperclip.gxl.ai/install.sh | bash
paperclip login
paperclip install          # teaches Cursor / Claude Code the Paperclip skill
paperclip config           # Auth should show your email
```

Docs: https://paperclip.gxl.ai/docs

Smoke test:

```bash
paperclip search "CRISPR base editing efficiency" -n 5
```

## 5. Claude Code (optional)

```bash
# macOS / Linux / WSL
curl -fsSL https://claude.ai/install.sh | bash

# Windows PowerShell
irm https://claude.ai/install.ps1 | iex

claude --version
cd /path/to/re-agent-bio-hackathon
claude                      # log in when prompted
```

Docs: https://code.claude.com/docs/en/quickstart

## 6. Cursor MCP (Paperclip + Proto)

This repo already ships `.cursor/mcp.json` with both servers.

1. Restart Cursor after cloning, or reopen the project folder.
2. **Cmd+Shift+P** (macOS) / **Ctrl+Shift+P** (Windows) → **Tools & MCPs**.
3. Enable **paperclip**. Authenticate in the browser when asked.
4. Enable **proto-bio**. It reads `PROTO_API_KEY` from your environment.

If Proto shows 0 tools or "authentication failed":

```bash
# make sure the key is in this shell, then restart Cursor
export PROTO_API_KEY="your-key"
```

On macOS, GUI apps often do not see keys that only live in `.zshrc`. Either:

- put `PROTO_API_KEY` in `.env` **and** launch Cursor from the project terminal: `cursor .`
- or add `export PROTO_API_KEY=...` to `~/.zshrc` and restart Cursor from Spotlight after opening a new terminal.

Hosted endpoints (no local install required for MCP):

- Paperclip: `https://paperclip.gxl.ai/mcp`
- Proto: `https://mcp.evodesign.org/mcp`

## 7. Proto Python SDK (when you start designing sequences)

```bash
uv sync --extra proto
uv run python -c "from proto_client import ProtoClient; print(ProtoClient())"
```

You need `PROTO_API_KEY` set first. Confirm credits with the MCP: *“Check my Proto workspace and remaining credits.”*

## 8. Confirm you are done

```bash
uv run python scripts/check_setup.py
```

Expected:

| Check | OK means |
| --- | --- |
| `python3` / `uv` | local env works |
| `claude` | Claude Code on PATH (optional) |
| `paperclip` | CLI installed and logged in |
| `proto-client` | only after `uv sync --extra proto` |
| `ANTHROPIC_API_KEY` / `PROTO_API_KEY` | present in `.env` |

`!!` on Claude or proto-client is fine until you need those tools. **Paperclip should be `ok` before you start literature work.**

## 9. How we work in this repo

- Branch from `main`: `git checkout -b yourname/short-topic`
- Do not `git add .env` or data dumps. `data/` and `results/` are gitignored except folder placeholders.
- Put durable code in `src/re_agent/`. Scratch work in `notebooks/` or `scripts/`.
- Open a PR or push the branch and tell the team in Discord.
- Ask Vikas to add you as a collaborator if `git push` is denied.

```bash
git checkout main
git pull
git checkout -b yourname/idea
# ... work ...
git add -p
git commit -m "Why this change exists."
git push -u origin HEAD
```

## 10. Common failures

**`uv: command not found`**  
Restart the terminal after installing uv, or run `source $HOME/.local/bin/env`.

**`paperclip: command not found`**  
The installer puts a wrapper in `~/.local/bin`. Add that to PATH:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

**Paperclip auth failed**  
`paperclip login` again. Event credits are tied to the account you used to register.

**Cursor cannot see Proto tools**  
Key missing from the GUI process. Export `PROTO_API_KEY` in `~/.zshrc`, or start Cursor from the repo: `cursor .`

**Python 3.14 on PATH**  
Ignore it. The repo pins 3.12 via `.python-version`. Always use `uv run ...`, not system `python3`.

**Need help**  
Partner teams are on-site and in Discord all weekend. Also: https://paperclip.gxl.ai/docs and https://proto.evodesign.org
