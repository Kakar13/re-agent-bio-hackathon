# Teammate setup

You should be ready in about 10–15 minutes. Each person does this on their own laptop. **Do not copy someone else's `.env`.** Keys are per-person from Discord / lightning talks.

**Fast path:** [START_HERE.md](../START_HERE.md) · **Repo:** https://github.com/Kakar13/re-agent-bio-hackathon

## 0. What you need

| Need | Notes |
| --- | --- |
| macOS, Linux, or WSL | Proto's local stack is macOS/Linux. Windows: use WSL. |
| GitHub account | Repo is public — clone is enough. Ask Vikas (`@Kakar13`) for write access. |
| Cursor, Claude Code, or **Pi** | **Pi is our primary harness** ([pi.dev](https://pi.dev/)). MCP via [pi-mcp-adapter](https://pi.dev/packages/pi-mcp-adapter) + `.mcp.json`. |
| Event Discord | Paperclip login, Proto key, Claude credits, Modal credits ($100 re:AGENT). |

Install these if you do not already have them:

```bash
# GitHub CLI (optional, for PRs)
brew install gh          # macOS
# or: https://cli.github.com

# uv — Python package manager we use
curl -LsSf https://astral.sh/uv/install.sh | sh

# Pi coding agent (primary harness)
curl -fsSL https://pi.dev/install.sh | sh
# or: npm install -g --ignore-scripts @earendil-works/pi-coding-agent
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

That copies `.env.example` → `.env` (if missing) and runs `uv sync` on Python 3.12 (base deps only — Proto is opt-in in §7).

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
| `PROTO_API_KEY` | Proto workspace (Arc Institute) | Hosted Proto MCP + `proto-client` SDK |
| `HF_TOKEN` | huggingface.co → Settings → Access Tokens | Gated models only (ESM3, AlphaFold3, AlphaGenome) |
| `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET` | Optional fallback | Prefer `uv run modal setup` (writes `~/.modal.toml`) |

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

## 5. Pi coding agent (primary harness)

We use [Pi](https://pi.dev/) — a minimal agent harness. Docs: [quickstart](https://pi.dev/docs/latest/quickstart).

```bash
# macOS / Linux / WSL
curl -fsSL https://pi.dev/install.sh | sh

# or npm
npm install -g --ignore-scripts @earendil-works/pi-coding-agent

pi --version   # or just: which pi
```

### MCP for Pi (Paperclip + Proto)

Pi core does not bake MCP in; install the adapter:

```bash
pi install npm:pi-mcp-adapter
```

Docs: [pi.dev/packages/pi-mcp-adapter](https://pi.dev/packages/pi-mcp-adapter). Restart Pi after install.

This repo ships shared MCP config in **`.mcp.json`** (Paperclip + Proto). The adapter loads it automatically (lazy servers — connect on first use).

```bash
cd /path/to/re-agent-bio-hackathon
set -a && source .env && set +a   # PROTO_API_KEY for proto-bio
pi
# /mcp
# /mcp-auth paperclip    # if prompted
# /mcp-auth proto-bio
```

Already had Cursor-only MCP? Run `/mcp setup` inside Pi to import host configs, or just use `.mcp.json` (preferred shared file).

### Start Pi

```bash
cd /path/to/re-agent-bio-hackathon
set -a && source .env && set +a
pi
```

Authenticate:

- **API key:** `ANTHROPIC_API_KEY` in `.env` / exported shell (event Claude credits), then start `pi`
- **Subscription:** inside Pi run `/login` (Claude Pro/Max, ChatGPT, Copilot, etc.)

Pi loads project instructions from **`AGENTS.md`** (and `CLAUDE.md`) at startup. After editing those files, run `/reload` or restart Pi.

Useful Pi habits:

```bash
pi                              # interactive TUI in this repo
pi -p "Summarize this repo"     # one-shot / scripts
pi -c                           # continue last session
!paperclip search "…" -n 5      # CLI fallback
```

Switch models with `/model` or `Ctrl+L`. Docs: [pi.dev](https://pi.dev/).

### Claude Code (optional)

Still fine if you prefer it; not required when using Pi.

```bash
curl -fsSL https://claude.ai/install.sh | bash
claude --version
```

## 6. Cursor MCP (optional — same servers)

Use this if you work in **Cursor**. Pi users get Paperclip + Proto via `.mcp.json` + `pi-mcp-adapter`.

This repo also ships `.cursor/mcp.json` with the same servers.

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

## 7. Proto + Modal (when you design / fold / optimize sequences)

There are **two** Proto paths. You can use either or both.

| Path | What it is | Auth |
| --- | --- | --- |
| **Hosted `proto-client` / Cursor MCP** | Agent or SDK calls tools over the network | `PROTO_API_KEY` in `.env` (+ Cursor MCP `proto-bio` if in Cursor) |
| **Local proto-tools + Modal** | Run AlphaFold, Boltz, Evo2, ESMC, … on Modal GPUs ($100 re:AGENT credits) | `uv run modal setup` → `~/.modal.toml` |

### Install the official stack (opt-in)

Same as the event handout (`pip install git+https://github.com/evo-design/proto-tools.git` and proto-language), via our uv extra:

```bash
uv sync --extra proto
```

That installs:

- **proto-tools** — ready-to-run models (AlphaFold, Boltz, Evo2, ESMC, Protenix, OpenDDE, AlphaGenome, and more)
- **proto-language** — design loops with segments, constructs, generators, constraints, optimizers ([docs](https://proto.evodesign.org))
- **modal** — remote GPU compute
- **proto-client** — hosted SDK (uses `PROTO_API_KEY`)

Git install is enough to *run* design loops. Only clone proto-language and `pip install -e .` if you need to write custom generators/constraints later.

### Modal compute ($100 re:AGENT credits)

Credits are provided day-of. Authenticate in **your** terminal (opens a browser):

```bash
uv sync --extra proto
uv run modal setup
uv run modal environment create proto-env
uv run proto-tools deploy --list
```

Deploy a tool when you need it, then call it with `device="modal"`:

```bash
uv run proto-tools deploy --apps esmc --env proto-env
```

```python
from proto_tools import run_esmc_embeddings, ESMCEmbeddingsInput, ESMCEmbeddingsConfig

output = run_esmc_embeddings(
    ESMCEmbeddingsInput(sequences=["MKTAYLLIGLLAIAAFSPQVLA"]),
    ESMCEmbeddingsConfig(device="modal"),
)
print(len(output.results[0].mean_embedding))
```

Full walkthrough: [proto.evodesign.org/docs/tools/modal-integration](https://proto.evodesign.org/docs/tools/modal-integration)

### Hosted SDK smoke (optional)

```bash
# PROTO_API_KEY must be set in .env
uv run python -c "from proto_client import ProtoClient; print(ProtoClient())"
```

Confirm credits via MCP: *“Check my Proto workspace and remaining credits.”*

## 8. Confirm you are done

```bash
uv run python scripts/check_setup.py
```

Expected:

| Check | OK means |
| --- | --- |
| `python3` / `uv` | local env works |
| `pi` | Pi coding agent on PATH (primary harness) |
| `claude` | Claude Code on PATH (optional) |
| `paperclip` | CLI installed and logged in |
| `proto-client` | after `uv sync --extra proto` (hosted SDK) |
| `proto-tools` / `proto-language` | after `uv sync --extra proto` |
| `modal` | after `uv sync --extra proto` (CLI on PATH via uv) |
| `ANTHROPIC_API_KEY` | present in `.env` (Pi API auth + Claude) |
| `PROTO_API_KEY` | present in `.env` (hosted Proto SDK / Cursor MCP) |

`!!` on Claude or Proto packages is fine until you need those tools. **`pi` and Paperclip should be `ok` before you start building.** Modal auth (`modal setup`) is separate from the import check — run it before deploying tools.

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

**`pi: command not found`**  
Install with `curl -fsSL https://pi.dev/install.sh | sh`, then restart the terminal. Confirm with `which pi`. Auth via `set -a && source .env && set +a && pi` or `/login` inside Pi.

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

**`proto_tools` / `proto_language` import fails**  
Run `uv sync --extra proto`. Default `./scripts/setup.sh` only installs base deps.

**Modal auth / deploy fails**  
Run `uv run modal setup` again. Claim $100 re:AGENT credits day-of. See [Modal integration docs](https://proto.evodesign.org/docs/tools/modal-integration).

**Python 3.14 on PATH**  
Ignore it. The repo pins 3.12 via `.python-version`. Always use `uv run ...`, not system `python3`.

**Need help**  
Partner teams are on-site and in Discord all weekend. Also: https://paperclip.gxl.ai/docs and https://proto.evodesign.org
