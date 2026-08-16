#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "==> Creating .env if missing"
if [[ ! -f .env ]]; then
  cp .env.example .env
  echo "    wrote .env — fill in API keys from Discord / lightning talks"
fi

echo "==> Syncing Python 3.12 env with uv (base deps)"
uv python pin 3.12
uv sync

echo "==> K-Dense scientific-agent-skills (162 skills, optional submodule)"
./scripts/install_scientific_skills.sh

echo
echo "Python env is ready. Finish the rest yourself (opens a browser):"
echo "  Full guide: docs/SETUP.md · Pi harness: harness/README.md"
echo
echo "  1. Put YOUR keys in .env (ANTHROPIC_API_KEY, PROTO_API_KEY)"
echo "     Do not copy a teammate's .env"
echo "  2. cd harness && npm install && npx pi install -l --approve npm:pi-mcp-adapter"
echo "  3. Fill harness/TASK.md, then ./run.sh → /denovo"
echo "  4. curl -fsSL https://paperclip.gxl.ai/install.sh | bash"
echo "  5. paperclip login && paperclip install"
echo "  6. uv run python scripts/check_setup.py"
echo
echo "  MCP: harness/.mcp.json · Adapter: https://pi.dev/packages/pi-mcp-adapter"
echo "  Node ≥ 22.19 recommended for Pi 0.84+"
echo
echo "Proto + Modal (when you need fold/design/optimize — opt-in):"
echo "  uv sync --extra proto"
echo "  uv run modal setup"
echo "  uv run modal environment create proto-env"
echo "  Details: docs/SETUP.md §7 + https://proto.evodesign.org/docs/tools/modal-integration"
echo
echo "Need push access? Ask Vikas (@Kakar13) to add you as a collaborator."
