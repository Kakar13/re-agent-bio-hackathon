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

echo
echo "Python env is ready. Finish the rest yourself (opens a browser):"
echo "  Full guide: docs/SETUP.md"
echo
echo "  1. curl -fsSL https://pi.dev/install.sh | sh"
echo "  2. curl -fsSL https://paperclip.gxl.ai/install.sh | bash"
echo "  3. paperclip login && paperclip install"
echo "  4. Put YOUR keys in .env (ANTHROPIC_API_KEY, PROTO_API_KEY)"
echo "     Do not copy a teammate's .env"
echo "  5. set -a && source .env && set +a && pi"
echo "  6. uv run python scripts/check_setup.py"
echo
echo "  Optional Cursor MCP: Cmd+Shift+P → Tools & MCPs → enable paperclip"
echo
echo "Proto + Modal (when you need fold/design/optimize — opt-in):"
echo "  uv sync --extra proto"
echo "  uv run modal setup"
echo "  uv run modal environment create proto-env"
echo "  Details: docs/SETUP.md §7 + https://proto.evodesign.org/docs/tools/modal-integration"
echo
echo "Need push access? Ask Vikas (@Kakar13) to add you as a collaborator."
