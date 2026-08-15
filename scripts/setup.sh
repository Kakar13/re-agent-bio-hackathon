#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "==> Creating .env if missing"
if [[ ! -f .env ]]; then
  cp .env.example .env
  echo "    wrote .env — fill in API keys from Discord / lightning talks"
fi

echo "==> Syncing Python 3.12 env with uv"
uv python pin 3.12
uv sync

echo
echo "Next (run these in your own terminal — they open a browser):"
echo "  1. curl -fsSL https://paperclip.gxl.ai/install.sh | bash"
echo "  2. paperclip login && paperclip install"
echo "  3. Fill PROTO_API_KEY and ANTHROPIC_API_KEY in .env"
echo "  4. uv run python scripts/check_setup.py"
echo
echo "Cursor MCP: Paperclip + Proto are already in .cursor/mcp.json"
echo "  Cmd+Shift+P → Tools & MCPs → enable paperclip, then authenticate"
