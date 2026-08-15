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
echo "Python env is ready. Finish the rest yourself (opens a browser):"
echo "  Full guide: docs/SETUP.md"
echo
echo "  1. curl -fsSL https://paperclip.gxl.ai/install.sh | bash"
echo "  2. paperclip login && paperclip install"
echo "  3. Put YOUR keys in .env (ANTHROPIC_API_KEY, PROTO_API_KEY)"
echo "     Do not copy a teammate's .env"
echo "  4. Cursor: Cmd+Shift+P → Tools & MCPs → enable paperclip → authenticate"
echo "  5. uv run python scripts/check_setup.py"
echo
echo "Need push access? Ask Vikas (@Kakar13) to add you as a collaborator."
