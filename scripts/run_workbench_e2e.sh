#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# shellcheck disable=SC1091
source "$ROOT/scripts/load_workbench_env.sh"
load_workbench_env "$ROOT"

GRAPH_PORT="${PLAYWRIGHT_GRAPH_PORT:-2124}"
WEB_PORT="${PLAYWRIGHT_WEB_PORT:-3100}"
if [[ -z "${LANGSMITH_API_KEY:-}" ]]; then
  export LANGSMITH_TRACING=false
else
  export LANGSMITH_TRACING="${LANGSMITH_TRACING:-true}"
fi
export LANGSMITH_PROJECT="${LANGSMITH_PROJECT:-re-agent-playwright}"
export NEXT_PUBLIC_LANGGRAPH_API_URL="http://127.0.0.1:${GRAPH_PORT}"
export RE_AGENT_FORCE_KEYLESS=true

cleanup() {
  local listeners
  listeners="$(lsof -tiTCP:"$WEB_PORT" -sTCP:LISTEN 2>/dev/null || true)"
  [[ -n "$listeners" ]] && kill $listeners 2>/dev/null || true
  listeners="$(lsof -tiTCP:"$GRAPH_PORT" -sTCP:LISTEN 2>/dev/null || true)"
  [[ -n "$listeners" ]] && kill $listeners 2>/dev/null || true
}
trap cleanup EXIT INT TERM

npm run build --prefix workbench
uv run langgraph dev \
  --host 127.0.0.1 \
  --port "$GRAPH_PORT" \
  --no-reload \
  --no-browser &
npm run start --prefix workbench -- \
  --hostname 127.0.0.1 \
  --port "$WEB_PORT" &

wait
