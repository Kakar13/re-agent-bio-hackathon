#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# shellcheck disable=SC1091
source "$ROOT/scripts/load_workbench_env.sh"
load_workbench_env "$ROOT"

export LANGSMITH_TRACING="${LANGSMITH_TRACING:-true}"
export LANGSMITH_PROJECT="${LANGSMITH_PROJECT:-re-agent-scientific-workbench}"
export NEXT_PUBLIC_LANGGRAPH_API_URL="${NEXT_PUBLIC_LANGGRAPH_API_URL:-http://localhost:2024}"

if [[ ! -d workbench/node_modules ]]; then
  npm install --prefix workbench
fi

npm run build --prefix workbench
uv run langgraph dev --host 127.0.0.1 --port 2024 --no-reload --no-browser &
GRAPH_PID=$!
npm run start --prefix workbench -- --hostname 127.0.0.1 --port 3000 &
WEB_PID=$!

terminate_tree() {
  local parent="$1"
  local child
  while IFS= read -r child; do
    [[ -n "$child" ]] && terminate_tree "$child"
  done < <(pgrep -P "$parent" 2>/dev/null || true)
  kill "$parent" 2>/dev/null || true
}

cleanup() {
  terminate_tree "$GRAPH_PID"
  terminate_tree "$WEB_PID"
}
trap cleanup EXIT INT TERM

echo "Scientific workbench: http://localhost:3000"
echo "LangGraph API:       http://localhost:2024"
echo "LangSmith project:   $LANGSMITH_PROJECT"
wait
