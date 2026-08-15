#!/usr/bin/env bash
# Launch the project-local Pi harness with parent .env loaded.
set -euo pipefail
cd "$(dirname "$0")"

if [[ ! -d node_modules/@earendil-works/pi-coding-agent ]]; then
  echo "Missing deps. Run: npm install" >&2
  exit 1
fi

load_env() {
  local f="$1"
  if [[ -f "$f" ]]; then
    set -a
    # shellcheck disable=SC1090
    source "$f"
    set +a
  fi
}

load_env ../.env
load_env .env

if [[ -z "${ANTHROPIC_API_KEY:-}" ]]; then
  echo "Warning: ANTHROPIC_API_KEY not set — use /login inside Pi or add it to ../.env" >&2
fi

# Pi's LangSmith extension reads TRACE_TO_LANGSMITH (not LANGSMITH_TRACING).
# Map the standard LangChain/LangSmith SDK flag if that's what you set.
if [[ -z "${TRACE_TO_LANGSMITH:-}" && -n "${LANGSMITH_TRACING:-}" ]]; then
  export TRACE_TO_LANGSMITH="$LANGSMITH_TRACING"
fi

# Fallbacks the extension already supports:
#   LANGSMITH_PI_API_KEY  ← LANGSMITH_API_KEY
#   LANGSMITH_PI_ENDPOINT ← LANGSMITH_ENDPOINT
#   LANGSMITH_PI_PROJECT  ← LANGSMITH_PROJECT

if [[ -n "${LANGSMITH_API_KEY:-}${LANGSMITH_PI_API_KEY:-}" ]]; then
  export TRACE_TO_LANGSMITH="${TRACE_TO_LANGSMITH:-true}"
  local_project="${LANGSMITH_PI_PROJECT:-${LANGSMITH_PROJECT:-reAgent-hackathon}}"
  echo "LangSmith: on → project=${local_project} (inside Pi: /langsmith-tracing)"
else
  echo "LangSmith: set LANGSMITH_API_KEY in ../.env to enable traces"
fi

exec npm run pi -- "$@"
