#!/usr/bin/env bash
# Interactive immuno-risk REPL (chat context + tools).
set -euo pipefail
cd "$(dirname "$0")"

if [[ ! -d node_modules/@anthropic-ai/sdk ]]; then
  echo "Installing deps…" >&2
  npm install
fi

if [[ -f ../.env ]]; then
  set -a
  # shellcheck disable=SC1091
  source ../.env
  set +a
fi
if [[ -f .env ]]; then
  set -a
  # shellcheck disable=SC1091
  source .env
  set +a
fi

if [[ -z "${ANTHROPIC_API_KEY:-}" ]]; then
  echo "Missing ANTHROPIC_API_KEY in ../.env or harness/.env" >&2
  exit 1
fi

exec npx tsx repl.ts "$@"
