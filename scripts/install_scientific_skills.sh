#!/usr/bin/env bash
# Fetch K-Dense scientific-agent-skills (162 skills) into vendor/.
# Used by setup and after clone without --recurse-submodules.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

VENDOR="vendor/scientific-agent-skills"
REPO="https://github.com/K-Dense-AI/scientific-agent-skills.git"

if [[ -d "$VENDOR/skills" ]]; then
  echo "==> $VENDOR already present ($(ls "$VENDOR/skills" | wc -l | tr -d ' ') skills)"
  exit 0
fi

if [[ -f .gitmodules ]] && grep -q scientific-agent-skills .gitmodules 2>/dev/null; then
  echo "==> Initializing git submodule: $VENDOR"
  git submodule update --init --depth 1 "$VENDOR"
else
  echo "==> Shallow clone (no submodule metadata): $VENDOR"
  mkdir -p vendor
  git clone --depth 1 "$REPO" "$VENDOR"
fi

mkdir -p .agents/skills
ln -sfn "../../$VENDOR" .agents/skills/scientific-agent-skills

echo "==> Done. Pi loads from harness/.pi/settings.json; Cursor from .agents/skills/"
