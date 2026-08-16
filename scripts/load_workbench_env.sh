#!/usr/bin/env bash

load_workbench_env() {
  local root="$1"
  local fallback="$root/harness/.env"

  if [[ -f "$root/.env" ]]; then
    set -a
    # shellcheck disable=SC1091
    source "$root/.env"
    set +a
  fi

  if [[ -f "$fallback" ]]; then
    local key
    local value
    for key in \
      ANTHROPIC_API_KEY \
      LANGSMITH_API_KEY \
      PROTO_API_KEY \
      PAPERCLIP_API_KEY \
      PAPERCLIP_MCP_BEARER_TOKEN \
      HF_TOKEN \
      MODAL_TOKEN_ID \
      MODAL_TOKEN_SECRET; do
      if [[ -z "${!key:-}" ]]; then
        value="$(
          set -a
          # shellcheck disable=SC1090
          source "$fallback"
          printf '%s' "${!key:-}"
        )"
        if [[ -n "$value" ]]; then
          printf -v "$key" '%s' "$value"
          export "$key"
        fi
      fi
    done
  fi

  if [[ -x "$HOME/.local/bin/paperclip" ]]; then
    export PATH="$HOME/.local/bin:$HOME/.paperclip/bin:$PATH"
  fi

  local sibling_proto_repo
  sibling_proto_repo="$(dirname "$root")/proto-tools"
  if [[ -z "${PROTO_TOOLS_REPO:-}" && -f "$sibling_proto_repo/pyproject.toml" ]]; then
    export PROTO_TOOLS_REPO="$sibling_proto_repo"
  fi
}
