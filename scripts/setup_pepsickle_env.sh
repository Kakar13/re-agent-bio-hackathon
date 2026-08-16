#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Pepsickle (proteasomal cleavage teacher, used by src/re_agent/e2e_pls/label.py)
# pins torch==1.13.1 / scikit-learn==1.2.0, which conflicts with this repo's
# main env (torch>=2.2 for ESM3/Track 2). It gets its own isolated venv here,
# invoked as a CLI subprocess -- see label.py's run_pepsickle_on_proteins().

echo "==> Creating isolated venv at .tools/pepsickle"
uv venv .tools/pepsickle --python 3.11

echo "==> Installing pepsickle"
uv pip install --python .tools/pepsickle/bin/python pepsickle

echo
echo "Done. Verify with:"
echo "  .tools/pepsickle/bin/pepsickle -s VSGLEQLESIINFEKLTEWTSSNV"
