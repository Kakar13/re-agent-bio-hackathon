# re:AGENT Scientific Workbench

Artifact-first Next.js client for the shared LangGraph scientific runtime.

## Start

From the repository root:

```bash
cp .env.example .env  # once; add ANTHROPIC_API_KEY and LANGSMITH_API_KEY
./scripts/run_workbench.sh
```

Open <http://localhost:3000>. The script starts:

- LangGraph API on port `2024`
- Next.js workbench on port `3000`
- LangSmith tracing under `LANGSMITH_PROJECT`
- local hashed run manifests under `results/workbench/runs/`

The launchers load the canonical root `.env` first, then fill only missing supported
credentials from `harness/.env`; keys are not copied or printed. See
[`START_HERE.md`](../START_HERE.md) for canonical setup.
The complete design-to-screen contract and credential boundary are documented in
[`docs/ORCHESTRATION.md`](../docs/ORCHESTRATION.md).

## Workbench contract

- Chat is a controller, not the final artifact.
- Tool, reviewer, and artifact events stream from the same graph state.
- NetMHCIIpan EL and BA remain separate evidence lanes.
- The team HLA-A*02:01 MHC-I processing surrogate is optional and excluded from MHC-II fusion.
- The composer defaults to `Chao1 + MHC-II`; reviewers can switch to MHC-II-only screening.
- A missing response artifact visibly withholds the combined rank.
- RFdiffusion3, ProteinMPNN, and AlphaFold2 require an explicit UI approval event.
- Structural metrics fail closed before generated candidates reach screening.
- The 95 historical IL-7Ralpha sequences are a no-GPU preflight, not generation input.
- Mol* renders attached structure artifacts; residue heatmaps render spatial evidence tracks.

Run the frozen team checkpoint on a leakage-audited PDA example:

```bash
uv run python scripts/smoke_team_model_pda.py
```

The generated artifact appears in the workbench run list with the MHC-I window table and
residue tracks. See [`models/chao1/MODEL_CARD.md`](../models/chao1/MODEL_CARD.md) for scope.

## UI and system validation

The suite uses the installed Google Chrome channel, so no duplicate browser
download is required. Run:

```bash
cd workbench
npm run test:e2e
```

On a machine without Google Chrome, install Playwright Chromium once with
`npx playwright install chromium` and remove `channel: "chrome"` from
`playwright.config.ts`.

Playwright launches an isolated production build on port `3100` and a LangGraph
runtime on port `2124`, leaving normal development ports untouched. The suite checks:

- workbench shell and LangGraph health
- streamed architecture artifact and deterministic reviewer
- review-tab annotations
- cached NetMHCIIpan 4.3 EL and BA evidence across the frozen 18-allele panel
- withheld fusion when the response model is unavailable
- Proto campaign planning without GPU execution
- historical-campaign screening preflight without GPU generation
- mobile responsive behavior

Failure traces and screenshots are written under `test-results/`. Open the
HTML report with `npm run test:e2e:report`; use `npm run test:e2e:ui` for interactive
debugging.
