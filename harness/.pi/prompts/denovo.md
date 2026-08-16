---
description: Start or continue the de novo design and MHC-II immunogenicity workflow
argument-hint: "[focus note]"
---
We are running the inspectable Track A workflow for de novo binder design and immunogenicity screening.

Pipeline: **design specification → RFdiffusion3 → ProteinMPNN → AlphaFold2 validation → response-model adapters → NetMHCIIpan EL and BA → processing/tolerance evidence → transparent candidate rank**.

1. Read `TASK.md` in this harness folder and follow it.
2. Read `../skills/reagent/SKILL.md`, then `../skills/paperclip/SKILL.md`, and Proto/Boltz skills only as needed for structure or complex checks.
3. Use the registered tools: `immuno_architecture_status` to inspect readiness and `run_immuno_pipeline` for the real Python adapter pipeline. Never call the retired placeholder tools.
4. Work from the **repo root** (`..`): durable code in `../src/re_agent/`, artifacts in `../results/immuno_risk/`.
5. Keep response propensity, NetMHCIIpan EL, NetMHCIIpan BA, processing, and tolerance separate. If the teammate response artifact is unavailable, preserve that status and withhold fusion.
6. Every claim must be cited through Paperclip or measured through versioned tool output. De novo candidates have no direct immune-response ground truth, so report uncertainty and claim boundaries.
7. Use LangSmith for graph/model/tool observability when configured, and retain local run manifests for reproducibility.
8. Prefer MCP (`/mcp`) or CLI; do not download whole papers by hand.

Focus for this turn: ${@:-Inspect architecture status, then run NetMHCIIpan EL and BA on a short natural control with the response adapter explicitly unavailable.}
