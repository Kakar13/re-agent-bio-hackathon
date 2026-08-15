---
description: Start or continue the late-stage immuno-risk pipeline (cleavage → MHC → risk)
argument-hint: "[focus note]"
---
We are running the **last stage** of the de novo binder project for re:AGENT — not backbone/MPNN design.

Pipeline: **sequence (+ structure) → accessibility features → cleavage vs ~10 catalytic sites → peptide pool → MHC I (priority) / MHC II (optional) → tolerance (HLA ligand / self refs) → risk score**. Optional Boltz/AF check that flagged peptides can engage the recognition complex.

1. Read `TASK.md` in this harness folder and follow it.
2. Read `../skills/reagent/SKILL.md`, then `../skills/paperclip/SKILL.md`, and Proto/Boltz skills only as needed for structure or complex checks.
3. Use the **custom immuno-risk tools** (registered by the harness extension): `list_catalytic_sites`, `structure_features`, `predict_cleavage`, `score_mhc`, `check_tolerance`, `score_immuno_risk`, or end-to-end `run_immuno_pipeline`. MHC and tolerance are stubs until real predictors/Atlas are wired — say so in outputs.
4. Work from the **repo root** (`..`): durable code in `../src/re_agent/`, artifacts in `../results/immuno_risk/`.
5. Every claim cited (Paperclip) or measured (tool JSON + versions/paths). Naturals can have holdout ground truth; de novo does not — report uncertainty.
6. Prefer MCP (`/mcp`) or CLI; do not download whole papers by hand.

Focus for this turn: ${@:-Run run_immuno_pipeline on a short natural control sequence and inspect results/immuno_risk/.}
