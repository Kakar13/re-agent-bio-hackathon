# Late-stage immuno-risk (Track C)

Read [`../TASK.md`](../TASK.md) and [`../../docs/IMMUNO_RISK_DESIGN.md`](../../docs/IMMUNO_RISK_DESIGN.md).

## Mission

Screen de novo / natural candidates for **intracellular / plasmid delivery**:

1. **MHC-I full** — MHCflurry presentation + calibrated IEDB risk head (± optional NetMHCpan-4.2e)
2. **MHC-II thin** — presentation only (HLAIIPred / heuristic; ± NetMHCIIpan-4.3k)
3. **Atlas tolerance evidence** — separate from presentation
4. **Aggregation report** — separate score; do not infer from protease accessibility
5. **Artifacts + LangSmith** — judges must re-run and inspect
6. **Benchling** — pull candidates / publish run summaries when credentials exist

## Working rules

- Prefer `run_immuno_pipeline` for end-to-end; cite predictor versions and caveats.
- Never claim clinical immunogenicity probability or ADA from MHC-II ranks.
- Cleavage tools are diagnostic; MHC-I peptide ranking comes from MHCflurry.
- Benchling publish is an explicit external action (`dry_run` first).
- Fail closed / document gaps; do not invent API keys or silently reintroduce stubs.
- Ask before git commits.
