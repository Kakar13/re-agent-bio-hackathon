# Design-to-screen agent orchestration

## Goal

The workbench accepts a natural-language protein-design objective and preserves an
inspectable artifact chain from literature evidence through candidate screening.
Historical IL-7Ralpha sequences are used only for a no-GPU preflight. They are not
generation inputs.

## Runtime flow

1. `research_design_objective`
   - Loads Paperclip's current skill instructions.
   - Routes the design intent when a routine is available.
   - Searches literature and protein records.
   - Writes a cited `target_research` artifact.
2. `inspect_proto_design_tools`
   - Uses a LangGraph-owned MCP client.
   - Checks the Proto workspace.
   - Discovers RFdiffusion3, ProteinMPNN, AlphaFold2, and PyMOL RMSD alignment.
   - Records the exact input, config, and output schemas.
3. `draft_cited_design_spec`
   - Accepts evidence-backed target metadata and source URLs.
   - Pins a local PDB by SHA-256.
   - Verifies every hotspot against the requested PDB chain.
   - Bounds binder length and generation budget.
   - Writes a versioned design specification under `results/design_specs/`.
4. `plan_design_campaign`
   - Converts the approved specification into exact Proto payloads.
   - Writes a dry-run manifest without starting compute.
5. `execute_design_campaign`
   - Pauses at a LangGraph `interrupt` before any GPU side effect.
   - Shows backbone, sequence, and AlphaFold2 call counts in the approval payload.
   - Resumes only with an explicit approval object.
   - Dispatches typed RFdiffusion3, ProteinMPNN, AlphaFold2, and RMSD calls to
     the authenticated Modal backend through runtime MCP.
   - Computes binder-chain RMSD with Proto and an explicit inter-chain clash count locally.
   - Applies structural validation before immunogenicity screening.
6. Immunogenicity handoff
   - Screens only candidates whose validation status is `pass`.
   - Preserves response, NetMHCIIpan EL, NetMHCIIpan BA, processing, tolerance,
     and the optional team MHC-I processing surrogate as separate channels.
   - Never fuses the HLA-A*02:01 9-mer MHC-I surrogate into the MHC-II score.
   - Withholds the combined rank when the calibrated response-model artifact is
     unavailable.
   - Maps residue tracks to Mol* only after exact PDB-chain sequence validation.
7. Deterministic review
   - Checks citations and claim boundaries.
   - Rejects any screened candidate that did not pass structural validation.
   - Reconciles candidate and assessment counts.
   - Enforces null combined scores when required evidence is unavailable.

## Structural validation policy

Every generated candidate records six checks:

- unique sequence
- ProteinMPNN perplexity present
- monomer pLDDT present and at least 80
- designed-versus-refolded binder RMSD present and at most 1 angstrom
- interface PAE present and at most 10
- clash check present and passed

Missing metrics fail closed. A campaign with generated candidates but no structural
passes is reported as `validation_blocked`, not as a successful screening run.

These thresholds are computational triage priors. They are not proof of folding,
binding, or immune response.

## Reference preflight

`run_reference_campaign_preflight` reads a bounded subset of:

- `data/raw/IL7Ra_binders_sequences.fasta`
- `data/raw/il7ra/IL7Ra_binders_screening_results.csv`
- `data/raw/il7ra/bli_binding_data.csv`

It joins historical monomer, binder, and BLI metadata to the screening artifact, then
runs the same NetMHCIIpan EL/BA and response-adapter path used for generated candidates.
No RFdiffusion3, ProteinMPNN, or AlphaFold2 call occurs in this lane.

The primary `models/chao1/cv5_heads.pkl 2` checkpoint can also be attached to this preflight or a
single-candidate screen. Its restricted adapter runs the matching ESM-2 650M encoder,
preserves cleavage, TAP, and MHC-I outputs, and adds two residue tracks without changing
the MHC-II fusion result.

## MCP and credential boundary

Cursor, Pi, and the LangGraph API are separate processes with separate MCP sessions.
The graph therefore owns its MCP connections in
`src/re_agent/agent/mcp_runtime.py`.

- Proto is available through the hosted bearer-token endpoint when
  `PROTO_API_KEY` is set, or through the local `proto-tools-mcp` stdio server.
- Paperclip uses its hosted endpoint when `PAPERCLIP_MCP_BEARER_TOKEN` is set.
- Without that token, the graph can use an installed, authenticated `paperclip` CLI.
- Cursor's Paperclip OAuth session is not copied or scraped. If neither graph
  authentication path exists, research fails with an explicit configuration error.

No credential is written to an artifact, trace, or UI payload.

## Current validation

- Python orchestration, reviewer, MCP status, and structural-gate tests: passing.
- Next.js lint and production build: passing.
- Playwright desktop and mobile system suite: passing.
- LangGraph-owned Proto MCP schema discovery: passing for all three campaign tools.
- Modal profile `mihirjoe` is authenticated, `proto-env` exists, and RFdiffusion3,
  ProteinMPNN, and AlphaFold2 are deployed and available.
- One-sequence reference preflight: NetMHCIIpan status `ok`; combined rank correctly
  withheld without a response artifact.
- PDA chao1 smoke: chao1 scored 91 windows for `pda:9s14:0`, with top-five mean
  processing composite 0.6975 and maximum 0.7053 at `WVLATEIYR`.
  NetMHCIIpan EL/BA also returned `ok`, while the combined rank remained null because
  the CD4 response adapter is still unavailable. The parent and all 91 9-mers were
  absent from the supplied teammate training parquet. This is an inference smoke test,
  not a reproduced held-out performance estimate.
- Deterministic review of that PDA artifact: all seven scientific gates passed,
  including MHC-I/MHC-II separation and residue-track alignment.
- Paperclip CLI fallback: installed and authenticated from the merged workbench
  environment; target research artifacts retain the returned search sets and citations.

## Remaining live-run gates

1. Run one low-budget RFdiffusion3 backbone after explicit UI approval.
2. Confirm the returned AlphaFold2 metric keys match the supported validation aliases.
3. Supply a frozen calibrated response-model artifact, or keep the combined rank null.
4. Add MixMHC2pred only as a separate second-opinion provider. It must not be merged
   into the NetMHCIIpan EL/BA record.
5. Run the bounded campaign, inspect the Mol* mapping, and retain the LangSmith trace
   plus local hashed manifests for the demo.
