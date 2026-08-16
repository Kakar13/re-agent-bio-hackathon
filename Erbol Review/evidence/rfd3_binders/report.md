# RFD3 binder backbones — undruggable targets

De novo binder **backbones** generated with Proto `rfdiffusion3-design` on Modal,
using the Baker-lab / RosettaCommons foundry protein-binder settings.

## Protocol

| Setting | Value | Source |
| --- | --- | --- |
| Tool | `rfdiffusion3-design` | proto-tools |
| Device | Modal (`proto-env`, A10) | this workspace |
| `infer_ori_strategy` | `hotspots` | foundry protein_binder_design.md |
| `is_non_loopy` | `true` | same |
| `step_scale` | `3.0` | low-temperature PPI preset |
| `gamma_0` | `0.2` | low-temperature PPI preset |
| `num_timesteps` | `200` | RFD3 default |
| Binder length bins | `55-85` and `85-120` | this run |
| Contig pattern | `<len>,/0,<target_chain><lo>-<hi>` | foundry binder example |

## Targets

| Target | PDB | Hotspots | Backbones | Mean Rg (Å) | GPU time (s) |
| --- | --- | --- | ---: | ---: | ---: |
| PD-L1 (positive control) | 5O45 | A56, A115, A123 | 100 | 16.51 | 847.1 |
| KRAS4b (effector face) | 6VJJ | A37, A38, A40 | 100 | 17.72 | 1024.1 |
| MYC leucine zipper (vs MAX) | 1NKP | A927, A957, A967 | 100 | 25.44 | 738.2 |
| STAT3 SH2 dimer interface | 1BG1 | A609, A611, A613 | 100 | 18.15 | 1157.6 |
| beta-catenin ARM groove (TCF site) | 1JDH | A312, A435, A306 | 100 | 27.66 | 2380.7 |

## Per-target notes

### PD-L1 (positive control) (`pdl1`)

- Residue range: `[17, 131]`
- Contigs: `{'len_55_85': '55-85,/0,A17-131', 'len_85_120': '85-120,/0,A17-131'}`
- `select_hotspots`: `{"A56": "CG,OH", "A115": "CG,SD", "A123": "CD2,OH"}`
- Backbones: `results/rfd3_binders/pdl1/backbones/` (100 PDBs)
- Caveat: Uses foundry doc hotspots A56/A115/A123 verbatim
- Note: no partner chain — using forced/doc hotspots only
- Note: using forced/doc hotspots

### KRAS4b (effector face) (`kras`)

- Residue range: `[1, 167]`
- Contigs: `{'len_55_85': '55-85,/0,A1-167', 'len_85_120': '85-120,/0,A1-167'}`
- `select_hotspots`: `{"A37": "OE1,OE2", "A38": "OD1,OD2", "A40": "CD2,OH"}`
- Backbones: `results/rfd3_binders/kras/backbones/` (100 PDBs)
- Caveat: 6VJJ is wild-type KRAS4b–RAF1 RBD, not G12D; residue 12 is not part of the binder epitope
- Note: literature residues present at interface: [37, 38, 40]

### MYC leucine zipper (vs MAX) (`myc_max`)

- Residue range: `[897, 984]`
- Contigs: `{'len_55_85': '55-85,/0,A897-984', 'len_85_120': '85-120,/0,A897-984'}`
- `select_hotspots`: `{"A927": "OE1,NE2", "A957": "OE1,OE2", "A967": "CD1,CD2"}`
- Backbones: `results/rfd3_binders/myc_max/backbones/` (100 PDBs)
- Caveat: Cropped to MYC chain A + MAX chain B; DNA and duplicate copies removed

### STAT3 SH2 dimer interface (`stat3`)

- Residue range: `[500, 688]`
- Contigs: `{'len_55_85': '55-85,/0,A500-688', 'len_85_120': '85-120,/0,A500-688'}`
- `select_hotspots`: `{"A609": "NH1,NH2", "A611": "OG", "A613": "OG"}`
- Backbones: `results/rfd3_binders/stat3/backbones/` (100 PDBs)
- Caveat: Uses biological assembly 1 (two models merged to A/B) to recover the SH2 dimer
- Caveat: Forced literature SH2 pocket hotspots R609/S611/S613 — computed contacts favor the reciprocal pY705 face
- Note: using forced/doc hotspots
- Note: contig cropped to continuous segment A500-688 (observed 203 residues in 500-716; gaps omitted)

### beta-catenin ARM groove (TCF site) (`beta_catenin`)

- Residue range: `[151, 548]`
- Contigs: `{'len_55_85': '55-85,/0,A151-548', 'len_85_120': '85-120,/0,A151-548'}`
- `select_hotspots`: `{"A312": "NZ", "A435": "NZ", "A306": "CD2,OH"}`
- Backbones: `results/rfd3_binders/beta_catenin/backbones/` (100 PDBs)
- Caveat: Cropped ARM repeats 134-664; hotspots from beta-catenin–hTcf-4 contacts
- Note: literature residues present at interface: [312, 435]
- Note: contig cropped to continuous segment A151-548 (observed 508 residues in 135-663; gaps omitted)

## Control gate (PD-L1)

First 16 PD-L1 RFdiffusion3 complexes were checked for min heavy-atom distance
from the generated binder to the foundry hotspots Y56 / M115 / Y123
(auth numbering; output target chain is renumbered 1..N).
15 of 16 contacted all three hotspots within 4.5 Å; all 16 did so within 6 Å.
The miss is `pdl1_bb_0015.pdb` at 5.87 Å from Y56 (M115 3.69 Å, Y123 2.98 Å).
This is geometry on generated backbone complexes, not AlphaFold or ESMFold.
See `results/rfd3_binders/pdl1/control_qc.json`.

## Deploy note

Modal H100/A100-80GB required a payment method on this workspace, so the
`proto-tools-rfdiffusion3` service was deployed with `gpu=['A10:1','L4:1','T4:1']`.

## Out of scope

Sequence design (ProteinMPNN), AF2/ESMFOLD2 validation, and disordered targets
(AR-NTD, tau, α-synuclein) are deferred.

## Citation

Butcher et al. (2025) De novo Design of All-atom Biomolecular Interactions with RFdiffusion3. bioRxiv. doi:10.1101/2025.09.18.676967

Foundry binder example:
https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/examples/protein_binder_design.md

## Reproduce

```bash
uv run python scripts/prep_binder_targets.py
proto-tools deploy --apps rfdiffusion3   # once; this workspace uses A10
uv run python scripts/run_rfd3_backbones.py --target all --n 100
uv run python scripts/summarize_rfd3_backbones.py
```
