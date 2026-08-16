# Binder target inputs for RFdiffusion3

Prepared by `scripts/prep_binder_targets.py`. Large PDBs are gitignored; regenerate locally.

| Target | PDB | Epitope | Hotspots |
| --- | --- | --- | --- |
| `pdl1` | 5O45 | IgV domain A17–131 | A56, A115, A123 (RFD3 doc) |
| `kras` | 6VJJ | Switch I/II vs RAF RBD | A37, A38, A40 |
| `myc_max` | 1NKP | MYC LZ vs MAX | computed zipper contacts |
| `stat3` | 1BG1 asm1 | SH2 pocket (dimer) | A609, A611, A613 |
| `beta_catenin` | 1JDH | ARM groove vs TCF | A312, A435, A306 |

Registry: [`configs/rfd3_targets.json`](../configs/rfd3_targets.json).
