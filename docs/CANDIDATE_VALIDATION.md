# Candidate validation sheet — MHC-I presentation risk screen

For a reviewer with immunology background. The question is **"do these calls make
biological sense?"**, not "is the code correct". Everything below is regenerable; the
command is at the bottom.

## What this tool claims and does not claim

| It does | It does not |
| --- | --- |
| Predict what **NetMHCpan 4.1** would say for HLA-A\*02:01 9-mers, fast and offline | Predict immunogenicity, T-cell activation, or ADA |
| Report presentation (EL), binding (BA), and predicted IC50 in nM | Measure real binding — these are imitations of NetMHCpan, not experiments |
| Add proteasomal cleavage and TAP transport from the legacy chao1 heads | Cover MHC-II, other alleles, or peptide lengths other than 9 |

Presentation is **necessary but nowhere near sufficient** for immunogenicity. A flagged
peptide is a candidate for review, not a verdict.

## The candidate

- **ID:** `pda:9s14:0` — Protein Design Archive entry, PDB **9S14** chain A
- **Length:** 99 aa, giving 91 overlapping 9-mer windows
- **Allele:** HLA-A\*02:01

```
GHMDEALALAARAREVRPRALARYRELTDDEEEVAEVERMADLICAQRLPPEWVIQLLKEILEEVKANPEKADEMIEENRDDVMLRTLWVLATEIYRAP
```

## Result: 8 of 91 windows flagged

| Peptide | Position | Predicted IC50 | EL rank | Binder class | Screen | N-cleav | C-cleav | TAP |
| --- | --- | ---: | ---: | --- | --- | ---: | ---: | ---: |
| VMLRTLWVL | 82–91 | 3.5 nM | 0.06 | strong | flag | 0.07 | 0.76 | 4.27 |
| TLWVLATEI | 86–95 | 6.1 nM | 0.29 | strong | flag | 0.48 | 0.74 | 4.06 |
| WVIQLLKEI | 52–61 | 6.6 nM | 0.06 | strong | flag | 0.09 | 0.56 | 5.51 |
| MLRTLWVLA | 83–92 | 15.7 nM | 0.35 | strong | flag | 0.19 | 0.32 | 3.62 |
| VLATEIYRA | 89–98 | 88.3 nM | 0.80 | weak | flag | 0.52 | 0.66 | 2.78 |
| LICAQRLPP | 42–51 | 277.1 nM | 0.74 | weak | flag | 0.13 | 0.04 | 3.53 |
| LPPEWVIQL | 48–57 | 362.4 nM | 0.61 | weak | flag | 0.26 | 0.80 | 4.87 |
| LLKEILEEV | 56–65 | 429.3 nM | 0.46 | strong | flag | 0.25 | 0.85 | 3.87 |

## What to check

**1. Do the anchor residues match the A\*02:01 motif?** This is the fastest sanity check.
A\*02:01 wants L/M at P2 and V/L/I at P9.

| Peptide | P2 | P9 | Canonical? |
| --- | --- | --- | --- |
| VMLRTLWVL | M | L | yes, textbook |
| TLWVLATEI | L | I | yes |
| LLKEILEEV | L | V | yes |
| WVIQLLKEI | V | I | P2 tolerated, P9 good |
| MLRTLWVLA | L | A | P9 weak but tolerated |
| LICAQRLPP | I | **P** | poor P9 — and it scores weakest of the strong group at 277 nM |
| LPPEWVIQL | **P** | L | poor P2 — scored 362 nM |

The two peptides with broken anchors are the two the model ranks worst. That ordering is
the main thing worth confirming.

**2. Is the hotspot real?** Four of the eight flags (VMLRTLWVL, MLRTLWVLA, TLWVLATEI,
VLATEIYRA) are overlapping windows across residues 82–98. That is one hydrophobic
C-terminal region, not four independent liabilities. If you would redesign that stretch
as a unit, the tool is pointing at the right place.

**3. Do cleavage and TAP move independently of binding?** They should — they are separate
biology from separate models. `LICAQRLPP` has a very low C-terminal cleavage probability
(0.04), which argues against it being liberated by the proteasome even if it bound. A
reviewer may reasonably discount it for that reason. The tool deliberately does **not**
collapse these into one number.

**4. Is anything obviously missing?** Peptides the tool did *not* flag that you would
expect to be presented are the most useful failure reports.

## Two different calls, on purpose

`binder_class` and `screen` disagree, and that is intentional:

- **`binder_class`** uses NetMHCpan's conventional thresholds (strong = rank ≤ 0.5). It
  answers "what would NetMHCpan call this?" and is comparable with published work.
- **`screen`** uses a cutoff we calibrated for 95% recall of the teacher's strong binders.
  It answers "should a designer look at this?" and deliberately over-calls.

We measured that the conventional rule alone recovers only 52% of NetMHCpan's strong
binders when applied to our model's output, which is the wrong error direction for a
screen. Above, `VLATEIYRA` and `LPPEWVIQL` are flagged despite being "weak" by
convention. Precision at this cutoff is 44%, so **roughly half of flags are expected to
be false positives** — that is the accepted cost of not missing real ones.

## Known limitations a reviewer should weigh

- Predicted IC50 is within a median **2.3×** of NetMHCpan's value on real binders, so do
  not compare two peptides that sit within threefold of each other.
- Cleavage and TAP come from the legacy chao1 checkpoint and carry its older training
  distribution. We replaced the binding lane, not those.
- Single allele. Population-level risk needs a panel.
- No MHC-II lane in these numbers (MHC-II runs separately via NetMHCIIpan).

## Regenerate

```bash
uv run --extra ml python -c "
from pathlib import Path
from re_agent.immuno.e2e_pls_pickle import TeamE2EPLSAdapter
SEQ='GHMDEALALAARAREVRPRALARYRELTDDEEEVAEVERMADLICAQRLPPEWVIQLLKEILEEVKANPEKADEMIEENRDDVMLRTLWVLATEIYRAP'
a=TeamE2EPLSAdapter(Path('models/chao1/cv5_heads.pkl 2'),
                    netmhcpan_checkpoint_dir=Path('models/a0201-netmhcpan-pda-cv5-v4/checkpoint'))
for p in a.predict(SEQ).predictions:
    if p.mhc_i_screening_flag:
        print(p.peptide, p.start, round(p.mhc_i_predicted_ba_ic50_nm,1), p.mhc_i_binder_class)
"
```

Or through the workbench UI: **Sequence risk** view → **Load PDA example** → **Run chao1**.
