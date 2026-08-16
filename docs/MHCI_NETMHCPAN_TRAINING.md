# NetMHCpan MHC-I student workflow

## Scientific boundary

NetMHCpan 4.1 is used only as an offline teacher. The deployed student requires
local ESM-2 weights and its trained head, but it does not call NetMHCpan.

The full-PDA corpus build deliberately:

1. labels every 9-mer from all PDA parents before deduplication;
2. retains every unique peptide plus occurrence and parent counts;
3. connects parents sharing any exact 9-mer and keeps each connected component
   in one cross-validation fold;
4. balances five folds by unique-peptide count;
5. retains NetMHCpan EL and BA as separate channels;
6. derives strong/weak/nonbinder classes from EL percentile rank.

Teacher agreement measures distillation fidelity, not biological accuracy.

## Detached cloud corpus build

```bash
uv run --extra proto modal run --detach scripts/build_netmhcpan_corpus_modal.py \
  --run-name a0201-netmhcpan-pda-full-v1 \
  --pda-only-full \
  --parent-batch-size 20 \
  --pda-challenge-rows 0
```

This CPU job survives laptop sleep or disconnection. It writes to the persistent
Modal volume `mhci-netmhcpan-artifacts`:

```text
corpora/a0201-netmhcpan-pda-full-v1/
├── pda_training.parquet
└── manifest.json
```

Inspect progress and outputs:

```bash
uv run --extra proto modal app list
uv run --extra proto modal volume ls \
  mhci-netmhcpan-artifacts corpora/a0201-netmhcpan-pda-full-v1
```

The public IEDB API is called serially and cached in
`mhci-netmhcpan-cache`. Do not fan out public API requests from many workers.

## Detached embeddings and training

Run only after the corpus manifest and `pda_training.parquet` exist:

```bash
uv run --extra proto modal run --detach scripts/train_mhci_student_modal.py \
  --corpus-run-name a0201-netmhcpan-pda-full-v1 \
  --run-name a0201-netmhcpan-pda-cv5-v4 \
  --epochs 60 \
  --hidden-dim 512
```

The student emits three channels: EL rank propensity, BA rank propensity, and a BA
affinity score on NetMHCpan's `1 - log(IC50)/log(50000)` scale, which
`affinity_to_ic50_nm` inverts back to nanomolar. Use `--hidden-dim 512`: at the older
256 the three tasks compete and the rank channels lose about 0.014 Spearman.

Re-running training is cheap because embeddings are cached in the run directory, so
sweeping the width costs about three minutes per configuration.

The A10G job reads the cloud corpus directly, creates frozen ESM-2
`esm2_t33_650M_UR50D` embeddings, and writes:

```text
models/a0201-netmhcpan-pda-cv5-v4/
├── embeddings.float16.npy
├── embedding_manifest.json
├── rows.parquet
└── checkpoint/
    ├── deployment_manifest.json
    ├── metrics.json
    ├── fold_0/
    ├── fold_1/
    ├── fold_2/
    ├── fold_3/
    └── fold_4/
```

Each fold is held out once for testing. The next fold is used for validation and
the remaining three folds fit that head. Metrics report the mean and population
standard deviation across the five untouched test folds.

## Agent-facing processing profile

At inference, the existing checkpoint retains its independently supervised
cleavage and TAP heads. The PDA-trained student replaces only its former
MHCflurry-derived MHC lane:

- processing: N- and C-terminal cleavage probabilities;
- transport: TAP score and uncertainty;
- MHC-I: EL presentation propensity, BA binding propensity, predicted ranks,
  and strong/weak/nonbinder class;
- overall: the geometric processing score from N-cleavage, C-cleavage, and BA.

TAP and EL remain explanatory fields and are not averaged into the overall
score. Agent tools accept the legacy checkpoint through
`mhci_surrogate_checkpoint` and the CV student directory through
`mhci_netmhcpan_checkpoint`.

## Materialized occurrence-level profile

The training corpus holds one row per unique peptide, but cleavage depends on
flanking context, so the same 9-mer in two parents has two cleavage values. The
profile job re-expands the corpus to one row per occurrence:

```bash
uv run --extra proto modal run --detach scripts/build_mhci_profile_modal.py \
  --source-corpus-run a0201-netmhcpan-pda-full-v1 \
  --student-run a0201-netmhcpan-pda-cv5-v4 \
  --run-name a0201-pda-mhci-profile-v4
```

It writes `profiles/a0201-pda-mhci-profile-v4/pda_mhci_profile.parquet` plus a
manifest, 289,335 rows over 1,602 parents.

Two pooling conventions are involved and they are not interchangeable. The
student was trained on bare 9-mers and is fed the stored corpus embeddings; the
chao1 cleavage, TAP, and MHC heads were trained on flanked segments and are fed
freshly embedded `n_flank + peptide + c_flank` contexts.

Student columns are out-of-fold: each peptide is scored by the single head that
held its parent component out, so the table carries no in-sample optimism. The
deployment ensemble mean is emitted alongside it as
`student_*_propensity_ensemble` for comparison.

## Required evaluation

- Report mean and standard deviation for held-out five-fold teacher imitation.
- Compare the student, direct NetMHCpan, and the old chao1 MHC head on an
  independent experimental HLA-A\*02:01 benchmark.
- Do not call PDA an external OOD benchmark after using the full PDA corpus for
  training. An independent non-PDA cohort is required for OOD claims.
- Keep EL separate from cleavage/TAP fusion because EL already contains
  presentation information; BA is the binding-only replacement channel.

All three are now satisfied; see [`SUNDIAL_WRITEUP.md`](SUNDIAL_WRITEUP.md) for
the full results and their limits.

```bash
uv run --extra ml python scripts/compare_mhc_heads_pda.py     # old head vs new head
uv run --extra ml python scripts/benchmark_external_nyeso.py  # external cohort
```

Three findings constrain how the numbers may be used:

1. On 159,038 unique PDA 9-mers the student reaches 0.943 Spearman against the
   NetMHCpan teacher versus 0.354 for the chao1 MHCflurry head. Read the
   average-precision column, not AUROC: strong binders are 2.5% of the corpus,
   so chao1's 0.731 AUROC is only 0.060 AP against the student's 0.820.
2. The external NY-ESO cohort is pre-screened on NetMHCpan, with 93 of 95
   peptides already at EL rank <= 2.0. It supports the out-of-corpus imitation
   claim (0.909 EL Spearman and 0.898 BA Spearman). Its T-cell activation labels
   are downstream of the modeled endpoint and are not used as the primary
   benchmark.
3. BA represents peptide-HLA binding only. It does not validate cleavage,
   transport, presentation, TCR recognition, or danger. Direct BA validation
   requires measured human HLA-A*02:01 IC50 or \(K_d\); EL and the processing
   stack require independent human HLA-A*02:01 eluted-ligand or monoallelic
   immunopeptidomics evidence.

The student is systematically conservative at the strong-binder threshold, in
both the external cohort (73 calls against the teacher's 88) and the PDA profile.
On 159,038 out-of-fold PDA peptides, the conventional predicted-rank cutoff
recovers only 52.3% of 3,956 teacher strong binders. The deployed design-time
screen therefore uses a separate 95%-recall operating point rather than
redefining the standard strong-binder class.
