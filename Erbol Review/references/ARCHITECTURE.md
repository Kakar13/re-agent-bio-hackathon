# Architecture: what is trained, what is frozen, and what changed

Reference for the CHAO / re:AGENT MHC-I lane. The short version: one large frozen
encoder is shared by every head in the system, and the only thing we trained this
weekend is a 660K-parameter readout that replaces one of four heads.

## 1. The shared frozen encoder

Every MHC-I score in this system starts from the same place:

- **Model:** `esm2_t33_650M_UR50D` (Meta, ESM-2), 650 million parameters
- **Status:** frozen. It never receives a gradient, from us or from chao1's authors.
- **Output:** final-layer (33) per-residue embeddings, 1,280 dimensions per residue

Fine-tuning it was considered and explicitly excluded in `docs/DATA_PLAN.md`: 160,012
unique peptides cannot support updating 650M parameters without severe overfitting, and
the cost is not justified when the representation already encodes the hydrophobicity and
positional chemistry that A\*02:01 anchor preference depends on.

**Pooling differs by head, and the difference is load-bearing.** The same encoder output
is pooled three ways:

| Consumer | Pooling window | Why |
| --- | --- | --- |
| Cleavage head | ±3 residues around the N-terminal site, and ±3 around the C-terminal site | Cleavage is a property of the junction, not the peptide |
| TAP head | mean over the 9 peptide residues, embedded **with 4-residue flanks** | Transport depends on the peptide, but the model was trained in flanked context |
| NetMHCpan student | mean over the 9 peptide residues of a **bare 9-mer** | The teacher labels bare 9-mers, so the student must see what the teacher saw |

The student and chao1 therefore require two separate encoder passes over different
inputs. They are not interchangeable embeddings, and the profile builder runs both lanes
deliberately rather than reusing one.

## 2. What chao1 was predicting

`chao1` (`models/chao1/cv5_heads.pkl 2`, bundle format `e2e_pls-heads-pickle-v1`) is a
**four-output antigen-processing surrogate** for HLA-A\*02:01 9-mers. It models the
MHC-I presentation pathway as a sequence of biological steps:

```
protein → [proteasome cleaves it] → [TAP transports the fragment] → [MHC-I binds it] → surface
              chao1 cleavage head        chao1 TAP head              chao1 MHC head
```

| Head | Architecture | Output | Kept? |
| --- | --- | --- | --- |
| **Cleavage** | MLP: `Linear(2560→128) → ReLU → Linear(128→64) → ReLU → Linear(64→2)`, sigmoid, then isotonic calibration | N- and C-terminal cleavage probability | **Kept** |
| **TAP** | Linear regression on the pooled 9-mer, plus a bootstrap ensemble of coefficients | Relative TAP log-IC50 and its uncertainty | **Kept** |
| **MHC presentation** | Learned linear projection → L2 normalize → cosine similarity against a per-allele centroid → isotonic calibration | Presentation propensity | **Replaced** |

Note the third row. The chao1 MHC head is not a regression head — it is a
**metric-learning / centroid-cosine** design. It projects the peptide embedding into a
space, then asks "how close is this to the learned A\*02:01 binder centroid?" A single
centroid and a monotonic calibrator is very low capacity for a binding-specificity
problem, which is a structural reason, independent of its training data, that it
struggles to resolve the rare strong binders.

Its labels also came from **MHCflurry**, and its training corpus had been sampled by
MHCflurry quantile — so the peptide distribution was shaped by the very model being
imitated.

## 3. What we replaced it with

The NetMHCpan student is a plain regression head, trained by knowledge distillation.
NetMHCpan 4.1 is the **teacher**: we called it on every 9-mer in the PDA corpus to
generate labels, then trained this head to reproduce them from the frozen embedding.

```
LayerNorm(1280) → Linear(1280→512) → GELU → Dropout(0.1) → Linear(512→3) → sigmoid
```

- **659,971 trainable parameters** — roughly 0.1% the size of the encoder in front of it
- 2.6 MB per fold on disk
- Deployed as **five heads**, one per cross-validation fold, averaged into an ensemble

Three output channels:

| Channel | Target | Interpretation |
| --- | --- | --- |
| `netmhcpan_el_rank_propensity` | monotonic transform of EL percentile rank | presentation |
| `netmhcpan_ba_rank_propensity` | monotonic transform of BA percentile rank | binding |
| `netmhcpan_ba_affinity_score` | NetMHCpan's `1 - log(IC50)/log(50000)` | inverts to nM: `IC50 = 50000^(1-score)` |

Ranks are used for the first two because raw EL scores pile up at zero and give the
optimizer almost no gradient where the ordering matters. The affinity channel uses the
raw bounded scale because it is well conditioned — only 0.02% of the corpus is censored
at the 50,000 nM ceiling.

## 4. Head-to-head

Both scored on the same 159,038 unique PDA 9-mers, against the same NetMHCpan EL teacher,
neither in sample:

| | chao1 MHC head | NetMHCpan student |
| --- | --- | --- |
| Trained to imitate | MHCflurry | NetMHCpan 4.1 |
| Head type | centroid cosine + calibrator | 2-layer MLP regression |
| Outputs | 1 (presentation propensity) | 3 (EL, BA, affinity in nM) |
| Spearman vs teacher | 0.354 | **0.943** |
| Strong-binder AUROC | 0.731 | **0.993** |
| Strong-binder average precision | 0.060 | **0.820** |

Average precision is the honest column. Strong binders are 2.5% of the corpus, so
chao1's 0.731 AUROC corresponds to an average precision of 0.060 — near-unusable for the
actual task of surfacing the rare risky peptide. The student is roughly 14× better on the
metric that matters under this class imbalance.

**This is a cross-teacher disagreement, not an error rate against experimental data.**
chao1 is being measured against a teacher it was never trained on. The claim is that
NetMHCpan is the better teacher and that we now imitate it faithfully — not that chao1's
authors built a broken model.

## 5. The deployed composite

What actually ships is a hybrid: chao1 for the two upstream biological steps, the student
for the binding step.

```
peptide window (9-mer, with 4-residue flanks for context)
  │
  ├─ frozen ESM-2 650M ── flanked pooling ──┬─ chao1 cleavage head → N/C probabilities
  │                                          └─ chao1 TAP head      → transport + uncertainty
  │
  └─ frozen ESM-2 650M ── bare 9-mer pooling ─ NetMHCpan student ×5 ┬─ EL propensity
                                                                     ├─ BA propensity
                                                                     └─ affinity → IC50 nM
                                                                          │
                        derived: binder_class (NetMHCpan convention), screening_flag (95% recall)
```

Runtime adapter ID: `team-e2e-pls-chao1-netmhcpan`. When the student checkpoint is absent
the adapter degrades to `team-e2e-pls-chao1` and the chao1 MHC head serves the binding
lane, so a fresh clone still runs.

`overall_mhci_risk` is the EL propensity alone. The composite processing score combines
N-cleavage, C-cleavage, and BA only — EL is deliberately excluded, because it already
encodes cleavage and transport information and averaging it in would count the same
evidence twice.

## 6. Where the runtime cost sits

In the 0.76 s measurement for scanning a 156-residue design (147 windows), essentially
all of it is the frozen ESM-2 forward pass. The 660K-parameter head is negligible. This
is why cached embeddings make rescoring nearly free, and why the deployed model could be
retrained and swapped in about three minutes without re-embedding the corpus.

## 7. Boundaries

- The student imitates NetMHCpan and **cannot exceed it**. 0.943 Spearman is a ceiling
  being approached, not a bar being cleared.
- Predicted IC50 is within a median 2.3× of the teacher on real binders — do not compare
  two peptides within threefold of each other.
- Cleavage and TAP are inherited from chao1 and carry its older training distribution.
  We replaced the binding lane, not the whole model.
- Single allele (HLA-A\*02:01), 9-mers only.
- Nothing here is trained on TCR activation or ADA labels, so no immunogenicity axis is
  reported.

## Source files

| Component | Path |
| --- | --- |
| Student model, training, calibration | `src/re_agent/e2e_pls/netmhcpan_student.py` |
| chao1 heads, encoder, deployed adapter | `src/re_agent/immuno/e2e_pls_pickle.py` |
| Corpus construction and fold assignment | `src/re_agent/e2e_pls/netmhcpan_corpus.py` |
| Teacher client (IEDB NetMHCpan 4.1) | `src/re_agent/immuno/netmhcpan.py` |
| Deployed manifest and CV metrics | `models/a0201-netmhcpan-pda-cv5-v4/checkpoint/` |
