"""Two biochemistry sanity checks on the trained head.

Nothing in training mentions binding-groove chemistry: the model sees only ESM-2
embeddings and a binary assay outcome.

1. Does pooling *attention* concentrate on large hydrophobics (class II P1 anchors)?
2. Does the predicted *risk* track hydrophobic content the way the true labels do?

These ask different things. Attention is where the head aggregates evidence, which need
not be where the peptide contacts the groove — ESM-2 embeddings are contextual, so a
single position already carries its neighbours. The risk-vs-hydrophobicity check is the
one that tests whether the model reproduces the real biological trend.
"""

from __future__ import annotations

import json

import numpy as np
import torch

from re_agent.immuno.config import PATHS
from re_agent.immuno.data import build_labeled
from re_agent.immuno.embed import load_embeddings
from re_agent.immuno.explain import HYDROPHOBIC, anchor_enrichment, window_attention
from re_agent.immuno.model import predict
from re_agent.immuno.score import load_bundle


def main(n_windows: int = 20_000, seed: int = 0) -> None:
    bundle = load_bundle("mean_teacher")
    labeled = build_labeled().reset_index(drop=True)
    emb, mask = load_embeddings("labeled")

    rng = np.random.default_rng(seed)
    rows = rng.choice(len(labeled), size=min(n_windows, len(labeled)), replace=False)
    rows.sort()
    windows = labeled["seq"].to_numpy()[rows].tolist()

    x = torch.from_numpy(np.array(emb[rows], dtype=np.float32)).to(bundle.device)
    m = torch.from_numpy(np.array(mask[rows], dtype=np.float32)).to(bundle.device)
    attention = window_attention(bundle.teacher, x, m)

    report = {"attention": anchor_enrichment(attention, windows)}
    labels = labeled["label"].to_numpy()[rows].astype(float)

    probs, _ = predict(bundle.teacher, x, m)
    risks = np.asarray(probs.cpu(), dtype=np.float64)
    hydro = np.array([sum(c in HYDROPHOBIC for c in w) / len(w) for w in windows])
    aromatic = np.array([sum(c in "FWY" for c in w) / len(w) for w in windows])

    report["biochemistry"] = {
        "n": int(len(rows)),
        "risk_vs_hydrophobic_fraction_r": round(float(np.corrcoef(risks, hydro)[0, 1]), 4),
        "label_vs_hydrophobic_fraction_r": round(float(np.corrcoef(labels, hydro)[0, 1]), 4),
        "risk_vs_aromatic_fraction_r": round(float(np.corrcoef(risks, aromatic)[0, 1]), 4),
        "label_vs_aromatic_fraction_r": round(float(np.corrcoef(labels, aromatic)[0, 1]), 4),
        "mean_risk_hydrophobic_top_quartile": round(
            float(risks[hydro >= np.quantile(hydro, 0.75)].mean()), 4
        ),
        "mean_risk_hydrophobic_bottom_quartile": round(
            float(risks[hydro <= np.quantile(hydro, 0.25)].mean()), 4
        ),
    }

    path = PATHS.reports / "anchor_enrichment.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2))

    attn = report["attention"]
    bio = report["biochemistry"]
    print(f"windows analysed: {len(rows)}\n")
    print("[1] attention vs class II anchor chemistry")
    print(f"  top-5 attended residues: {attn['top5']}")
    print(f"  hydrophobic (FWYLIVM) mean attention ratio: {attn['hydrophobic_mean_ratio']}")
    print(f"  all other residues    mean attention ratio: {attn['other_mean_ratio']}")
    for aa, ratio in attn["per_residue_attention"].items():
        flag = " <- hydrophobic" if aa in HYDROPHOBIC else ""
        print(f"    {aa}  {ratio:.3f}{flag}")

    print("\n[2] predicted risk vs hydrophobic content (does the model track the real trend?)")
    print(
        f"  hydrophobic fraction:  risk r={bio['risk_vs_hydrophobic_fraction_r']:+.3f}   "
        f"true label r={bio['label_vs_hydrophobic_fraction_r']:+.3f}"
    )
    print(
        f"  aromatic fraction:     risk r={bio['risk_vs_aromatic_fraction_r']:+.3f}   "
        f"true label r={bio['label_vs_aromatic_fraction_r']:+.3f}"
    )
    print(
        f"  mean risk, most hydrophobic quartile {bio['mean_risk_hydrophobic_top_quartile']:.3f} "
        f"vs least {bio['mean_risk_hydrophobic_bottom_quartile']:.3f}"
    )
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
