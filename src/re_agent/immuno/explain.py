"""Per-residue heatmaps: pooling attention and an integrated-gradients cross-check.

Attention says where the model looked inside a window; integrated gradients say
what actually moved the logit. Reporting both means a hotspot that only one method
supports is visibly weaker evidence than one both agree on.
"""

from __future__ import annotations

import numpy as np
import torch

from re_agent.immuno.config import WINDOW
from re_agent.immuno.model import saliency


@torch.no_grad()
def window_attention(model, x: torch.Tensor, mask: torch.Tensor, batch: int = 2048) -> np.ndarray:
    """(n_windows, WINDOW) pooling-attention weights."""
    model.eval()
    out = []
    for i in range(0, len(x), batch):
        _, weights = model(x[i : i + batch], mask[i : i + batch], return_attention=True)
        out.append(weights.cpu().numpy())
    return np.concatenate(out) if out else np.zeros((0, WINDOW))


def _scatter(values: np.ndarray, starts: np.ndarray, length: int) -> np.ndarray:
    """Average overlapping window contributions back onto the protein."""
    total = np.zeros(length, dtype=np.float64)
    counts = np.zeros(length, dtype=np.float64)
    for row, start in zip(values, starts):
        span = min(len(row), length - int(start))
        if span <= 0:
            continue
        total[int(start) : int(start) + span] += row[:span]
        counts[int(start) : int(start) + span] += 1.0
    return total / np.clip(counts, 1.0, None)


def per_residue_heatmap(
    model,
    x: torch.Tensor,
    mask: torch.Tensor,
    risks: np.ndarray,
    starts: np.ndarray,
    length: int,
) -> np.ndarray:
    """Risk-weighted attention: a residue is hot if risky windows leaned on it.

    Attention alone would highlight the focal residue of a benign window just as
    strongly as that of a dangerous one, so each window's attention profile is
    scaled by that window's risk before being scattered back.
    """
    attn = window_attention(model, x, mask)
    weighted = attn * risks[:, None]
    heat = _scatter(weighted, starts, length)
    peak = heat.max()
    return heat / peak if peak > 0 else heat


def per_residue_saliency(
    model, x: torch.Tensor, mask: torch.Tensor, starts: np.ndarray, length: int, batch: int = 512
) -> np.ndarray:
    """Integrated gradients, aggregated per residue and normalized to [0, 1]."""
    chunks = []
    for i in range(0, len(x), batch):
        attribution = saliency(model, x[i : i + batch], mask[i : i + batch])
        chunks.append(attribution.detach().cpu().numpy())
    values = np.concatenate(chunks) if chunks else np.zeros((0, WINDOW))
    heat = _scatter(values, starts, length)
    peak = np.abs(heat).max()
    return heat / peak if peak > 0 else heat


def top_hotspots(heat: np.ndarray, sequence: str, k: int = 5, min_gap: int = 8) -> list[dict]:
    """Local maxima of the heatmap, separated so one broad peak is not counted twice."""
    order = np.argsort(-heat)
    picked: list[int] = []
    for i in order:
        if all(abs(int(i) - p) >= min_gap for p in picked):
            picked.append(int(i))
        if len(picked) >= k:
            break
    return [
        {
            "position": p + 1,
            "residue": sequence[p],
            "score": round(float(heat[p]), 4),
            "context": sequence[max(0, p - 7) : p + 8],
        }
        for p in sorted(picked)
    ]


def plot_heatmap(score, path, width: float = 14.0):
    """Sequence-track figure: per-residue risk, saliency, and the called regions."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    seq = score.sequence
    n = len(seq)
    fig, axes = plt.subplots(
        3, 1, figsize=(width, 6.0), sharex=True, height_ratios=[1.0, 1.0, 0.7]
    )

    axes[0].imshow(
        score.per_residue[None, :], aspect="auto", cmap="inferno", vmin=0, vmax=1,
        extent=(0.5, n + 0.5, 0, 1),
    )
    axes[0].set_yticks([])
    axes[0].set_ylabel("attention\nrisk", rotation=0, ha="right", va="center")
    axes[0].set_title(
        f"{score.name} — immunogenicity risk {score.risk:.3f} "
        f"({score.risk_percentile:.0f}th percentile vs natural), "
        f"confidence {score.confidence:.2f}"
    )

    sal = score.per_residue_saliency
    limit = max(float(np.abs(sal).max()), 1e-9)
    axes[1].imshow(
        sal[None, :], aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit,
        extent=(0.5, n + 0.5, 0, 1),
    )
    axes[1].set_yticks([])
    axes[1].set_ylabel("integrated\ngradients", rotation=0, ha="right", va="center")

    axes[2].plot(np.arange(1, n + 1), score.per_residue, color="black", lw=1.0)
    for region in score.regions[:8]:
        axes[2].axvspan(region.start + 1, region.end, color="crimson", alpha=0.18)
        axes[2].text(
            (region.start + region.end) / 2, 1.02, f"{region.risk:.2f}",
            ha="center", va="bottom", fontsize=7, color="crimson",
        )
    axes[2].set_ylim(0, 1.15)
    axes[2].set_ylabel("risk", rotation=0, ha="right", va="center")
    axes[2].set_xlabel("residue position")
    axes[2].set_xlim(0.5, n + 0.5)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
