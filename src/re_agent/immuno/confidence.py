"""Decomposed confidence: stability, student/teacher agreement, and familiarity.

A single opaque number is not useful for a de novo binder, because the honest
answer is often "this sequence is unlike anything the model was trained on".
Each component is reported separately so that case stays visible.
"""

from __future__ import annotations

import numpy as np
import torch
from sklearn.neighbors import NearestNeighbors

from re_agent.immuno.model import mc_dropout_predict

# A Bernoulli's standard deviation maxes out at 0.5, which is the worst possible
# MC-dropout spread; normalizing by it puts stability on a 0-1 scale.
MAX_BERNOULLI_STD = 0.5


def temperature_scale(logits: np.ndarray, labels: np.ndarray, max_iter: int = 200) -> float:
    """Fit a single scalar temperature by minimizing validation NLL."""
    log_t = torch.zeros(1, requires_grad=True)
    z = torch.from_numpy(logits.astype(np.float32))
    y = torch.from_numpy(labels.astype(np.float32))
    opt = torch.optim.LBFGS([log_t], lr=0.1, max_iter=max_iter)

    def closure():
        opt.zero_grad()
        loss = torch.nn.functional.binary_cross_entropy_with_logits(z / log_t.exp(), y)
        loss.backward()
        return loss

    opt.step(closure)
    return float(log_t.exp().item())


def apply_temperature(logits: np.ndarray, temperature: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-logits / max(temperature, 1e-6)))


class OODIndex:
    """kNN distance to labeled training windows in ESM-2 space.

    De novo binders are systematically far from natural training peptides; this
    turns that distance into an explicit familiarity score.
    """

    def __init__(
        self,
        train_emb: np.ndarray,
        mask: np.ndarray,
        n_ref: int = 20_000,
        k: int = 10,
        seed: int = 0,
    ):
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(train_emb), size=min(n_ref, len(train_emb)), replace=False)
        idx.sort()
        self.k = k
        self.reference = self._pool(np.asarray(train_emb[idx]), np.asarray(mask[idx]))
        self.nn = NearestNeighbors(n_neighbors=k, metric="cosine").fit(self.reference)
        # Self-distances of the reference set define what "familiar" looks like.
        self_dist, _ = self.nn.kneighbors(self.reference, n_neighbors=k + 1)
        self.baseline = np.sort(self_dist[:, 1:].mean(axis=1))

    @staticmethod
    def _pool(emb: np.ndarray, mask: np.ndarray) -> np.ndarray:
        weights = mask[..., None].astype(np.float32)
        total = np.clip(weights.sum(axis=1), 1e-6, None)
        return (emb.astype(np.float32) * weights).sum(axis=1) / total

    def distance(self, emb: np.ndarray, mask: np.ndarray) -> np.ndarray:
        dist, _ = self.nn.kneighbors(self._pool(emb, mask), n_neighbors=self.k)
        return dist.mean(axis=1)

    def familiarity(self, emb: np.ndarray, mask: np.ndarray) -> np.ndarray:
        """Fraction of training windows at least as isolated as this one."""
        d = self.distance(emb, mask)
        return 1.0 - np.searchsorted(self.baseline, d) / len(self.baseline)


def window_confidence(
    student,
    teacher,
    x: torch.Tensor,
    mask: torch.Tensor,
    ood: OODIndex | None = None,
    emb_np: np.ndarray | None = None,
    mask_np: np.ndarray | None = None,
    passes: int = 20,
) -> dict[str, np.ndarray]:
    """Per-window confidence components, each in [0, 1] where 1 means confident."""
    _, std_prob = mc_dropout_predict(teacher, x, mask, passes=passes)
    with torch.no_grad():
        teacher.eval()
        student.eval()
        p_teacher = torch.sigmoid(teacher(x, mask))
        p_student = torch.sigmoid(student(x, mask))

    stability = 1.0 - (std_prob / MAX_BERNOULLI_STD).clamp(0, 1)
    agreement = 1.0 - (p_student - p_teacher).abs().clamp(0, 1)

    out = {
        "risk": p_teacher.cpu().numpy(),
        "stability": stability.cpu().numpy(),
        "agreement": agreement.cpu().numpy(),
    }
    if ood is not None and emb_np is not None:
        out["familiarity"] = ood.familiarity(emb_np, mask_np).astype(np.float32)
        out["ood_distance"] = ood.distance(emb_np, mask_np).astype(np.float32)
    else:
        out["familiarity"] = np.ones_like(out["risk"])
        out["ood_distance"] = np.zeros_like(out["risk"])

    out["confidence"] = (out["stability"] + out["agreement"] + out["familiarity"]) / 3.0
    return out


def aggregate_confidence(
    components: dict[str, np.ndarray], weights: np.ndarray | None = None
) -> dict:
    """Collapse per-window confidence to protein level.

    Windows are weighted by risk: confidence in a benign stretch is not what
    matters when the headline number is driven by the risky ones.
    """
    if weights is None:
        weights = np.ones_like(components["risk"])
    total = float(weights.sum())
    w = weights / total if total > 1e-9 else np.full_like(weights, 1.0 / len(weights))

    parts = ("stability", "agreement", "familiarity")
    summary = {key: float((components[key] * w).sum()) for key in parts}
    summary["confidence"] = float(np.mean(list(summary.values())))
    summary["mean_ood_distance"] = float(components["ood_distance"].mean())
    return summary
