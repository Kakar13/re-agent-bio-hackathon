"""Attention-pooling immunogenicity head over frozen ESM-2 embeddings.

Small enough to retrain in seconds, and the pooling attention doubles as the
per-residue explanation: the weight a residue receives is literally how much it
contributed to the pooled vector the classifier saw.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from re_agent.immuno.config import ESM_DIM, TrainConfig


class ImmunoHead(nn.Module):
    """Project -> self-attention -> attention-pool -> MLP -> logit."""

    def __init__(self, cfg: TrainConfig, in_dim: int = ESM_DIM):
        super().__init__()
        h = cfg.hidden
        self.input_norm = nn.LayerNorm(in_dim)
        self.project = nn.Linear(in_dim, h)
        # Attention-internal dropout is left at zero: MPS cannot run dropout inside
        # scaled_dot_product_attention under no_grad, which is exactly the path the
        # EMA teacher and MC-dropout take. The explicit Dropout layers below supply
        # the stochasticity those two need.
        self.attn = nn.MultiheadAttention(h, cfg.heads, dropout=0.0, batch_first=True)
        self.attn_norm = nn.LayerNorm(h)
        self.ffn = nn.Sequential(
            nn.Linear(h, h * 2), nn.GELU(), nn.Dropout(cfg.dropout), nn.Linear(h * 2, h)
        )
        self.ffn_norm = nn.LayerNorm(h)
        self.pool_score = nn.Linear(h, 1)
        self.dropout = nn.Dropout(cfg.dropout)
        self.classifier = nn.Sequential(
            nn.Linear(h, h // 2), nn.GELU(), nn.Dropout(cfg.dropout), nn.Linear(h // 2, 1)
        )

    def forward(
        self, x: torch.Tensor, mask: torch.Tensor, return_attention: bool = False
    ) -> torch.Tensor | tuple[torch.Tensor, torch.Tensor]:
        """`x`: (B, L, in_dim) embeddings. `mask`: (B, L) with 1 for real residues."""
        pad = mask < 0.5
        h = self.project(self.input_norm(x))
        attended, _ = self.attn(h, h, h, key_padding_mask=pad, need_weights=False)
        h = self.attn_norm(h + self.dropout(attended))
        h = self.ffn_norm(h + self.dropout(self.ffn(h)))

        scores = self.pool_score(h).squeeze(-1).masked_fill(pad, float("-inf"))
        weights = torch.softmax(scores, dim=-1)
        pooled = torch.einsum("bl,blh->bh", weights, h)
        logit = self.classifier(self.dropout(pooled)).squeeze(-1)
        if return_attention:
            return logit, weights
        return logit


@torch.no_grad()
def predict(
    model: nn.Module, x: torch.Tensor, mask: torch.Tensor, return_attention: bool = False
) -> tuple[torch.Tensor, torch.Tensor | None]:
    model.eval()
    if return_attention:
        logit, weights = model(x, mask, return_attention=True)
        return torch.sigmoid(logit), weights
    return torch.sigmoid(model(x, mask)), None


@torch.no_grad()
def mc_dropout_predict(
    model: nn.Module, x: torch.Tensor, mask: torch.Tensor, passes: int = 20
) -> tuple[torch.Tensor, torch.Tensor]:
    """Sample the predictive distribution by keeping dropout active at inference.

    The spread across passes is the model's own uncertainty about this window.
    """
    was_training = model.training
    model.train()  # re-enables dropout; there are no batchnorm layers to corrupt
    probs = torch.stack([torch.sigmoid(model(x, mask)) for _ in range(passes)])
    model.train(was_training)
    return probs.mean(0), probs.std(0)


def saliency(
    model: nn.Module, x: torch.Tensor, mask: torch.Tensor, steps: int = 32
) -> torch.Tensor:
    """Integrated gradients over the embedding, reduced to one value per residue.

    Cross-checks the pooling attention: attention says where the model looked,
    integrated gradients say what actually moved the logit.
    """
    model.eval()
    baseline = torch.zeros_like(x)
    total = torch.zeros_like(x)
    for step in range(1, steps + 1):
        point = (baseline + (x - baseline) * (step / steps)).detach().requires_grad_(True)
        logit = model(point, mask)
        grad = torch.autograd.grad(logit.sum(), point)[0]
        total += grad
    attribution = ((x - baseline) * total / steps).sum(-1)
    return attribution * mask


def count_parameters(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def bce_loss(logit: torch.Tensor, target: torch.Tensor, pos_weight: float | None = None):
    weight = torch.tensor(pos_weight, device=logit.device) if pos_weight else None
    return F.binary_cross_entropy_with_logits(logit, target, pos_weight=weight)
