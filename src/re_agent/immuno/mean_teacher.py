"""Mean Teacher: EMA teacher, embedding perturbations, and the consistency term.

The teacher is an exponential moving average of the student's own weights, so no
label is ever invented and no external predictor is imitated. The unlabeled de novo
windows only ever contribute the requirement that the model answer them the same
way twice under different perturbations.
"""

from __future__ import annotations

import copy
import math

import torch
import torch.nn as nn

from re_agent.immuno.config import TrainConfig


def sigmoid_rampup(current: float, length: float) -> float:
    """Gaussian ramp from 0 to 1, as in the original Mean Teacher schedule."""
    if length <= 0:
        return 1.0
    phase = 1.0 - max(0.0, min(current, length)) / length
    return float(math.exp(-5.0 * phase * phase))


def build_teacher(student: nn.Module) -> nn.Module:
    teacher = copy.deepcopy(student)
    for param in teacher.parameters():
        param.requires_grad_(False)
    return teacher


@torch.no_grad()
def update_teacher(student: nn.Module, teacher: nn.Module, decay: float, step: int) -> None:
    """EMA step. Early on, `step` caps the decay so the teacher is not stuck at init."""
    alpha = min(1.0 - 1.0 / (step + 1), decay)
    for t_param, s_param in zip(teacher.parameters(), student.parameters()):
        t_param.mul_(alpha).add_(s_param.detach(), alpha=1.0 - alpha)
    for t_buf, s_buf in zip(teacher.buffers(), student.buffers()):
        t_buf.copy_(s_buf)


def augment(
    x: torch.Tensor,
    mask: torch.Tensor,
    cfg: TrainConfig,
    generator: torch.Generator | None = None,
) -> torch.Tensor:
    """Label-preserving perturbations of a frozen embedding window.

    Noise is scaled by the embedding's own spread so the strength does not depend
    on which ESM checkpoint produced it.
    """
    out = x
    if cfg.input_noise > 0:
        scale = cfg.input_noise * out.detach().std()
        noise = torch.empty_like(out).normal_(0.0, 1.0, generator=generator)
        out = out + noise * scale
    if cfg.feature_dropout > 0:
        keep = torch.rand(out.shape[0], 1, out.shape[2], device=out.device, generator=generator)
        out = out * (keep > cfg.feature_dropout).to(out.dtype)
    if cfg.residue_mask_prob > 0:
        keep = torch.rand(out.shape[0], out.shape[1], 1, device=out.device, generator=generator)
        out = out * (keep > cfg.residue_mask_prob).to(out.dtype)
    return out * mask.unsqueeze(-1)


def consistency_loss(student_logit: torch.Tensor, teacher_logit: torch.Tensor) -> torch.Tensor:
    """MSE between probabilities — bounded, so a confused teacher cannot dominate."""
    student_prob = torch.sigmoid(student_logit)
    teacher_prob = torch.sigmoid(teacher_logit.detach())
    return ((student_prob - teacher_prob) ** 2).mean()


def consistency_weight(epoch: float, cfg: TrainConfig) -> float:
    return cfg.consistency_weight * sigmoid_rampup(epoch, cfg.rampup_epochs)


def ema_decay_for(epoch: float, cfg: TrainConfig) -> float:
    """Track the student closely while it is still improving fast, then slow down."""
    return cfg.ema_decay_rampup if epoch < cfg.rampup_epochs else cfg.ema_decay_final
