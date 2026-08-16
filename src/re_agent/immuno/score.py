"""Window risks -> distinct epitope regions -> protein risk and percentile.

Naively combining stride-1 windows with a noisy-or saturates every protein at 1.0,
because neighbouring windows are near-duplicates rather than independent tests.
Overlapping windows are first collapsed to disjoint candidate epitope regions, and
only those are combined.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import torch

from re_agent.immuno.confidence import OODIndex, aggregate_confidence, window_confidence
from re_agent.immuno.config import PATHS, WINDOW, TrainConfig, ensure_dirs, torch_device
from re_agent.immuno.data import windows
from re_agent.immuno.embed import embed_sequences, load_embeddings
from re_agent.immuno.model import ImmunoHead


@dataclass
class Region:
    start: int
    end: int
    risk: float
    peptide: str
    confidence: float = 0.0
    stability: float = 0.0
    agreement: float = 0.0
    familiarity: float = 0.0


@dataclass
class ProteinScore:
    name: str
    sequence: str
    risk: float
    risk_percentile: float
    # Length-independent companion to `risk`: the single most dangerous window.
    peak_window_risk: float
    confidence: float
    confidence_parts: dict
    regions: list[Region]
    per_residue: np.ndarray = field(repr=False)
    per_residue_saliency: np.ndarray = field(repr=False)
    window_risks: np.ndarray = field(repr=False)
    window_starts: np.ndarray = field(repr=False)


def select_regions(
    risks: np.ndarray, starts: np.ndarray, width: int = WINDOW, max_regions: int = 25
) -> list[int]:
    """Greedy non-maximum suppression over window spans.

    Returns indices of disjoint windows, highest risk first, so each retained
    window stands for one distinct candidate epitope.
    """
    order = np.argsort(-risks)
    chosen: list[int] = []
    claimed: list[tuple[int, int]] = []
    for i in order:
        s, e = int(starts[i]), int(starts[i]) + width
        if any(s < ce and cs < e for cs, ce in claimed):
            continue
        chosen.append(int(i))
        claimed.append((s, e))
        if len(chosen) >= max_regions:
            break
    return chosen


def noisy_or(risks: np.ndarray) -> float:
    """Probability that at least one region drives a response."""
    if len(risks) == 0:
        return 0.0
    return float(1.0 - np.prod(1.0 - np.clip(risks, 0.0, 1.0 - 1e-9)))


def protein_risk(
    risks: np.ndarray, starts: np.ndarray, width: int = WINDOW
) -> tuple[float, list[int]]:
    idx = select_regions(risks, starts, width)
    return noisy_or(risks[idx]), idx


def percentile_of(
    value: float,
    reference: np.ndarray,
    length: int | None = None,
    reference_lengths: np.ndarray | None = None,
    tolerance: float = 0.5,
    min_matched: int = 30,
) -> float:
    """Percentile against the natural cohort, matched on length where possible.

    A noisy-or over regions grows with sequence length, so a 90-residue binder
    compared against 400-residue proteins would look artificially safe. Restricting
    the comparison to similar-length natural proteins removes that bias.
    """
    if len(reference) == 0:
        return float("nan")
    pool = reference
    has_lengths = reference_lengths is not None and len(reference_lengths) == len(reference)
    if length is not None and has_lengths:
        low, high = length * (1 - tolerance), length * (1 + tolerance)
        matched = reference[(reference_lengths >= low) & (reference_lengths <= high)]
        if len(matched) >= min_matched:
            pool = matched
    return float(100.0 * np.searchsorted(np.sort(pool), value) / len(pool))


# --------------------------------------------------------------------------- #
# Runtime bundle
# --------------------------------------------------------------------------- #


@dataclass
class Bundle:
    student: ImmunoHead
    teacher: ImmunoHead
    ood: OODIndex | None
    device: str
    reference_risks: np.ndarray
    reference_lengths: np.ndarray = field(default_factory=lambda: np.array([]))


def load_bundle(
    name: str = "mean_teacher", with_ood: bool = True, device: str | None = None
) -> Bundle:
    device = device or torch_device()
    ckpt = torch.load(PATHS.models / f"{name}.pt", map_location=device, weights_only=False)
    cfg = TrainConfig(**ckpt["config"])
    student = ImmunoHead(cfg).to(device)
    teacher = ImmunoHead(cfg).to(device)
    student.load_state_dict(ckpt["student"])
    teacher.load_state_dict(ckpt["teacher"])
    student.eval()
    teacher.eval()

    ood = None
    if with_ood:
        from re_agent.immuno.train import prepare_labeled

        labeled, _ = prepare_labeled(seed=cfg.seed)
        emb, mask = load_embeddings("labeled")
        rows = labeled[labeled.split == "train"]["row"].to_numpy()
        ood = OODIndex(emb[rows], mask[rows])

    ref_path = PATHS.models / f"{name}_reference_risks.npz"
    if ref_path.exists():
        cached = np.load(ref_path)
        reference, lengths = cached["risks"], cached["lengths"]
    else:
        reference, lengths = np.array([]), np.array([])
    return Bundle(student, teacher, ood, device, reference, lengths)


@torch.no_grad()
def _window_probs(
    model, emb: np.ndarray, mask: np.ndarray, device: str, batch: int = 4096
) -> np.ndarray:
    out = []
    for i in range(0, len(emb), batch):
        x = torch.from_numpy(np.array(emb[i : i + batch], dtype=np.float32)).to(device)
        m = torch.from_numpy(np.array(mask[i : i + batch], dtype=np.float32)).to(device)
        out.append(torch.sigmoid(model(x, m)).cpu().numpy())
    return np.concatenate(out) if out else np.array([])


def build_reference_risks(bundle: Bundle, name: str = "mean_teacher") -> np.ndarray:
    """Protein-level risk for every natural reference protein, for the percentile scale."""
    from re_agent.immuno.data import build_reference

    ensure_dirs()
    ref = build_reference()
    emb, mask = load_embeddings("reference")
    probs = _window_probs(bundle.teacher, emb, mask, bundle.device)
    ref = ref.copy()
    ref["risk"] = probs

    risks, lengths = [], []
    for _, grp in ref.groupby("parent", sort=False):
        starts = grp["start"].to_numpy()
        r, _ = protein_risk(grp["risk"].to_numpy(), starts)
        risks.append(r)
        lengths.append(int(starts.max()) + WINDOW)
    arr, lens = np.array(risks), np.array(lengths)
    np.savez(PATHS.models / f"{name}_reference_risks.npz", risks=arr, lengths=lens)
    bundle.reference_risks = arr
    bundle.reference_lengths = lens
    return arr


def score_sequence(
    sequence: str,
    bundle: Bundle,
    name: str = "query",
    stride: int = 1,
    mc_passes: int = 20,
) -> ProteinScore:
    """Full assessment of one protein: risk, percentile, confidence, heatmap."""
    from re_agent.immuno.explain import per_residue_heatmap, per_residue_saliency

    tiles = windows(sequence, stride=stride)
    starts = np.array([s for s, _ in tiles])
    seqs = [w for _, w in tiles]
    emb, mask = embed_sequences(seqs, device=bundle.device)

    x = torch.from_numpy(emb.astype(np.float32)).to(bundle.device)
    m = torch.from_numpy(mask).to(bundle.device)

    comps = window_confidence(
        bundle.student,
        bundle.teacher,
        x,
        m,
        ood=bundle.ood,
        emb_np=emb,
        mask_np=mask,
        passes=mc_passes,
    )
    risks = comps["risk"]

    total_risk, region_idx = protein_risk(risks, starts)
    regions = [
        Region(
            start=int(starts[i]),
            end=int(starts[i]) + len(seqs[i]),
            risk=float(risks[i]),
            peptide=seqs[i],
            confidence=float(comps["confidence"][i]),
            stability=float(comps["stability"][i]),
            agreement=float(comps["agreement"][i]),
            familiarity=float(comps["familiarity"][i]),
        )
        for i in region_idx
    ]

    conf_parts = aggregate_confidence(comps, weights=risks)
    heat = per_residue_heatmap(bundle.teacher, x, m, risks, starts, len(sequence))
    sal = per_residue_saliency(bundle.teacher, x, m, starts, len(sequence))

    return ProteinScore(
        name=name,
        sequence=sequence,
        risk=total_risk,
        risk_percentile=percentile_of(
            total_risk,
            bundle.reference_risks,
            length=len(sequence),
            reference_lengths=bundle.reference_lengths,
        ),
        peak_window_risk=float(risks.max()) if len(risks) else 0.0,
        confidence=conf_parts["confidence"],
        confidence_parts=conf_parts,
        regions=regions,
        per_residue=heat,
        per_residue_saliency=sal,
        window_risks=risks,
        window_starts=starts,
    )


def score_to_dict(score: ProteinScore) -> dict:
    return {
        "name": score.name,
        "length": len(score.sequence),
        "risk": round(score.risk, 4),
        "risk_percentile_vs_natural": round(score.risk_percentile, 1),
        "peak_window_risk": round(score.peak_window_risk, 4),
        "confidence": round(score.confidence, 4),
        "confidence_breakdown": {k: round(v, 4) for k, v in score.confidence_parts.items()},
        "n_regions": len(score.regions),
        "top_regions": [
            {
                "start": r.start + 1,  # 1-indexed for biologists
                "end": r.end,
                "peptide": r.peptide,
                "risk": round(r.risk, 4),
                "confidence": round(r.confidence, 4),
                "familiarity": round(r.familiarity, 4),
            }
            for r in score.regions[:10]
        ],
        "per_residue_risk": [round(float(v), 4) for v in score.per_residue],
    }


def write_report(score: ProteinScore, path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or PATHS.reports / f"{score.name}.json"
    path.write_text(json.dumps(score_to_dict(score), indent=2))
    return path
