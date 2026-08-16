"""Leakage-safe IEDB dataset construction and calibrated MHC-I risk head."""

from __future__ import annotations

import json
import logging
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from re_agent.immuno_risk.peptides import peptide_descriptors
from re_agent.immuno_risk.reference_data import (
    ROOT,
    ensure_fixtures,
    load_overlap_8mers,
    load_tcell_rows,
)

log = logging.getLogger(__name__)

MODEL_DIR = ROOT / "data" / "processed" / "immuno" / "models"


def _normalize_allele(a: str) -> str:
    a = a.strip().replace("HLA-", "")
    if a.startswith("A") or a.startswith("B") or a.startswith("C"):
        if "*" not in a and len(a) >= 3:
            # A0201 -> A*02:01
            locus, rest = a[0], a[1:]
            if len(rest) >= 4 and rest.isdigit():
                return f"HLA-{locus}*{rest[:2]}:{rest[2:4]}"
        return f"HLA-{a}" if not a.startswith("HLA-") else a
    return a if a.startswith("HLA-") else f"HLA-{a}"


def _is_positive(outcome: str) -> bool | None:
    o = outcome.strip().lower()
    if o.startswith("pos") or o in {"positive", "positive-high", "positive-intermediate", "positive-low"}:
        return True
    if o.startswith("neg") or o == "negative":
        return False
    return None


def _has_8mer_overlap(peptide: str, ledger: set[str]) -> bool:
    p = peptide.upper()
    if len(p) < 8:
        return p in ledger
    for i in range(len(p) - 7):
        if p[i : i + 8] in ledger:
            return True
    return False


@dataclass
class LabeledExample:
    peptide: str
    allele: str
    label: int
    source_protein: str
    reference_id: str
    publication_year: str
    overlap_flag: bool


def build_iedb_examples(path: Path | None = None) -> list[LabeledExample]:
    ensure_fixtures()
    ledger = load_overlap_8mers()
    rows = load_tcell_rows(path)
    examples: list[LabeledExample] = []
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        peptide = row["peptide"].upper().strip()
        if not (8 <= len(peptide) <= 11) or not peptide.isalpha():
            continue
        allele = _normalize_allele(row["allele"])
        if not any(x in allele for x in ("A*", "B*", "C*")):
            continue
        label = _is_positive(row.get("qualitative_outcome", ""))
        if label is None:
            continue  # never invent negatives from untested
        key = (peptide, allele, row.get("reference_id", ""))
        if key in seen:
            continue
        seen.add(key)
        examples.append(
            LabeledExample(
                peptide=peptide,
                allele=allele,
                label=int(label),
                source_protein=row.get("source_protein", ""),
                reference_id=row.get("reference_id", ""),
                publication_year=row.get("publication_year", ""),
                overlap_flag=_has_8mer_overlap(peptide, ledger),
            )
        )
    return examples


def grouped_split(
    examples: list[LabeledExample],
    *,
    val_frac: float = 0.2,
    test_frac: float = 0.2,
) -> dict[str, list[LabeledExample]]:
    """Split by publication + source protein — not random rows."""
    groups: dict[str, list[LabeledExample]] = defaultdict(list)
    for ex in examples:
        gid = f"{ex.reference_id}|{ex.source_protein}"
        groups[gid].append(ex)
    keys = sorted(groups.keys())
    n = len(keys)
    n_test = max(1, int(n * test_frac)) if n > 2 else (1 if n > 1 else 0)
    n_val = max(1, int(n * val_frac)) if n > 3 else (1 if n > 2 else 0)
    test_keys = set(keys[:n_test])
    val_keys = set(keys[n_test : n_test + n_val])
    train_keys = set(keys) - test_keys - val_keys
    if not train_keys and keys:
        train_keys = {keys[-1]}
        test_keys -= train_keys
        val_keys -= train_keys
    return {
        "train": [e for k in train_keys for e in groups[k]],
        "val": [e for k in val_keys for e in groups[k]],
        "test": [e for k in test_keys for e in groups[k]],
    }


def _features_from_mhc(
    peptide: str,
    mhc: dict[str, Any] | None,
) -> list[float]:
    desc = peptide_descriptors(peptide)
    affinity = float(mhc.get("affinity_nm") or 5000.0) if mhc else 5000.0
    presentation = float(mhc.get("presentation_score") or 0.0) if mhc else 0.0
    processing = float(mhc.get("processing_score") or 0.0) if mhc else 0.0
    pct = float(mhc.get("percentile_rank") if mhc and mhc.get("percentile_rank") is not None else 50.0)
    return [
        math.log10(max(affinity, 1e-3)),
        presentation,
        processing,
        pct / 100.0,
        desc["hydrophobic_fraction"],
        desc["charge_proxy"],
        desc["aromatic_fraction"],
        desc["cysteine_fraction"],
        desc["length"] / 11.0,
    ]


class LogisticRiskModel:
    """Tiny inspectable logistic regression with optional Platt-like temperature."""

    def __init__(self, weights: list[float], bias: float = 0.0, temperature: float = 1.0):
        self.weights = weights
        self.bias = bias
        self.temperature = temperature

    def predict_proba(self, x: list[float]) -> float:
        z = self.bias + sum(w * xi for w, xi in zip(self.weights, x, strict=False))
        z /= max(self.temperature, 1e-6)
        return 1.0 / (1.0 + math.exp(-max(min(z, 40), -40)))

    def to_dict(self) -> dict[str, Any]:
        return {
            "type": "logistic_risk_v1",
            "weights": self.weights,
            "bias": self.bias,
            "temperature": self.temperature,
            "feature_names": [
                "log10_affinity",
                "presentation",
                "processing",
                "percentile_frac",
                "hydrophobic_fraction",
                "charge_proxy",
                "aromatic_fraction",
                "cysteine_fraction",
                "length_frac",
            ],
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> LogisticRiskModel:
        return cls(list(d["weights"]), float(d["bias"]), float(d.get("temperature", 1.0)))


def _sigmoid(z: float) -> float:
    return 1.0 / (1.0 + math.exp(-max(min(z, 40), -40)))


def train_logistic(
    X: list[list[float]],
    y: list[int],
    *,
    lr: float = 0.05,
    epochs: int = 400,
    l2: float = 0.01,
) -> LogisticRiskModel:
    if not X:
        # Uninformative prior
        return LogisticRiskModel([0.0] * 9, bias=-1.0)
    dim = len(X[0])
    w = [0.0] * dim
    b = 0.0
    n = len(X)
    for _ in range(epochs):
        gw = [0.0] * dim
        gb = 0.0
        for xi, yi in zip(X, y, strict=True):
            p = _sigmoid(b + sum(wj * xj for wj, xj in zip(w, xi, strict=True)))
            err = p - yi
            for j in range(dim):
                gw[j] += err * xi[j]
            gb += err
        for j in range(dim):
            w[j] -= lr * ((gw[j] / n) + l2 * w[j])
        b -= lr * (gb / n)
    return LogisticRiskModel(w, b)


def calibrate_temperature(
    model: LogisticRiskModel,
    X: list[list[float]],
    y: list[int],
) -> LogisticRiskModel:
    if not X:
        return model
    best_t, best_nll = 1.0, float("inf")
    for t in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]:
        nll = 0.0
        for xi, yi in zip(X, y, strict=True):
            p = LogisticRiskModel(model.weights, model.bias, t).predict_proba(xi)
            p = min(max(p, 1e-6), 1 - 1e-6)
            nll += -(yi * math.log(p) + (1 - yi) * math.log(1 - p))
        if nll < best_nll:
            best_nll, best_t = nll, t
    return LogisticRiskModel(model.weights, model.bias, best_t)


def _mhc_features_for_examples(
    examples: list[LabeledExample],
) -> tuple[list[list[float]], list[int]]:
    X: list[list[float]] = []
    y: list[int] = []
    try:
        from re_agent.immuno_risk.mhcflurry_backend import score_peptide_allele_pairs

        pairs = [(e.peptide, e.allele) for e in examples]
        scored = score_peptide_allele_pairs(pairs)
        by_key = {(s["peptide"], s["allele"]): s for s in scored}
    except Exception as exc:  # noqa: BLE001
        log.warning("MHCflurry features unavailable (%s); using descriptors only", exc)
        by_key = {}

    for e in examples:
        mhc = by_key.get((e.peptide, e.allele))
        X.append(_features_from_mhc(e.peptide, mhc))
        y.append(e.label)
    return X, y


def train_and_save(path: Path | None = None) -> Path:
    examples = build_iedb_examples(path)
    splits = grouped_split(examples)
    X_train, y_train = _mhc_features_for_examples(splits["train"] or examples)
    X_val, y_val = _mhc_features_for_examples(splits["val"] or splits["train"] or examples)
    model = train_logistic(X_train, y_train)
    model = calibrate_temperature(model, X_val, y_val)

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    out = MODEL_DIR / "iedb_logistic_v1.json"
    payload = {
        "model": model.to_dict(),
        "n_train": len(X_train),
        "n_val": len(X_val),
        "n_test": len(splits["test"]),
        "n_overlap_flagged": sum(1 for e in examples if e.overlap_flag),
        "split": "grouped_by_reference_and_source_protein",
        "label_rule": "explicit IEDB qualitative T-cell outcomes only",
    }
    out.write_text(json.dumps(payload, indent=2))
    return out


def load_model(path: Path | None = None) -> LogisticRiskModel:
    path = path or (MODEL_DIR / "iedb_logistic_v1.json")
    if not path.exists():
        # Train from fixtures on first use
        path = train_and_save()
    data = json.loads(path.read_text())
    return LogisticRiskModel.from_dict(data["model"])


def score_iedb_risk(
    peptide: str,
    *,
    affinity_nm: float | None,
    presentation_score: float | None,
    processing_score: float | None,
    percentile_rank: float | None,
    model: LogisticRiskModel | None = None,
) -> float:
    model = model or load_model()
    mhc = {
        "affinity_nm": affinity_nm,
        "presentation_score": presentation_score,
        "processing_score": processing_score,
        "percentile_rank": percentile_rank,
    }
    return model.predict_proba(_features_from_mhc(peptide, mhc))


def baseline_from_presentation(
    presentation_score: float | None,
    percentile_rank: float | None,
) -> float:
    """Fixed MHCflurry baseline in [0,1] — not blended with IEDB head."""
    if presentation_score is not None:
        return float(min(max(presentation_score, 0.0), 1.0))
    if percentile_rank is not None:
        return float(min(max(1.0 - percentile_rank / 100.0, 0.0), 1.0))
    return 0.0
