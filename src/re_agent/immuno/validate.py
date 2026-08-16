"""Experiments that decide whether the unlabeled de novo data actually helped.

There are no de novo immunogenicity labels, so the claim is supported three ways:
a controlled low-label experiment where the answer is known, an organism-holdout
that simulates the natural-to-designed shift with labels on the far side, and
descriptive robustness/OOD measurements on the real de novo pool.

Every comparison holds the architecture, seed, step count and augmentation fixed
so the only moving part is the consistency term over unlabeled windows.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

from re_agent.immuno import mean_teacher as mt
from re_agent.immuno.confidence import OODIndex
from re_agent.immuno.config import PATHS, TrainConfig, ensure_dirs, torch_device
from re_agent.immuno.embed import load_embeddings
from re_agent.immuno.train import prepare_labeled, train_model


def _arm(
    cfg: TrainConfig,
    labeled: pd.DataFrame,
    emb_l,
    mask_l,
    use_unlabeled: bool,
    emb_u=None,
    mask_u=None,
    device: str = "cpu",
) -> dict:
    result = train_model(
        cfg, labeled, emb_l, mask_l, use_unlabeled, emb_u, mask_u, device=device, verbose=False
    )
    result.pop("_models", None)
    return result


class MainPair:
    """The full-data baseline/Mean Teacher pair, trained once and shared.

    Three of the five experiments need exactly this pair; retraining it each time
    would triple the runtime for identical models.
    """

    def __init__(self, cfg, labeled, emb_l, mask_l, emb_u, mask_u, device):
        self._args = (cfg, labeled, emb_l, mask_l, emb_u, mask_u, device)
        self._pair = None

    def get(self) -> tuple[dict, dict]:
        if self._pair is None:
            cfg, labeled, emb_l, mask_l, emb_u, mask_u, device = self._args
            print("  training shared baseline / mean-teacher pair...")
            base = train_model(cfg, labeled, emb_l, mask_l, False, device=device, verbose=False)
            ssl = train_model(
                cfg, labeled, emb_l, mask_l, True, emb_u, mask_u, device=device, verbose=False
            )
            self._pair = (base, ssl)
        return self._pair


def experiment_in_distribution(pair: MainPair) -> dict:
    """Does adding unlabeled de novo data damage natural-peptide performance?"""
    print("\n[1] in-distribution preservation (unlabeled = real de novo windows)")
    base, ssl = pair.get()
    out = {
        "baseline": base["test_teacher"],
        "mean_teacher": ssl["test_teacher"],
        "delta_auroc": ssl["test_teacher"]["auroc"] - base["test_teacher"]["auroc"],
    }
    print(
        f"    baseline AUROC {out['baseline']['auroc']:.4f} -> "
        f"mean teacher {out['mean_teacher']['auroc']:.4f} ({out['delta_auroc']:+.4f})"
    )
    return out


def experiment_label_fraction(
    cfg, labeled, emb_l, mask_l, device, fractions=(0.05, 0.1, 0.2), seeds=(0, 1, 2)
) -> dict:
    """The controlled proof: hide most labels and let SSL use them unlabeled.

    Here the ground truth is known, so any gain is measurable rather than argued.
    """
    print("\n[2] label-fraction experiment (unlabeled = IEDB train windows, labels discarded)")
    results = []
    train_mask = labeled.split == "train"

    for frac in fractions:
        for seed in seeds:
            rng = np.random.default_rng(1000 * seed + int(frac * 100))
            groups = labeled.loc[train_mask, "group"].unique()
            rng.shuffle(groups)
            keep = set(groups[: max(1, int(len(groups) * frac))])

            subset = labeled.copy()
            keep_rows = subset["group"].isin(keep) & train_mask
            # Everything not kept as labeled becomes the unlabeled pool.
            hidden_rows = train_mask & ~keep_rows
            subset.loc[train_mask & ~keep_rows, "split"] = "unused"

            hidden_index = labeled.loc[hidden_rows, "row"].to_numpy()
            emb_hidden = emb_l[hidden_index]
            mask_hidden = mask_l[hidden_index]

            arm_cfg = TrainConfig(**{**cfg.__dict__, "seed": seed})
            base = _arm(arm_cfg, subset, emb_l, mask_l, False, device=device)
            ssl = _arm(
                arm_cfg, subset, emb_l, mask_l, True, emb_hidden, mask_hidden, device=device
            )
            row = {
                "fraction": frac,
                "seed": seed,
                "n_labeled": int(keep_rows.sum()),
                "n_unlabeled": int(len(hidden_index)),
                "baseline_auroc": base["test_teacher"]["auroc"],
                "mean_teacher_auroc": ssl["test_teacher"]["auroc"],
            }
            row["delta"] = row["mean_teacher_auroc"] - row["baseline_auroc"]
            results.append(row)
            print(
                f"    frac {frac:>5.2f} seed {seed}: "
                f"labeled {row['n_labeled']:>6d} | baseline {row['baseline_auroc']:.4f} -> "
                f"MT {row['mean_teacher_auroc']:.4f} ({row['delta']:+.4f})"
            )

    df = pd.DataFrame(results)
    summary = (
        df.groupby("fraction")
        .agg(
            baseline=("baseline_auroc", "mean"),
            mean_teacher=("mean_teacher_auroc", "mean"),
            delta=("delta", "mean"),
            delta_std=("delta", "std"),
        )
        .reset_index()
    )
    print("\n    mean over seeds:")
    for r in summary.itertuples(index=False):
        print(
            f"      frac {r.fraction:>5.2f}: {r.baseline:.4f} -> {r.mean_teacher:.4f} "
            f"(delta {r.delta:+.4f} +/- {0.0 if pd.isna(r.delta_std) else r.delta_std:.4f})"
        )
    return {"runs": results, "summary": summary.to_dict(orient="records")}


def experiment_domain_shift(cfg, labeled, emb_l, mask_l, device, min_windows: int = 3000) -> dict:
    """Train on some source organisms, adapt to a held-out one with its labels hidden.

    This is the closest labeled analogue of the natural-to-de-novo shift we care
    about: the target domain's sequences are available, its labels are not.
    """
    print("\n[3] organism-holdout domain shift")
    counts = labeled["source_organism"].value_counts()
    candidates = [o for o in counts.index if counts[o] >= min_windows and o != "unknown"][:4]
    print(f"    candidate holdout organisms: {candidates}")

    results = []
    for organism in candidates:
        is_target = labeled["source_organism"] == organism
        if is_target.sum() < min_windows:
            continue
        target = labeled[is_target]
        # Half the target domain is the unlabeled adaptation pool, half is the test set.
        rng = np.random.default_rng(0)
        groups = target["group"].unique()
        rng.shuffle(groups)
        adapt_groups = set(groups[: len(groups) // 2])

        subset = labeled.copy()
        subset.loc[is_target, "split"] = "unused"
        eval_rows = is_target & ~labeled["group"].isin(adapt_groups)
        subset.loc[eval_rows, "split"] = "test"
        # Source domain keeps its own train/val, target contributes nothing labeled.
        subset.loc[~is_target & (labeled.split == "test"), "split"] = "train"

        adapt_index = labeled.loc[is_target & labeled["group"].isin(adapt_groups), "row"].to_numpy()
        if len(adapt_index) < 500 or int(eval_rows.sum()) < 500:
            continue
        y_eval = labeled.loc[eval_rows, "label"]
        if y_eval.nunique() < 2:
            continue

        base = _arm(cfg, subset, emb_l, mask_l, False, device=device)
        ssl = _arm(
            cfg, subset, emb_l, mask_l, True, emb_l[adapt_index], mask_l[adapt_index], device=device
        )
        row = {
            "organism": organism,
            "n_adapt_unlabeled": int(len(adapt_index)),
            "n_test": int(eval_rows.sum()),
            "baseline_auroc": base["test_teacher"]["auroc"],
            "mean_teacher_auroc": ssl["test_teacher"]["auroc"],
        }
        row["delta"] = row["mean_teacher_auroc"] - row["baseline_auroc"]
        results.append(row)
        print(
            f"    {organism[:38]:<38} baseline {row['baseline_auroc']:.4f} -> "
            f"MT {row['mean_teacher_auroc']:.4f} ({row['delta']:+.4f})"
        )
    return {"runs": results}


def experiment_denovo_behavior(cfg, labeled, emb_l, mask_l, emb_u, mask_u, device, pair) -> dict:
    """Descriptive: how the two models behave on real de novo windows.

    No labels exist here, so we report robustness under perturbation, predictive
    entropy, and how far de novo windows sit from the training distribution.
    """
    print("\n[4] de novo behaviour and OOD gap")
    base, ssl = pair.get()

    rng = np.random.default_rng(0)
    sample = np.sort(rng.choice(len(emb_u), size=min(8000, len(emb_u)), replace=False))
    xu = torch.from_numpy(np.asarray(emb_u[sample], dtype=np.float32)).to(device)
    mu = torch.from_numpy(np.asarray(mask_u[sample], dtype=np.float32)).to(device)

    out = {}
    for name, run in (("baseline", base), ("mean_teacher", ssl)):
        model = run["_models"]["teacher"]
        model.eval()
        with torch.no_grad():
            probs = torch.sigmoid(model(xu, mu))
            # Prediction spread under repeated perturbation = robustness on this domain.
            perturbed = torch.stack(
                [torch.sigmoid(model(mt.augment(xu, mu, cfg), mu)) for _ in range(10)]
            )
        entropy = -(
            probs * torch.log(probs.clamp_min(1e-9))
            + (1 - probs) * torch.log((1 - probs).clamp_min(1e-9))
        )
        out[name] = {
            "mean_risk": float(probs.mean()),
            "perturbation_std": float(perturbed.std(0).mean()),
            "mean_entropy": float(entropy.mean()),
            "frac_overconfident": float(((probs > 0.9) | (probs < 0.1)).float().mean()),
        }
        print(
            f"    {name:<13} perturbation std {out[name]['perturbation_std']:.4f} | "
            f"entropy {out[name]['mean_entropy']:.4f} | "
            f"overconfident {out[name]['frac_overconfident']:.3f}"
        )

    train_rows = labeled[labeled.split == "train"]["row"].to_numpy()
    ood = OODIndex(emb_l[train_rows], mask_l[train_rows])
    test_rows = labeled[labeled.split == "test"]["row"].to_numpy()
    nat = ood.distance(np.asarray(emb_l[test_rows[:4000]]), np.asarray(mask_l[test_rows[:4000]]))
    dn = ood.distance(np.asarray(emb_u[sample[:4000]]), np.asarray(mask_u[sample[:4000]]))
    out["ood_gap"] = {
        "natural_test_mean_knn_distance": float(nat.mean()),
        "denovo_mean_knn_distance": float(dn.mean()),
        "ratio": float(dn.mean() / max(nat.mean(), 1e-9)),
    }
    print(
        f"    kNN distance to training peptides: natural {nat.mean():.4f} vs "
        f"de novo {dn.mean():.4f} ({out['ood_gap']['ratio']:.2f}x)"
    )
    return out


def experiment_epitope_recovery(
    labeled, emb_l, mask_l, device, pair, min_peptides: int = 25
) -> dict:
    """Does the per-window risk localize real epitopes within held-out antigens?

    For each test-split antigen with enough mapped peptides, rank its peptides by
    predicted risk and score against the measured positives. This is the heatmap
    sanity check, run across many antigens rather than a single case study.
    """
    print("\n[5] epitope recovery on held-out antigens")
    base, ssl = pair.get()

    test = labeled[labeled.split == "test"]
    per_antigen = []
    for antigen, grp in test.groupby("group"):
        if len(grp) < min_peptides or grp["label"].nunique() < 2:
            continue
        rows = grp["row"].to_numpy()
        y = grp["label"].to_numpy()
        x = torch.from_numpy(np.asarray(emb_l[rows], dtype=np.float32)).to(device)
        m = torch.from_numpy(np.asarray(mask_l[rows], dtype=np.float32)).to(device)
        entry = {"antigen": str(antigen)[:60], "n_peptides": len(grp), "n_positive": int(y.sum())}
        for name, run in (("baseline", base), ("mean_teacher", ssl)):
            with torch.no_grad():
                p = torch.sigmoid(run["_models"]["teacher"](x, m)).cpu().numpy()
            entry[f"{name}_auc"] = float(roc_auc_score(y, p))
        per_antigen.append(entry)

    df = pd.DataFrame(per_antigen)
    summary = {}
    if len(df):
        summary = {
            "n_antigens": len(df),
            "baseline_mean_auc": float(df["baseline_auc"].mean()),
            "mean_teacher_mean_auc": float(df["mean_teacher_auc"].mean()),
            "mean_teacher_wins": int((df["mean_teacher_auc"] > df["baseline_auc"]).sum()),
        }
        print(
            f"    {summary['n_antigens']} antigens | baseline {summary['baseline_mean_auc']:.4f} "
            f"-> MT {summary['mean_teacher_mean_auc']:.4f} | "
            f"MT wins on {summary['mean_teacher_wins']}/{summary['n_antigens']}"
        )
    return {"per_antigen": per_antigen, "summary": summary}


def main() -> None:
    parser = argparse.ArgumentParser(description="Validation experiments")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument(
        "--only",
        nargs="*",
        choices=["indist", "fraction", "shift", "denovo", "epitope"],
        default=None,
    )
    args = parser.parse_args()

    ensure_dirs()
    cfg = TrainConfig(epochs=args.epochs)
    device = torch_device()
    labeled, _ = prepare_labeled(seed=cfg.seed)
    emb_l, mask_l = load_embeddings("labeled")
    emb_u, mask_u = load_embeddings("unlabeled")
    print(f"device={device} labeled={len(labeled)} unlabeled={len(emb_u)}")

    want = set(args.only) if args.only else {"indist", "fraction", "shift", "denovo", "epitope"}
    pair = MainPair(cfg, labeled, emb_l, mask_l, emb_u, mask_u, device)
    report = {}
    if "indist" in want:
        report["in_distribution"] = experiment_in_distribution(pair)
    if "fraction" in want:
        report["label_fraction"] = experiment_label_fraction(cfg, labeled, emb_l, mask_l, device)
    if "shift" in want:
        report["domain_shift"] = experiment_domain_shift(cfg, labeled, emb_l, mask_l, device)
    if "denovo" in want:
        report["denovo_behavior"] = experiment_denovo_behavior(
            cfg, labeled, emb_l, mask_l, emb_u, mask_u, device, pair
        )
    if "epitope" in want:
        report["epitope_recovery"] = experiment_epitope_recovery(
            labeled, emb_l, mask_l, device, pair
        )

    path = PATHS.reports / "validation.json"
    path.write_text(json.dumps(report, indent=2, default=str))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
