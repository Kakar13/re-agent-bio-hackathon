"""Summary figures for the writeup, built from results/reports/validation.json."""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from re_agent.immuno.config import PATHS, ensure_dirs  # noqa: E402


def load_validation() -> dict:
    path = PATHS.reports / "validation.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} missing — run re_agent.immuno.validate first")
    return json.loads(path.read_text())


def plot_label_fraction(report: dict, path=None):
    """The controlled proof: SSL gain as a function of how many labels are hidden."""
    runs = report.get("label_fraction", {}).get("runs", [])
    if not runs:
        return None
    df = pd.DataFrame(runs)
    grouped = df.groupby("fraction").agg(["mean", "std"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    x = grouped.index.to_numpy()
    for col, label, color in (
        ("baseline_auroc", "labeled-only baseline", "#4C72B0"),
        ("mean_teacher_auroc", "mean teacher", "#C44E52"),
    ):
        mean = grouped[(col, "mean")].to_numpy()
        std = np.nan_to_num(grouped[(col, "std")].to_numpy())
        axes[0].plot(x, mean, "o-", color=color, label=label)
        axes[0].fill_between(x, mean - std, mean + std, color=color, alpha=0.18)
    axes[0].set_xlabel("fraction of IEDB labels kept")
    axes[0].set_ylabel("held-out AUROC")
    axes[0].set_title("Semi-supervised gain vs label budget")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    delta_mean = grouped[("delta", "mean")].to_numpy()
    delta_std = np.nan_to_num(grouped[("delta", "std")].to_numpy())
    axes[1].bar([str(v) for v in x], delta_mean, yerr=delta_std, color="#55A868", capsize=4)
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set_xlabel("fraction of labels kept")
    axes[1].set_ylabel("AUROC delta (MT - baseline)")
    axes[1].set_title("Improvement from unlabeled data")
    axes[1].grid(alpha=0.3, axis="y")

    fig.tight_layout()
    path = path or PATHS.figures / "label_fraction.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_denovo_gap(report: dict, path=None):
    """Why adaptation is needed at all: de novo windows sit outside training space."""
    behavior = report.get("denovo_behavior")
    if not behavior:
        return None
    gap = behavior["ood_gap"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].bar(
        ["natural\n(held-out IEDB)", "de novo\n(designed)"],
        [gap["natural_test_mean_knn_distance"], gap["denovo_mean_knn_distance"]],
        color=["#4C72B0", "#C44E52"],
    )
    axes[0].set_ylabel("mean kNN cosine distance to training peptides")
    axes[0].set_title(f"Distribution gap ({gap['ratio']:.2f}x further)")
    axes[0].grid(alpha=0.3, axis="y")

    metrics = ["perturbation_std", "mean_entropy", "frac_overconfident"]
    labels = [
        "prediction spread\nunder perturbation",
        "predictive\nentropy",
        "fraction\noverconfident",
    ]
    width = 0.36
    pos = np.arange(len(metrics))
    for offset, (name, color) in enumerate(
        (("baseline", "#4C72B0"), ("mean_teacher", "#C44E52"))
    ):
        values = [behavior[name][m] for m in metrics]
        axes[1].bar(pos + offset * width, values, width, label=name.replace("_", " "), color=color)
    axes[1].set_xticks(pos + width / 2)
    axes[1].set_xticklabels(labels, fontsize=8)
    axes[1].set_title("Behaviour on unlabeled de novo windows")
    axes[1].legend()
    axes[1].grid(alpha=0.3, axis="y")

    fig.tight_layout()
    path = path or PATHS.figures / "denovo_gap.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_epitope_recovery(report: dict, path=None):
    """Per-antigen epitope localization, baseline vs mean teacher."""
    rows = report.get("epitope_recovery", {}).get("per_antigen", [])
    if not rows:
        return None
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(df["baseline_auc"], df["mean_teacher_auc"], s=28, alpha=0.75, color="#4C72B0")
    lims = [
        min(df["baseline_auc"].min(), df["mean_teacher_auc"].min()) - 0.05,
        max(df["baseline_auc"].max(), df["mean_teacher_auc"].max()) + 0.05,
    ]
    ax.plot(lims, lims, "k--", lw=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("baseline AUC")
    ax.set_ylabel("mean teacher AUC")
    wins = int((df["mean_teacher_auc"] > df["baseline_auc"]).sum())
    ax.set_title(f"Epitope recovery per held-out antigen\n(MT better on {wins}/{len(df)})")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    path = path or PATHS.figures / "epitope_recovery.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> None:
    ensure_dirs()
    report = load_validation()
    for fn in (plot_label_fraction, plot_denovo_gap, plot_epitope_recovery):
        out = fn(report)
        print(f"{'wrote' if out else 'skipped'} {out or fn.__name__}")


if __name__ == "__main__":
    main()
