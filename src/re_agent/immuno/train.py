"""Training for the labeled-only baseline and the Mean Teacher model.

Both arms keep an EMA teacher and both see the same augmented labeled batches, so
the only difference between them is the consistency term over unlabeled de novo
windows. Any gap is attributable to the unlabeled data, not to weight averaging.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score

from re_agent.immuno import mean_teacher as mt
from re_agent.immuno.config import PATHS, TrainConfig, ensure_dirs, torch_device
from re_agent.immuno.data import build_labeled, group_split
from re_agent.immuno.embed import load_embeddings
from re_agent.immuno.model import ImmunoHead, bce_loss, count_parameters


def prepare_labeled(seed: int = 0) -> tuple[pd.DataFrame, np.ndarray]:
    """Deduplicate windows by sequence and assign group-wise splits.

    The same 15-mer can be tiled out of several parent peptides; collapsing them
    first stops an identical window from landing in both train and test.
    """
    df = build_labeled().reset_index(drop=True)
    df["row"] = np.arange(len(df))
    agg = (
        df.groupby("seq", sort=False)
        .agg(
            label=("label", "max"),
            group=("group", "first"),
            row=("row", "first"),
            source_organism=("source_organism", "first"),
            parent=("parent", "first"),
        )
        .reset_index()
    )
    agg["split"] = group_split(agg["group"], seed=seed)
    return agg, agg["row"].to_numpy()


def _to_tensor(emb: np.ndarray, mask: np.ndarray, rows: np.ndarray, device: str):
    x = torch.from_numpy(np.asarray(emb[rows], dtype=np.float32)).to(device)
    m = torch.from_numpy(np.asarray(mask[rows], dtype=np.float32)).to(device)
    return x, m


def calibration_error(probs: np.ndarray, labels: np.ndarray, bins: int = 15) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    idx = np.clip(np.digitize(probs, edges) - 1, 0, bins - 1)
    ece = 0.0
    for b in range(bins):
        sel = idx == b
        if not sel.any():
            continue
        ece += sel.mean() * abs(probs[sel].mean() - labels[sel].mean())
    return float(ece)


@torch.no_grad()
def evaluate(
    model, emb, mask, rows: np.ndarray, labels: np.ndarray, device: str, batch: int = 4096
) -> dict:
    model.eval()
    probs = []
    for i in range(0, len(rows), batch):
        x, m = _to_tensor(emb, mask, rows[i : i + batch], device)
        probs.append(torch.sigmoid(model(x, m)).cpu().numpy())
    probs = np.concatenate(probs)
    return {
        "auroc": float(roc_auc_score(labels, probs)),
        "auprc": float(average_precision_score(labels, probs)),
        "accuracy": float(((probs > 0.5).astype(int) == labels).mean()),
        "brier": float(np.mean((probs - labels) ** 2)),
        "ece": calibration_error(probs, labels),
        "n": int(len(labels)),
    }


def train_model(
    cfg: TrainConfig,
    labeled: pd.DataFrame,
    emb_l: np.ndarray,
    mask_l: np.ndarray,
    use_unlabeled: bool,
    emb_u: np.ndarray | None = None,
    mask_u: np.ndarray | None = None,
    device: str | None = None,
    log_every: int = 5,
    verbose: bool = True,
) -> dict:
    device = device or torch_device()
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    train_df = labeled[labeled.split == "train"]
    val_df = labeled[labeled.split == "val"]
    test_df = labeled[labeled.split == "test"]
    train_rows = train_df["row"].to_numpy()
    train_y = train_df["label"].to_numpy().astype(np.float32)

    student = ImmunoHead(cfg).to(device)
    teacher = mt.build_teacher(student).to(device)
    opt = torch.optim.AdamW(student.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)

    n_train = len(train_rows)
    steps_per_epoch = max(1, n_train // cfg.batch_labeled)
    n_unlabeled = len(emb_u) if use_unlabeled and emb_u is not None else 0
    rng = np.random.default_rng(cfg.seed)

    history = []
    best = {"auroc": -1.0}
    best_state = None
    global_step = 0
    start = time.time()

    for epoch in range(cfg.epochs):
        student.train()
        order = rng.permutation(n_train)
        sup_total = cons_total = 0.0
        w_cons = mt.consistency_weight(epoch, cfg) if use_unlabeled else 0.0
        decay = mt.ema_decay_for(epoch, cfg)

        for step in range(steps_per_epoch):
            batch_idx = order[step * cfg.batch_labeled : (step + 1) * cfg.batch_labeled]
            if len(batch_idx) == 0:
                continue
            x, m = _to_tensor(emb_l, mask_l, train_rows[batch_idx], device)
            y = torch.from_numpy(train_y[batch_idx]).to(device)

            logit = student(mt.augment(x, m, cfg), m)
            loss_sup = bce_loss(logit, y)
            loss = loss_sup
            sup_total += float(loss_sup.detach())

            if use_unlabeled and n_unlabeled:
                u_idx = rng.integers(0, n_unlabeled, cfg.batch_unlabeled)
                u_idx.sort()
                xu, mu = _to_tensor(emb_u, mask_u, u_idx, device)
                # Independent perturbations: the student and teacher must agree on
                # the same window seen through different noise.
                s_logit = student(mt.augment(xu, mu, cfg), mu)
                with torch.no_grad():
                    teacher.train()
                    t_logit = teacher(mt.augment(xu, mu, cfg), mu)
                loss_cons = mt.consistency_loss(s_logit, t_logit)
                loss = loss + w_cons * loss_cons
                cons_total += float(loss_cons.detach())

            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(student.parameters(), 5.0)
            opt.step()
            mt.update_teacher(student, teacher, decay, global_step)
            global_step += 1

        val_rows = val_df["row"].to_numpy()
        val_y = val_df["label"].to_numpy()
        m_student = evaluate(student, emb_l, mask_l, val_rows, val_y, device)
        m_teacher = evaluate(teacher, emb_l, mask_l, val_rows, val_y, device)
        history.append(
            {
                "epoch": epoch,
                "loss_sup": sup_total / steps_per_epoch,
                "loss_cons": cons_total / steps_per_epoch,
                "w_cons": w_cons,
                "val_student_auroc": m_student["auroc"],
                "val_teacher_auroc": m_teacher["auroc"],
            }
        )
        if m_teacher["auroc"] > best["auroc"]:
            best = m_teacher
            best_state = {
                "student": {k: v.detach().cpu().clone() for k, v in student.state_dict().items()},
                "teacher": {k: v.detach().cpu().clone() for k, v in teacher.state_dict().items()},
                "epoch": epoch,
            }
        if verbose and (epoch % log_every == 0 or epoch == cfg.epochs - 1):
            print(
                f"  epoch {epoch:3d} sup {sup_total / steps_per_epoch:.4f} "
                f"cons {cons_total / steps_per_epoch:.4f} (w={w_cons:.2f}) "
                f"val auroc student {m_student['auroc']:.4f} teacher {m_teacher['auroc']:.4f}",
                flush=True,
            )

    if best_state:
        student.load_state_dict(best_state["student"])
        teacher.load_state_dict(best_state["teacher"])

    test_rows = test_df["row"].to_numpy()
    test_y = test_df["label"].to_numpy()
    result = {
        "config": asdict(cfg),
        "use_unlabeled": use_unlabeled,
        "n_train": int(n_train),
        "n_unlabeled": int(n_unlabeled),
        "params": count_parameters(student),
        "best_epoch": best_state["epoch"] if best_state else None,
        "val_teacher": best,
        "test_student": evaluate(student, emb_l, mask_l, test_rows, test_y, device),
        "test_teacher": evaluate(teacher, emb_l, mask_l, test_rows, test_y, device),
        "history": history,
        "seconds": round(time.time() - start, 1),
    }
    result["_models"] = {"student": student, "teacher": teacher}
    return result


def save_run(result: dict, name: str) -> None:
    ensure_dirs()
    models = result.pop("_models", None)
    if models:
        torch.save(
            {
                "student": models["student"].state_dict(),
                "teacher": models["teacher"].state_dict(),
                "config": result["config"],
            },
            PATHS.models / f"{name}.pt",
        )
    (PATHS.models / f"{name}.json").write_text(json.dumps(result, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Train immunogenicity models")
    parser.add_argument("--epochs", type=int, default=TrainConfig.epochs)
    parser.add_argument("--seed", type=int, default=TrainConfig.seed)
    parser.add_argument("--arm", choices=["baseline", "mean_teacher", "both"], default="both")
    args = parser.parse_args()

    cfg = TrainConfig(epochs=args.epochs, seed=args.seed)
    device = torch_device()
    print(f"device={device}")

    labeled, _ = prepare_labeled(seed=cfg.seed)
    emb_l, mask_l = load_embeddings("labeled")
    emb_u, mask_u = load_embeddings("unlabeled")
    print(
        f"labeled windows {len(labeled)} "
        f"(train {int((labeled.split == 'train').sum())}, "
        f"val {int((labeled.split == 'val').sum())}, "
        f"test {int((labeled.split == 'test').sum())}) | unlabeled {len(emb_u)}"
    )

    arms = ["baseline", "mean_teacher"] if args.arm == "both" else [args.arm]
    summary = {}
    for arm in arms:
        print(f"\n=== {arm} ===")
        result = train_model(
            cfg,
            labeled,
            emb_l,
            mask_l,
            use_unlabeled=(arm == "mean_teacher"),
            emb_u=emb_u,
            mask_u=mask_u,
            device=device,
        )
        summary[arm] = {
            "test_teacher": result["test_teacher"],
            "test_student": result["test_student"],
        }
        save_run(result, arm)
        print(f"  test teacher: {result['test_teacher']}")

    if len(summary) == 2:
        b = summary["baseline"]["test_teacher"]["auroc"]
        m = summary["mean_teacher"]["test_teacher"]["auroc"]
        print(f"\nin-distribution AUROC  baseline {b:.4f} -> mean teacher {m:.4f} ({m - b:+.4f})")


if __name__ == "__main__":
    main()
