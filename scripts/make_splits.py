#!/usr/bin/env python3
"""Cluster-split an affinity CSV. Never random-split by row.

Default: MMseqs2 easy-cluster on binder sequences (identity ≥ 0.9, coverage ≥ 0.8),
then assign whole clusters to train/val/test.

--by target: hold out entire antigens (generalization to unseen targets).
--by both: a row is out-of-train if its binder cluster OR its target is held out.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]


def _require_mmseqs() -> str:
    exe = shutil.which("mmseqs")
    if not exe:
        raise SystemExit(
            "mmseqs not found on PATH. Install MMseqs2 "
            "(e.g. `brew install mmseqs2` or `conda install -c bioconda mmseqs2`) "
            "and re-run. Random row splits are not a fallback."
        )
    return exe


def _write_fasta(sequences: pd.Series, path: Path) -> int:
    n = 0
    with path.open("w") as fh:
        for i, seq in enumerate(sequences):
            if not isinstance(seq, str) or not seq:
                continue
            fh.write(f">s{i}\n{seq}\n")
            n += 1
    return n


def _cluster_binders(sequences: pd.Series, identity: float, coverage: float) -> pd.Series:
    """Return a Series aligned to `sequences` with cluster representative ids."""
    mmseqs = _require_mmseqs()
    unique = sequences.fillna("").astype(str)
    # Map each distinct non-empty sequence to a cluster.
    uniq = pd.Series(sorted({s for s in unique if s}))
    if uniq.empty:
        return pd.Series([pd.NA] * len(sequences), index=sequences.index)

    with tempfile.TemporaryDirectory(prefix="mmseqs_split_") as tmp:
        tmp_path = Path(tmp)
        fasta = tmp_path / "binders.fasta"
        out_prefix = tmp_path / "clu"
        tsv = tmp_path / "clu_cluster.tsv"
        _write_fasta(uniq, fasta)
        cmd = [
            mmseqs,
            "easy-cluster",
            str(fasta),
            str(out_prefix),
            str(tmp_path / "tmp"),
            "--min-seq-id",
            str(identity),
            "-c",
            str(coverage),
            "--cov-mode",
            "1",
        ]
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)
        # easy-cluster writes <prefix>_cluster.tsv : representative, member
        cluster_tsv = tmp_path / "clu_cluster.tsv"
        if not cluster_tsv.exists():
            # some versions name it clu_cluster.tsv vs prefix_rep_seq
            candidates = list(tmp_path.glob("*cluster.tsv"))
            if not candidates:
                raise SystemExit(f"MMseqs2 produced no cluster TSV under {tmp_path}")
            cluster_tsv = candidates[0]
        pairs = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["rep", "member"])

    # FASTA ids were s0..sN matching uniq order
    id_to_seq = {f"s{i}": seq for i, seq in enumerate(uniq)}
    member_to_rep: dict[str, str] = {}
    for _, row in pairs.iterrows():
        member_seq = id_to_seq.get(str(row["member"]))
        rep_seq = id_to_seq.get(str(row["rep"]), str(row["rep"]))
        if member_seq:
            member_to_rep[member_seq] = rep_seq

    # sequences that never clustered (empty) stay NA
    return unique.map(lambda s: member_to_rep.get(s, s if s else pd.NA))


def _assign_groups(groups: pd.Series, train: float, val: float, seed: int) -> dict:
    uniq = pd.Series(groups.dropna().unique()).sample(frac=1.0, random_state=seed)
    n = len(uniq)
    n_train = int(n * train)
    n_val = int(n * val)
    split_of = {}
    for i, g in enumerate(uniq):
        if i < n_train:
            split_of[g] = "train"
        elif i < n_train + n_val:
            split_of[g] = "val"
        else:
            split_of[g] = "test"
    return split_of


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True, help="Affinity CSV from build_affinity_dataset.py")
    p.add_argument("--out", required=True, help="CSV with an added `split` column")
    p.add_argument(
        "--by",
        choices=("binder", "target", "both"),
        default="binder",
        help="binder = MMseqs2 cluster split; target = hold out antigens; both = either",
    )
    p.add_argument("--identity", type=float, default=0.9)
    p.add_argument("--coverage", type=float, default=0.8)
    p.add_argument("--train", type=float, default=0.8)
    p.add_argument("--val", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=13)
    args = p.parse_args()
    if args.train + args.val >= 1.0:
        raise SystemExit("--train + --val must be < 1 (remainder is test)")

    df = pd.read_csv(args.inp)
    if "binder_sequence" not in df.columns or "target_sequence" not in df.columns:
        raise SystemExit(f"Expected binder_sequence and target_sequence. Columns: {list(df.columns)}")

    binder_split = None
    target_split = None

    if args.by in {"binder", "both"}:
        print(f"Clustering {df['binder_sequence'].nunique():,} unique binders "
              f"(identity ≥ {args.identity}, coverage ≥ {args.coverage})")
        clusters = _cluster_binders(df["binder_sequence"], args.identity, args.coverage)
        df["binder_cluster_id"] = clusters
        binder_split = _assign_groups(clusters, args.train, args.val, args.seed)

    if args.by in {"target", "both"}:
        targets = df["target_sequence"].fillna("").astype(str)
        # empty target sequences (neg controls) stay together
        target_split = _assign_groups(targets, args.train, args.val, args.seed + 1)
        df["target_id"] = targets

    def _row_split(row: pd.Series) -> str:
        labels = []
        if binder_split is not None:
            labels.append(binder_split.get(row.get("binder_cluster_id"), "train"))
        if target_split is not None:
            labels.append(target_split.get(row.get("target_id"), "train"))
        if args.by == "both":
            # unseen if either axis is held out
            if "test" in labels:
                return "test"
            if "val" in labels:
                return "val"
            return "train"
        return labels[0]

    df["split"] = df.apply(_row_split, axis=1)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"\n=== splits → {out} ===")
    for name, n in df["split"].value_counts().sort_index().items():
        print(f"  {name:5} {n:,}  ({n / len(df):.1%})")
    if "binder_cluster_id" in df:
        print(f"  binder clusters {df['binder_cluster_id'].nunique():,}")
    if args.by in {"target", "both"}:
        for split in ("train", "val", "test"):
            n_t = df.loc[df["split"] == split, "target_sequence"].nunique(dropna=True)
            print(f"  unique targets in {split}: {n_t:,}")


if __name__ == "__main__":
    main()
