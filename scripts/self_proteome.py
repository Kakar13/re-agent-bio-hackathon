#!/usr/bin/env python3
"""Self-proteome similarity: the tolerance axis, computed on TCR-facing residues only.

The naive version of this feature — count exact 9-mer matches against the human
proteome — is the wrong one, and the literature says so specifically.

A CD4 T cell does not see the whole peptide. In a class II 9-mer binding core,
positions 1, 4, 6 and 9 point down into the HLA groove and positions 2, 3, 5, 7 and 8
point up at the T-cell receptor [1]. Two peptides that differ at anchors but match at
those five TCR-facing positions present the *same face* to the repertoire, so a
foreign peptide matching self on the TCR face meets T cells that were negatively
selected against it. Liu et al. built H7N9 influenza epitopes this way and measured the
consequence: across 18 naive donors, ELISpot stimulation index correlated
significantly negatively with TCR-face cross-conservation, and every peptide above
their threshold expanded FoxP3+ regulatory T cells rather than effectors [1].

So the match rule here is 100% identity at 9-mer positions 2,3,5,7,8 against UniProt
reviewed human, not whole-peptide identity. Exact 9-mer counts are computed too, purely
to show how much the naive feature would have missed.

Deviation from the published tool: JanusMatrix additionally requires the human match to
be predicted to bind at least one shared HLA supertype [1]. We have no MHC predictor in
this pipeline, so our counts are a permissive superset and carry a chance-match
background. `tcr_log2_enrichment` divides out that background under a composition null.

[1] Liu R, Moise L, Tassone R, et al. "H7N9 T-cell epitopes that mimic human sequences
    are less immunogenic and may induce Treg-mediated tolerance." Hum Vaccin Immunother
    11(9):2241-2252 (2015). PMC4635734

    uv run python scripts/self_proteome.py
"""

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from re_agent.immuno.config import AA  # noqa: E402

PROTEOME = ROOT / "data" / "raw" / "proteome" / "human_reviewed.fasta.gz"
PROCESSED = ROOT / "data" / "processed"
OUT = PROCESSED / "self_proteome.parquet"
OUT_JSON = ROOT / "results" / "reports" / "self_proteome_summary.json"

CORE = 9
# 1-indexed 2,3,5,7,8 from Liu et al.; the complement 1,4,6,9 are the HLA anchors.
TCR_FACING = np.array([1, 2, 4, 6, 7])
CODE = np.full(128, -1, dtype=np.int8)
for i, a in enumerate(AA):
    CODE[ord(a)] = i

TABLES = {
    "iedb_curated": "iedb_curated_windows.parquet",
    "denovo": "denovo_windows.parquet",
    "pda": "pda_windows.parquet",
    "reference": "reference_windows.parquet",
}


def encode(seq: str) -> np.ndarray:
    return CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]


def proteome_cores() -> np.ndarray:
    """Every human 9-mer as a (n, 9) integer matrix, dropping any with odd residues."""
    seqs = []
    cur: list[str] = []
    with gzip.open(PROTEOME, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    print(f"human proteome: {len(seqs):,} proteins")

    blocks = []
    for s in seqs:
        if len(s) < CORE:
            continue
        c = encode(s)
        idx = np.arange(len(c) - CORE + 1)[:, None] + np.arange(CORE)[None, :]
        w = c[idx]
        blocks.append(w[(w >= 0).all(1)])
    cores = np.concatenate(blocks)
    print(f"human 9-mers: {len(cores):,}")
    return cores


def tcr_key(cores: np.ndarray) -> np.ndarray:
    """Base-20 integer over the five TCR-facing positions."""
    k = np.zeros(len(cores), dtype=np.int64)
    for p in TCR_FACING:
        k = k * 20 + cores[:, p]
    return k


def full_key(cores: np.ndarray) -> np.ndarray:
    k = np.zeros(len(cores), dtype=np.int64)
    for p in range(CORE):
        k = k * 20 + cores[:, p]
    return k


def build_index(cores: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts = np.bincount(tcr_key(cores), minlength=20**5).astype(np.int32)
    exact = np.unique(full_key(cores))
    # Composition null: expected hits for a key drawn from human residue frequencies.
    freq = np.bincount(cores.ravel(), minlength=20) / cores.size
    print(f"TCR-face keyspace: {20**5:,} | occupied {int((counts > 0).sum()):,} "
          f"| mean count where occupied {counts[counts > 0].mean():.1f}")
    return counts, exact, freq


def score(seqs: list[str], counts: np.ndarray, exact: np.ndarray,
          freq: np.ndarray) -> pd.DataFrame:
    tcr_tot, tcr_max, frames_hit, exact_tot, expected = [], [], [], [], []
    logf = np.log(freq, out=np.full(20, -50.0), where=freq > 0)

    for s in seqs:
        c = encode(s)
        n = len(c) - CORE + 1
        if n < 1 or (c < 0).any():
            for acc in (tcr_tot, tcr_max, frames_hit, exact_tot):
                acc.append(0)
            expected.append(np.nan)
            continue
        idx = np.arange(n)[:, None] + np.arange(CORE)[None, :]
        w = c[idx]
        tk, fk = tcr_key(w), full_key(w)
        hits = counts[tk]
        tcr_tot.append(int(hits.sum()))
        tcr_max.append(int(hits.max()))
        frames_hit.append(int((hits > 0).sum()))
        pos = np.searchsorted(exact, fk)
        pos = np.clip(pos, 0, len(exact) - 1)
        exact_tot.append(int((exact[pos] == fk).sum()))
        # expected hits per frame under an independent-residue null
        p = np.exp(logf[w[:, TCR_FACING]].sum(1))
        expected.append(float((p * counts.sum()).sum()))

    df = pd.DataFrame({
        "seq": seqs,
        "self_tcr_matches": tcr_tot,
        "self_tcr_max_frame": tcr_max,
        "self_frames_matched": frames_hit,
        "self_exact_9mer": exact_tot,
        "expected_by_chance": expected,
    })
    df["tcr_log2_enrichment"] = np.log2(
        (df["self_tcr_matches"] + 1) / (df["expected_by_chance"] + 1)
    )
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description="Self-proteome TCR-face similarity")
    ap.add_argument("--tables", nargs="*", default=list(TABLES))
    args = ap.parse_args()

    if not PROTEOME.exists():
        raise SystemExit(f"missing {PROTEOME} -- download the human proteome first")

    counts, exact, freq = build_index(proteome_cores())

    frames = {}
    for name in args.tables:
        path = PROCESSED / TABLES[name]
        if not path.exists():
            print(f"  {name}: missing {path.name}, skipping")
            continue
        frames[name] = pd.read_parquet(path, columns=["seq"])["seq"]

    uniq = sorted(set().union(*[set(s) for s in frames.values()]))
    print(f"\nscoring {len(uniq):,} unique windows")
    scored = score(uniq, counts, exact, freq)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    scored.to_parquet(OUT, index=False)

    lut = scored.set_index("seq")
    summary = {}
    print(f"\n{'cohort':<14}{'windows':>10}{'TCR matches':>14}{'exact 9-mer':>13}"
          f"{'frames hit':>12}{'log2 enrich':>13}")
    for name, s in frames.items():
        d = lut.reindex(s.unique())
        summary[name] = {
            "unique_windows": int(len(d)),
            "median_tcr_matches": float(d["self_tcr_matches"].median()),
            "mean_tcr_matches": float(d["self_tcr_matches"].mean()),
            "frac_with_exact_9mer": float((d["self_exact_9mer"] > 0).mean()),
            "median_frames_matched": float(d["self_frames_matched"].median()),
            "median_log2_enrichment": float(d["tcr_log2_enrichment"].median()),
        }
        print(f"{name:<14}{len(d):>10,}{d['self_tcr_matches'].median():>14.0f}"
              f"{(d['self_exact_9mer'] > 0).mean():>13.1%}"
              f"{d['self_frames_matched'].median():>12.0f}"
              f"{d['tcr_log2_enrichment'].median():>13.2f}")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(f"\n-> {OUT}\n-> {OUT_JSON}")


if __name__ == "__main__":
    main()
