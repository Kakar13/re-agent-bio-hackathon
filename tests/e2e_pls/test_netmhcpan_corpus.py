from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.e2e_pls import fixtures
from re_agent.e2e_pls.netmhcpan_corpus import (
    _sample_by_split,
    attach_netmhcpan_labels,
    build_full_pda_corpus,
    build_pda_challenge_pool,
)
from re_agent.immuno.netmhcpan import NetMHCpanTeacher


class _FakeTeacher:
    def label(self, frame: pd.DataFrame) -> pd.DataFrame:
        out = frame.copy()
        out["netmhcpan_el_score"] = 0.4
        out["netmhcpan_el_rank"] = 3.0
        out["netmhcpan_ba_score"] = 0.5
        out["netmhcpan_ba_rank"] = 2.0
        out["netmhcpan_ba_ic50_nm"] = 125.0
        return out


class _FakeParentTeacher:
    def label_parent_sequences(
        self,
        parents: pd.DataFrame,
        *,
        id_column: str,
        sequence_column: str,
        parents_per_batch: int,
    ) -> pd.DataFrame:
        assert parents_per_batch > 0
        rows = []
        for parent in parents.to_dict(orient="records"):
            sequence = parent[sequence_column]
            metadata = {key: value for key, value in parent.items() if key != sequence_column}
            for start in range(len(sequence) - 8):
                peptide = sequence[start : start + 9]
                rank = 0.5 if peptide == "SLYNTVATL" else 3.0
                rows.append(
                    {
                        **metadata,
                        "start": start,
                        "end": start + 9,
                        "peptide": peptide,
                        "n_flank": sequence[max(0, start - 4) : start],
                        "c_flank": sequence[start + 9 : start + 13],
                        "length": 9,
                        "hla_allele": "HLA-A*02:01",
                        "netmhcpan_el_score": 0.9,
                        "netmhcpan_el_rank": rank,
                        "netmhcpan_ba_score": 0.8,
                        "netmhcpan_ba_rank": rank,
                        "netmhcpan_ba_ic50_nm": 25.0,
                        "netmhcpan_allele": "HLA-A*02:01",
                        "netmhcpan_version": "4.1",
                        "netmhcpan_source": "fake",
                    }
                )
        return pd.DataFrame(rows)


def _fake_iedb_response(data: dict[str, str]) -> str:
    is_ba = "netmhcpan_ba" in data["method"]
    values = "ic50\tpercentile_rank" if is_ba else "score\tpercentile_rank"
    header = f"allele\tseq_num\tstart\tend\tlength\tpeptide\tcore\ticore\t{values}"
    sequences = [
        line for line in data["sequence_text"].splitlines() if line and not line.startswith(">")
    ]
    rows = []
    for seq_num, peptide in enumerate(sequences, start=1):
        value = "50.0\t0.1" if is_ba else "0.8\t0.1"
        rows.append(
            f"HLA-A*02:01\t{seq_num}\t1\t9\t9\t{peptide}\t{peptide}\t{peptide}\t{value}"
        )
    return "\n".join([header, *rows, ""])


def test_neutral_sampling_preserves_all_optimization_splits() -> None:
    frame = pd.DataFrame(
        [
            {"row": row, "split": split}
            for split in ("train", "val", "test")
            for row in range(100)
        ]
    )

    sampled = _sample_by_split(frame, 100, seed=5)

    assert len(sampled) == 100
    assert sampled["split"].value_counts().to_dict() == {"train": 80, "val": 10, "test": 10}
    assert not sampled[["row", "split"]].duplicated().any()


def test_netmhcpan_labels_populate_binding_columns_without_collapsing_el() -> None:
    frame = fixtures.load_dev_fixture()

    labeled = attach_netmhcpan_labels(frame, _FakeTeacher())

    assert np.allclose(labeled["mhc_affinity_nm"], 125.0)
    assert np.allclose(labeled["mhc_percentile"], 2.0)
    assert np.allclose(labeled["netmhcpan_el_rank"], 3.0)
    assert labeled["label_model_version"].eq("netmhcpan=4.1;channels=el,ba").all()


def test_ba_api_rows_without_score_are_normalized(tmp_path: Path, monkeypatch) -> None:
    teacher = NetMHCpanTeacher(tmp_path)
    monkeypatch.setattr(teacher, "_post_with_retries", _fake_iedb_response)
    frame = pd.DataFrame({"peptide": ["SLYNTVATL", "GILGFVFTL"]})

    labeled = teacher.label(frame)

    assert labeled["netmhcpan_ba_ic50_nm"].eq(50.0).all()
    assert labeled["netmhcpan_ba_score"].between(0, 1).all()
    assert labeled["netmhcpan_el_score"].eq(0.8).all()


def test_pda_challenge_excludes_training_peptides(tmp_path: Path) -> None:
    path = tmp_path / "pda.parquet"
    pd.DataFrame(
        [
            {
                "parent": "pda:test:0",
                "seq": "ACDEFGHIKLMN",
                "release_date": "2025-01-01",
                "tm_natural": 0.4,
                "novelty_bin": "novel",
            }
        ]
    ).to_parquet(path, index=False)

    challenge = build_pda_challenge_pool(
        pda_designs_path=path,
        target_rows=10,
        seed=0,
        exclude_peptides={"CDEFGHIKL"},
    )

    assert "CDEFGHIKL" not in set(challenge["peptide"])


def test_full_pda_corpus_deduplicates_and_keeps_shared_parents_in_one_fold(
    tmp_path: Path,
) -> None:
    path = tmp_path / "pda.parquet"
    peptides = [
        "SLYNTVATL",
        "SLYNTVATL",
        "GILGFVFTL",
        "LLFGYPVYV",
        "NLVPMVATV",
        "GLCTLVAML",
    ]
    pd.DataFrame(
        [
            {
                "parent": f"pda:test:{index}",
                "seq": peptide,
                "release_date": "2025-01-01",
                "novelty_bin": "novel",
            }
            for index, peptide in enumerate(peptides)
        ]
    ).to_parquet(path, index=False)

    corpus, metadata = build_full_pda_corpus(
        pda_designs_path=path,
        teacher=_FakeParentTeacher(),
        parent_batch_size=2,
        seed=0,
    )

    shared = corpus.loc[corpus["peptide"] == "SLYNTVATL"].iloc[0]
    assert len(corpus) == 5
    assert shared["occurrence_count"] == 2
    assert shared["parent_count"] == 2
    assert shared["binder_class"] == "strong"
    assert set(corpus["cv_fold"]) == set(range(5))
    assert metadata["n_windows_before_deduplication"] == 6


def test_parent_table_validation_retries_transient_semantic_failure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    teacher = NetMHCpanTeacher(tmp_path)
    attempts = 0

    def flaky(*args, **kwargs):
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise RuntimeError("netmhcpan_ba-4.1 returned invalid seq_num=163")
        return pd.DataFrame({"peptide": ["SLYNTVATL"]})

    monkeypatch.setattr(teacher, "_label_parent_batch", flaky)
    monkeypatch.setattr("re_agent.immuno.netmhcpan.time.sleep", lambda _seconds: None)

    result = teacher._label_parent_batch_with_validation_retries(
        [{"parent_sequence_id": "p1", "sequence": "SLYNTVATL"}],
        "ba",
        id_column="parent_sequence_id",
        sequence_column="sequence",
    )

    assert attempts == 2
    assert result["peptide"].tolist() == ["SLYNTVATL"]
