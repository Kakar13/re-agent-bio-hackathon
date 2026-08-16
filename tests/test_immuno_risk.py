"""Unit tests for immuno-risk (offline-friendly with IMMUNO_ALLOW_HEURISTIC_MHC=1)."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

os.environ.setdefault("IMMUNO_ALLOW_HEURISTIC_MHC", "1")

from re_agent.immuno_risk.aggregation import aggregation_report
from re_agent.immuno_risk.peptides import clean_sequence, sliding_windows
from re_agent.immuno_risk.pipeline import run_immuno_risk
from re_agent.immuno_risk.reference_data import ensure_fixtures, load_overlap_8mers
from re_agent.immuno_risk.risk import project_residue_risk
from re_agent.immuno_risk.schemas import PeptideEvidence, PeptideHit


INSULIN = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)


def test_fixtures_and_manifests():
    paths = ensure_fixtures()
    assert paths["tcell"].exists()
    assert paths["atlas"].exists()
    assert paths["overlap_ledger"].exists()
    root = Path(__file__).resolve().parents[1]
    manifests = list((root / "data" / "processed" / "immuno").glob("*_manifest.json"))
    assert manifests


def test_overlap_ledger_8mer():
    ledger = load_overlap_8mers()
    assert "GILGFVFT" in ledger


def test_aggregation_deterministic():
    a = aggregation_report("x", INSULIN)
    b = aggregation_report("x", INSULIN)
    assert a.score0to100 == b.score0to100
    assert a.free_cysteine_count == INSULIN.count("C")


def test_residue_projection():
    seq = "GILGFVFTLAAA"
    ev = [
        PeptideEvidence(
            peptide="GILGFVFTL",
            allele="HLA-A*02:01",
            mhc_class="I",
            start=0,
            end=9,
            mhc=PeptideHit(
                peptide="GILGFVFTL",
                allele="HLA-A*02:01",
                mhc_class="I",
                length=9,
                binder=True,
                method="test",
            ),
            presentation_points=3,
            tolerance_points=1,
            point_score=4,
            contribution=0.8,
        )
    ]
    rr = project_residue_risk(seq, ev)
    assert rr[0].risk > 0
    assert rr[9].risk == 0


def test_end_to_end_pipeline(tmp_path: Path):
    result = run_immuno_risk(
        INSULIN,
        sequence_id="insulin_test",
        write=True,
        results_root=tmp_path,
    )
    assert result.risk.overall in {"low", "moderate", "high"}
    assert result.risk.method == "dual_arm_existing_predictor_points_v1"
    assert 0 <= result.risk.total_points <= result.risk.max_points
    assert result.risk.mhc_i is not None
    assert result.risk.mhc_ii is not None
    assert "iedb_risk_head" not in result.predictor_versions
    assert result.aggregation is not None
    assert result.confidence.score0to1 >= 0
    assert result.artifact_dir
    art = Path(result.artifact_dir)
    for name in [
        "manifest.json",
        "summary.json",
        "peptides.csv",
        "residue-risk.json",
        "aggregation.json",
        "report.md",
    ]:
        assert (art / name).exists()
    manifest = json.loads((art / "manifest.json").read_text())
    assert manifest["run_id"] == result.run_id
    # No stub_hash in methods
    methods = {e.mhc.method for e in result.peptides}
    assert "stub_hash_rank_v0" not in methods
    assert {e.mhc_class for e in result.peptides} == {"I", "II"}


def test_benchling_dry_run():
    from re_agent.immuno_risk.benchling import publish_run, pull_candidates

    rows = pull_candidates(dry_run=True)
    assert rows and rows[0].get("dry_run")
    # publish dry-run needs a fake run dir with manifest
    d = Path("/tmp")  # will fail without files — create temp via pipeline first


def test_benchling_publish_dry_run(tmp_path: Path):
    from re_agent.immuno_risk.benchling import publish_run

    result = run_immuno_risk("GILGFVFTLAAAAKKKLLL", sequence_id="tiny", write=True, results_root=tmp_path)
    out = publish_run(Path(result.artifact_dir), dry_run=True)
    assert out["status"] == "dry_run"
    assert out["run_id"] == result.run_id


def test_sliding_windows():
    wins = sliding_windows("ACDEFGHIKLMNPQRSTVWY", range(8, 10))
    assert all(len(w[2]) in (8, 9) for w in wins)
    assert clean_sequence("acdeX123") == "ACDE"


def test_pda_designed_chains_to_netmhcpan_fasta(tmp_path: Path):
    from re_agent.immuno_risk.pda import extract_designed_chains, ingest_pda, netmhcpan_allele

    entries = [
        {
            "pdb": "1qys",
            "release_date": "2003-01-01",
            "subtitle": "Top7 de novo",
            "tags": ["de novo", "binder"],
            "seq_max_sim_natural": {"partner": "1lfw", "sim": 40.0},
            "chains": [
                {
                    "chain_id": "A,B",
                    "chain_type": "D",
                    "chain_source": "synthetic construct",
                    "chain_seq_nat": "GILGFVFTLAAAAKKKLLLGGGGLEHHHHHH",
                },
                {
                    "chain_id": "Z",
                    "chain_type": "N",
                    "chain_source": "Homo sapiens",
                    "chain_seq_nat": "MALWMRLLPLLALLALWGPDPAAA",
                },
                {
                    "chain_id": "B",
                    "chain_type": "D",
                    "chain_source": "synthetic construct",
                    "chain_seq_nat": "GGTTTACCGATATACACCCCTAAGAG",
                },
            ],
        }
    ]
    rows = extract_designed_chains(entries, strip_his=True)
    assert len(rows) == 1
    assert rows[0]["sequence_id"] == "1qys_A"
    assert rows[0]["sequence"].endswith("GGG")
    assert "H" not in rows[0]["sequence"][-6:]
    assert netmhcpan_allele("HLA-A*02:01") == "HLA-A02:01"

    fixture = tmp_path / "pda_mini.json"
    fixture.write_text(json.dumps(entries))
    summary = ingest_pda(
        source="json",
        json_path=fixture,
        out_dir=tmp_path / "out",
        write_peptides=True,
    )
    fasta = Path(summary["fasta"])
    text = fasta.read_text()
    assert text.startswith(">1qys_A")
    assert "GILGFVFTL" in text
    peptides = (tmp_path / "out" / "netmhcpan_peptides.txt").read_text().splitlines()
    assert "GILGFVFTL" in peptides
    assert "HLA-A02:01" in Path(summary["alleles"]).read_text()
    assert "netMHCpan -f" in summary["netmhcpan_cmd"]


def test_netmhcpan_parser_unit():
    from re_agent.immuno_risk.netmhcpan import _parse_netmhcpan_table

    sample = """
# Rank Threshold for Strong binder   0.500
 Pos            MHC           Peptide      Core Of Gp Gl Ip Il             Icore        Identity Score_EL %Rank_EL BindLevel
  34    HLA-A*02:01         HLVEALYLV HLVEALYLV  0  0  0  0  0         HLVEALYLV sp_P01308_INS_H 0.939284    0.030 <= SB
"""
    hits = _parse_netmhcpan_table(sample)
    assert hits
    assert hits[0].peptide == "HLVEALYLV"
    assert hits[0].percentile_rank == pytest.approx(0.03)
    assert hits[0].provenance.get("identity") == "sp_P01308_INS_H"


def test_cleavage_subsite_and_cathepsin():
    from re_agent.immuno_risk.cleavage import predict_cleavage

    # CatS needs hydrophobic P2 (V/L/I/M/F). Sequence: ...L K... → P2=L, P1=K at index 6.
    seq = "AAAVALKAAA"  # indices: 5=L, 6=K
    events = predict_cleavage(seq, site_ids=["cathepsin_s", "trypsin_kr"])
    cats = [e for e in events if e.site_id == "cathepsin_s"]
    assert any(e.position == 6 and e.p1 == "K" and e.p2 == "L" for e in cats)
    # Trypsin should also cut at K
    assert any(e.site_id == "trypsin_kr" and e.position == 6 for e in events)

    # Furin R-X-K-R (use only standard AA — clean_sequence strips X)
    furin_seq = "AAARAKRYYY"
    fe = predict_cleavage(furin_seq, site_ids=["furin_rxkr"])
    assert any(e.p1 == "R" and e.n_terminal_product.endswith("RAKR") for e in fe)


def test_batch_peptide_index(tmp_path: Path):
    from re_agent.immuno_risk.batch import build_peptide_index, parse_fasta

    fasta = tmp_path / "mini.fasta"
    fasta.write_text(">a\nGILGFVFTLAAA\n>b\nGILGFVFTLBBB\n")
    records = parse_fasta(fasta)
    peptides, index = build_peptide_index(records)
    assert "GILGFVFTL" in peptides
    # Shared 9mer appears in both chains
    locs = index["GILGFVFTL"]
    assert {c for c, _, _ in locs} == {"a", "b"}


def test_join_creating_destroying(tmp_path: Path):
    from re_agent.immuno_risk.cleavage import join_epitope_protease

    mhc = tmp_path / "mhc.csv"
    mhc.write_text(
        "chain_id,start,end,peptide,allele,percentile_rank,presentation_score,"
        "affinity_nm,processing_score,el_score,binder,method\n"
        "x,2,11,ABCDEFGHI,HLA-A*02:01,0.5,0.9,,,1,test\n"
    )
    cleave = tmp_path / "cleave.csv"
    # creating_n: cut at start-1=1; creating_c: cut at end-1=10; destroying: cut at 5
    cleave.write_text(
        "chain_id,site_id,site_name,protease_class,position,p1,p1_prime,p2,p3,score\n"
        "x,cathepsin_s,CatS,cysteine,1,A,B,X,Y,1.0\n"
        "x,immunoproteasome_b5i,b5i,threonine,10,I,J,H,G,1.0\n"
        "x,trypsin_kr,Trypsin,serine,5,F,G,E,D,1.0\n"
    )
    out = tmp_path / "join.csv"
    short = tmp_path / "short.csv"
    summary = join_epitope_protease(mhc, cleave, out_csv=out, shortlist_csv=short, shortlist_n=10)
    assert summary["n_creating"] == 2
    assert summary["n_destroying"] == 1
    text = out.read_text()
    assert "creating_n" in text and "creating_c" in text and "destroying" in text
