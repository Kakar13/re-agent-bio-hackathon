import numpy as np

from re_agent.e2e_pls import score as score_mod
from re_agent.e2e_pls.score import (
    PeptideScore,
    aggregate_protein_risk,
    compute_peptide_risk,
    score_sequence,
    tile_sequence,
)


def test_tile_sequence_window_count_and_bounds():
    seq = "ABCDEFGHIJKLM"  # length 13 -> 13-9+1 = 5 windows
    windows = tile_sequence(seq, peptide_len=9, flank_len=4)
    assert len(windows) == 5
    assert windows[0].peptide == seq[0:9]
    assert windows[0].n_flank == ""  # nothing before position 0
    assert windows[-1].c_flank == ""  # nothing after the last window
    for w in windows:
        assert len(w.peptide) == 9
        assert seq[w.start : w.end] == w.peptide


def test_tile_sequence_flank_length_in_interior():
    seq = "A" * 30
    windows = tile_sequence(seq, peptide_len=9, flank_len=4)
    mid = windows[10]
    assert len(mid.n_flank) == 4
    assert len(mid.c_flank) == 4


def test_compute_peptide_risk_monotonic():
    low = compute_peptide_risk(0.1, 0.1, 0.1)
    high = compute_peptide_risk(0.9, 0.9, 0.9)
    assert 0 <= low <= 1
    assert 0 <= high <= 1
    assert high > low


def test_compute_peptide_risk_zero_component_pulls_down():
    with_zero = compute_peptide_risk(0.0, 0.9, 0.9)
    without_zero = compute_peptide_risk(0.9, 0.9, 0.9)
    assert with_zero < without_zero


def test_score_sequence_and_aggregate(trained_heads, mock_client):
    seq = "GSHMSEEELKEAVKLLKKAEELVKKGD"
    scores = score_sequence(seq, "HLA-A*02:01", mock_client, trained_heads)
    assert len(scores) == len(tile_sequence(seq))
    for s in scores:
        assert isinstance(s, PeptideScore)
        assert 0 <= s.composite_risk <= 1

    risk = aggregate_protein_risk(scores, top_k=5, threshold=0.5)
    risks = np.array([s.composite_risk for s in scores])
    assert risk.max_risk == risks.max()
    assert risk.max_risk_window.composite_risk == risks.max()
    assert risk.count_above_threshold == int((risks >= 0.5).sum())
    assert risk.n_windows == len(scores)


def test_aggregate_protein_risk_empty():
    risk = aggregate_protein_risk([], top_k=5, threshold=0.5)
    assert risk.n_windows == 0
    assert risk.max_risk_window is None
    assert risk.max_risk == 0.0


def test_score_window_matches_score_sequence(trained_heads, mock_client):
    seq = "GSHMSEEELKEAVKLLKKAEELVKKGD"
    window = tile_sequence(seq)[3]
    direct = score_mod.score_window(window, "HLA-A*02:01", mock_client, trained_heads)
    via_sequence = score_sequence(seq, "HLA-A*02:01", mock_client, trained_heads)[3]
    assert direct.peptide == via_sequence.peptide == window.peptide
    assert direct.composite_risk == via_sequence.composite_risk
