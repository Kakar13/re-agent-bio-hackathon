import numpy as np

from re_agent.e2e_pls import score as score_mod
from re_agent.e2e_pls.score import (
    PeptideScore,
    aggregate_protein_risk,
    compute_confidence,
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


def test_compute_confidence_bounds_and_monotone_tap():
    low_unc = compute_confidence(0.0, 0.9, 0.8, 0.8, "AAAA", "AAAA")
    high_unc = compute_confidence(2.0, 0.9, 0.8, 0.8, "AAAA", "AAAA")
    for key, value in low_unc.items():
        assert 0 <= value <= 1, key
    assert low_unc["confidence_tap"] > high_unc["confidence_tap"]
    assert low_unc["confidence"] > high_unc["confidence"]


def test_compute_confidence_missing_flanks_lowers_context():
    full = compute_confidence(0.2, 0.8, 0.7, 0.7, "AAAA", "AAAA")
    none = compute_confidence(0.2, 0.8, 0.7, 0.7, "", "")
    assert full["confidence_context"] == 1.0
    assert none["confidence_context"] == 0.0
    assert full["confidence"] > none["confidence"]


def test_compute_confidence_mhc_more_decisive_at_extremes():
    mid = compute_confidence(0.2, 0.5, 0.5, 0.5, "AAAA", "AAAA")
    extreme = compute_confidence(0.2, 0.95, 0.5, 0.5, "AAAA", "AAAA")
    assert extreme["confidence_mhc"] > mid["confidence_mhc"]


def test_score_sequence_and_aggregate(trained_heads, mock_client):
    seq = "GSHMSEEELKEAVKLLKKAEELVKKGD"
    scores = score_sequence(seq, "HLA-A*02:01", mock_client, trained_heads)
    assert len(scores) == len(tile_sequence(seq))
    for s in scores:
        assert isinstance(s, PeptideScore)
        assert 0 <= s.composite_risk <= 1
        assert 0 <= s.confidence <= 1

    risk = aggregate_protein_risk(scores, top_k=5, threshold=0.5)
    risks = np.array([s.composite_risk for s in scores])
    assert risk.max_risk == risks.max()
    assert risk.max_risk_window.composite_risk == risks.max()
    assert risk.count_above_threshold == int((risks >= 0.5).sum())
    assert risk.n_windows == len(scores)
    assert 0 <= risk.mean_confidence <= 1
    assert risk.max_risk_confidence == risk.max_risk_window.confidence


def test_aggregate_protein_risk_empty():
    risk = aggregate_protein_risk([], top_k=5, threshold=0.5)
    assert risk.n_windows == 0
    assert risk.max_risk_window is None
    assert risk.max_risk == 0.0


def test_residue_heatmap_rolls_window_risk_onto_covering_positions():
    seq = "ABCDEFGHI"  # one 9-mer window covering all positions
    scores = [
        PeptideScore(
            peptide=seq,
            start=0,
            end=9,
            cleave_n_prob=0.8,
            cleave_c_prob=0.8,
            tap_log_ic50_relative=0.0,
            tap_uncertainty=0.1,
            mhc_cosine_similarity=0.5,
            mhc_presentation_propensity=0.8,
            composite_risk=0.7,
            confidence=0.6,
            confidence_tap=0.6,
            confidence_mhc=0.6,
            confidence_context=1.0,
            confidence_agreement=1.0,
        )
    ]
    heat = score_mod.compute_residue_heatmap(seq, scores)
    assert len(heat) == 9
    assert [h.residue for h in heat] == list(seq)
    assert all(h.risk_max == 0.7 for h in heat)
    assert all(h.n_windows == 1 for h in heat)
    assert all(h.confidence_at_max == 0.6 for h in heat)


def test_residue_heatmap_max_uses_hottest_covering_window():
    seq = "ABCDEFGHIJ"  # two overlapping 9-mers
    cool = PeptideScore(
        peptide=seq[0:9],
        start=0,
        end=9,
        cleave_n_prob=0.2,
        cleave_c_prob=0.2,
        tap_log_ic50_relative=0.0,
        tap_uncertainty=0.1,
        mhc_cosine_similarity=0.1,
        mhc_presentation_propensity=0.2,
        composite_risk=0.2,
        confidence=0.4,
        confidence_tap=0.4,
        confidence_mhc=0.4,
        confidence_context=1.0,
        confidence_agreement=1.0,
    )
    hot = PeptideScore(
        peptide=seq[1:10],
        start=1,
        end=10,
        cleave_n_prob=0.9,
        cleave_c_prob=0.9,
        tap_log_ic50_relative=0.0,
        tap_uncertainty=0.1,
        mhc_cosine_similarity=0.9,
        mhc_presentation_propensity=0.9,
        composite_risk=0.9,
        confidence=0.8,
        confidence_tap=0.8,
        confidence_mhc=0.8,
        confidence_context=1.0,
        confidence_agreement=1.0,
    )
    heat = score_mod.compute_residue_heatmap(seq, [cool, hot])
    assert heat[0].risk_max == 0.2  # only the cool window covers position 0
    assert heat[0].n_windows == 1
    assert heat[1].risk_max == 0.9  # overlap: hottest wins
    assert heat[1].n_windows == 2
    assert heat[1].confidence_at_max == 0.8
    assert heat[-1].risk_max == 0.9  # only the hot window covers the last residue


def test_residue_heatmap_figure_has_two_tracks():
    from re_agent.e2e_pls.dashboard import build_residue_heatmap_figure

    heat = score_mod.compute_residue_heatmap(
        "ABCDEFGHI",
        [
            PeptideScore(
                peptide="ABCDEFGHI",
                start=0,
                end=9,
                cleave_n_prob=0.5,
                cleave_c_prob=0.5,
                tap_log_ic50_relative=0.0,
                tap_uncertainty=0.1,
                mhc_cosine_similarity=0.0,
                mhc_presentation_propensity=0.5,
                composite_risk=0.4,
                confidence=0.5,
                confidence_tap=0.5,
                confidence_mhc=0.5,
                confidence_context=1.0,
                confidence_agreement=1.0,
            )
        ],
    )
    for h in heat:
        h.mhc_importance = 0.3
    fig = build_residue_heatmap_figure(heat)
    assert len(fig.data) == 1
    assert list(fig.data[0].z[0]) == [0.4] * 9
    assert list(fig.data[0].y) == ["window risk", "MHC residue"]


def test_explain_sequence_returns_heatmap(trained_heads, mock_client):
    seq = "GSHMSEEELKEAVKLLKKAEELVKKGD"
    scores, protein, heatmap = score_mod.explain_sequence(
        seq, "HLA-A*02:01", mock_client, trained_heads
    )
    assert len(heatmap) == len(seq)
    assert len(scores) == len(tile_sequence(seq))
    assert protein.n_windows == len(scores)
    assert all(h.mhc_importance is not None for h in heatmap)
    assert all(0.0 <= h.mhc_importance <= 1.0 for h in heatmap)
    assert all(0.0 <= h.risk_max <= 1.0 for h in heatmap)


def test_score_window_matches_score_sequence(trained_heads, mock_client):
    seq = "GSHMSEEELKEAVKLLKKAEELVKKGD"
    window = tile_sequence(seq)[3]
    direct = score_mod.score_window(window, "HLA-A*02:01", mock_client, trained_heads)
    via_sequence = score_sequence(seq, "HLA-A*02:01", mock_client, trained_heads)[3]
    assert direct.peptide == via_sequence.peptide == window.peptide
    assert direct.composite_risk == via_sequence.composite_risk
