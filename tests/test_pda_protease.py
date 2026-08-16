"""Tests for the cleavage-first protease attribution pipeline.

These cover the logic that would silently produce wrong biology if it broke:
Schechter-Berger window alignment, the digest constraints that are the whole
point of the cleavage-first ordering, and the controls that make the structural
claim falsifiable. Nothing here touches the network.
"""

from __future__ import annotations

import random

import pytest

from re_agent.pda_protease import cleavage, controls, digest, merops
from re_agent.pda_protease.accessibility import (
    NATIVE_ACCESSIBLE,
    UNFOLDING_REQUIRED,
    ResidueStructure,
    span_features,
)


def _matrix(name: str = "test_protease", *, favour: str = "L", at: str = "P1") -> merops.SpecificityMatrix:
    """A matrix that strongly prefers one residue at one position."""
    counts = {p: dict.fromkeys(merops.STANDARD_AA, 10) for p in merops.POSITIONS}
    counts[at][favour] = 5000
    m = merops.SpecificityMatrix(
        name=name, merops_id="X00.000", family="test", counts=counts, n_cleavages=1000
    )
    m.build_pwm({aa: 1 / 20 for aa in merops.STANDARD_AA})
    return m


class TestWindows:
    def test_window_is_p4_to_p4prime_around_the_bond(self):
        seq = "AAAACCCCGGGG"
        # A bond before index 4 means P4-P1 = AAAA and P1'-P4' = CCCC.
        w = merops.window_at(seq, 4)
        assert w == "AAAACCCC"
        assert w[merops.N_NONPRIME - 1] == "A"  # P1
        assert w[merops.N_NONPRIME] == "C"  # P1'

    def test_window_none_when_it_would_run_off_either_end(self):
        seq = "ACDEFGHIKL"
        assert merops.window_at(seq, 0) is None
        assert merops.window_at(seq, 3) is None
        assert merops.window_at(seq, len(seq)) is None
        assert merops.window_at(seq, 4) is not None

    def test_scissile_bond_label_matches_the_window(self):
        site = cleavage.CleavageSite(
            design_id="d", cut_index=10, window="AAAALCCC", protease="p",
            pwm_score=1.0, p_value=0.01, z_score=1.0, netchop_score=None,
            mean_rsa=None, mean_bfactor=None, coil_fraction=0.0,
            observed_fraction=1.0, accessibility="x",
        )
        # P1 is the fourth residue of the window, P1' the fifth.
        assert site.scissile_bond == "A10|L11"


class TestScoring:
    def test_preferred_residue_scores_above_a_disfavoured_one(self):
        m = _matrix(favour="L", at="P1")
        good = m.score("AAALAAAA")  # L at P1
        bad = m.score("AAAGAAAA")  # G at P1
        assert good > bad

    def test_vectorised_scan_matches_scoring_one_window_at_a_time(self):
        m = _matrix(favour="W", at="P2")
        seq = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
        pwm = cleavage.pwm_array(m)
        vec = cleavage.score_all_windows(cleavage.encode(seq), pwm)
        for k in range(len(vec)):
            assert vec[k] == pytest.approx(m.score(seq[k : k + 8]), abs=1e-9)

    def test_p_value_never_zero_and_strong_site_is_significant(self):
        m = _matrix(favour="W", at="P1")
        # One tryptophan, positioned so its bond is scoreable.
        seq = "AAAAAAAAAAAAWAAAAAAAAAAAA"
        sites = cleavage.scan_design("d", seq, {"test": m}, n_shuffles=100)
        assert all(s.p_value > 0 for s in sites)
        best = min(sites, key=lambda s: s.p_value)
        assert best.window[merops.N_NONPRIME - 1] == "W"

    def test_background_reflects_the_corpus_given(self):
        bg = merops.background_from_sequences(["AAAAAAAAAC"])
        assert bg["A"] == pytest.approx(0.9)
        assert bg["W"] > 0  # floored, never zero


class TestAccessibility:
    def _features(self, n: int, rsa: float, sse: str) -> list[ResidueStructure]:
        return [
            ResidueStructure(seq_index=i, residue="A", observed=True, rsa=rsa, bfactor=30.0, sse=sse)
            for i in range(n)
        ]

    def test_exposed_span_is_native_accessible(self):
        f = self._features(20, rsa=0.6, sse="c")
        assert span_features(f, 10).classification == NATIVE_ACCESSIBLE

    def test_buried_helix_requires_unfolding(self):
        f = self._features(20, rsa=0.02, sse="a")
        assert span_features(f, 10).classification == UNFOLDING_REQUIRED

    def test_buried_but_coil_is_still_reachable(self):
        f = self._features(20, rsa=0.02, sse="c")
        assert span_features(f, 10).classification == NATIVE_ACCESSIBLE

    def test_disordered_region_counts_as_reachable(self):
        f = self._features(20, rsa=0.02, sse="a")
        for i in range(6, 14):
            f[i].observed = False
        assert span_features(f, 10).classification == NATIVE_ACCESSIBLE


class TestDigest:
    def test_class_ii_peptides_never_cross_a_confident_cut(self):
        seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(120))
        cuts = {30, 70}
        peptides = digest.digest_mhcii("d", seq, cuts, [])
        for p in peptides:
            for c in cuts:
                assert not (p.start < c < p.end), f"{p.peptide} spans cut {c}"

    def test_fragments_tile_the_whole_sequence(self):
        seq = "A" * 100
        frags = digest.build_fragments("d", seq, {25, 60}, [])
        assert [(f.start, f.end) for f in frags] == [(0, 25), (25, 60), (60, 100)]
        assert "".join(f.sequence for f in frags) == seq

    def test_class_i_c_terminus_is_always_a_proteasome_site(self):
        seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(60))
        cuts = {20, 40}
        peptides = digest.digest_mhci("d", seq, cuts)
        assert peptides
        for p in peptides:
            assert p.end - 1 in cuts
            assert p.c_term_source == "proteasome"
            assert 8 <= p.length <= 11

    def test_constrained_digest_is_smaller_than_the_naive_scan(self):
        seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(120))
        constrained = digest.digest_mhcii("d", seq, {30, 60, 90}, [])
        naive = digest.unconstrained_peptides("d", seq)
        assert len(constrained) < len(naive)

    def test_terminus_provenance_only_claims_a_protease_at_fragment_edges(self):
        seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(80))
        peptides = digest.digest_mhcii("d", seq, {40}, [])
        interior = [p for p in peptides if p.start not in (0, 40)]
        assert all(p.n_term_source == "internal" for p in interior)


class TestControls:
    def test_scramble_preserves_composition_and_length(self):
        rng = random.Random(0)
        seg = "ACDEFGHIKLMNPQR"
        out = controls.scramble_segment(seg, rng)
        assert sorted(out) == sorted(seg)
        assert len(out) == len(seg)

    def test_swap_picks_a_different_and_worse_scoring_protease(self):
        sites = [
            cleavage.CleavageSite(
                design_id="d", cut_index=10, window="AAAALAAA", protease=p,
                pwm_score=score, p_value=0.1, z_score=0.0, netchop_score=None,
                mean_rsa=None, mean_bfactor=None, coil_fraction=0.0,
                observed_fraction=1.0, accessibility="x",
            )
            for p, score in [("cathepsin_S", 9.0), ("cathepsin_L", 1.0), ("cathepsin_K", -4.0)]
        ]
        swap = controls.pick_swap_protease(sites, 10, "cathepsin_S")
        assert swap == "cathepsin_K"

    def test_decoys_avoid_the_real_sites(self):
        sites = [
            cleavage.CleavageSite(
                design_id="d", cut_index=i, window="AAAAAAAA", protease="cathepsin_S",
                pwm_score=-float(i), p_value=0.9, z_score=0.0, netchop_score=None,
                mean_rsa=None, mean_bfactor=None, coil_fraction=0.0,
                observed_fraction=1.0, accessibility="native-accessible",
            )
            for i in range(4, 30)
        ]
        real = {5, 6, 7}
        decoys = controls.pick_decoy_sites(sites, "cathepsin_S", real, n=3)
        assert decoys
        assert all(d.cut_index not in real for d in decoys)


class TestJobConstruction:
    def test_segment_centres_on_the_bond_and_keeps_the_offset(self):
        seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(100))
        from re_agent.pda_protease.structure import make_segment

        seg, start, cut_off = make_segment(seq, 50, flank=14)
        assert seg == seq[36:64]
        assert start == 36
        # The residue at the cut offset within the segment is the one that
        # follows the scissile bond in the parent sequence.
        assert seg[cut_off] == seq[50]

    def test_segment_clamps_at_the_terminus_without_losing_alignment(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        from re_agent.pda_protease.structure import make_segment

        seg, start, cut_off = make_segment(seq, 3, flank=14)
        assert start == 0
        assert seg[cut_off] == seq[3]

    def test_p1_is_the_residue_before_the_bond(self):
        from re_agent.pda_protease.structure import build_job

        seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(100))
        job = build_job(design_id="d", sequence=seq, cut_index=50, protease="cathepsin_S")
        assert job.segment[job.p1_offset] == seq[49]
        assert job.segment[job.cut_offset] == seq[50]
