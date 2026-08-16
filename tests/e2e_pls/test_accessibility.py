from __future__ import annotations

import numpy as np
import pytest

from re_agent.e2e_pls.accessibility import (
    BURIED_RSA_THRESHOLD,
    MAX_ACCESSIBLE_SURFACE_AREA,
    align_chain_to_design,
    project_onto_design,
    relative_accessibility,
    window_accessibility,
)


def test_relative_accessibility_normalizes_by_residue_maximum() -> None:
    codes = ["G", "W"]
    absolute = np.array([52.0, 285.0])

    rsa = relative_accessibility(codes, absolute)

    assert rsa[0] == pytest.approx(52.0 / MAX_ACCESSIBLE_SURFACE_AREA["G"])
    assert rsa[1] == pytest.approx(1.0)


def test_relative_accessibility_marks_unknown_residues_as_nan() -> None:
    rsa = relative_accessibility(["A", "X"], np.array([64.5, 100.0]))

    assert rsa[0] == pytest.approx(0.5)
    assert np.isnan(rsa[1])


@pytest.mark.parametrize(
    ("design", "chain", "expected"),
    [
        ("ABCDEFGH", "ABCDEFGH", 0),
        # Disordered termini are the common case: the chain is a substring.
        ("ABCDEFGH", "CDEF", 2),
        # An expression tag the design record omits gives a negative offset.
        ("CDEF", "ABCDEFGH", -2),
        ("ABCDEFGH", "WXYZ", None),
        ("ABCDEFGH", "", None),
    ],
)
def test_align_chain_to_design(design: str, chain: str, expected: int | None) -> None:
    assert align_chain_to_design(design, chain) == expected


def test_project_onto_design_leaves_unobserved_positions_missing() -> None:
    projected = project_onto_design(6, 2, np.array([0.4, 0.6]))

    assert np.isnan(projected[[0, 1, 4, 5]]).all()
    assert projected[2] == pytest.approx(0.4)
    assert projected[3] == pytest.approx(0.6)


def test_window_accessibility_skips_windows_without_enough_coverage() -> None:
    rsa = np.array([0.1] * 9 + [np.nan] * 7 + [0.9, 0.9])
    starts = np.array([0, 9])
    ends = np.array([9, 18])

    mean_rsa, buried_fraction, observed = window_accessibility(rsa, starts, ends)

    assert mean_rsa[0] == pytest.approx(0.1)
    assert buried_fraction[0] == pytest.approx(1.0)
    assert observed[0] == 9
    # Only two of nine residues are resolved, so the window reports nothing.
    assert np.isnan(mean_rsa[1])
    assert observed[1] == 2


def test_window_accessibility_reports_partial_burial() -> None:
    exposed = BURIED_RSA_THRESHOLD + 0.3
    rsa = np.array([0.05] * 3 + [exposed] * 6)

    mean_rsa, buried_fraction, observed = window_accessibility(
        rsa, np.array([0]), np.array([9])
    )

    assert observed[0] == 9
    assert buried_fraction[0] == pytest.approx(3 / 9)
    assert mean_rsa[0] > BURIED_RSA_THRESHOLD
