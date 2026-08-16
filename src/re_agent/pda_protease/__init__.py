"""Cleavage-first protease attribution for de novo designs.

Pipeline order follows the antigen-processing pathway rather than predictor
convenience: proteolysis is resolved first, the digest produces peptides, and
only those peptides are scored for MHC presentation.

    pda -> merops -> accessibility -> cleavage -> digest -> iedb
                                   \\-> structure -> geometry -> report
"""

__all__ = [
    "accessibility",
    "cleavage",
    "controls",
    "digest",
    "geometry",
    "iedb",
    "merops",
    "pda",
    "report",
    "structure",
]
