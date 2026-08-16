"""MEROPS cleavage-site specificity matrices as scoring models.

MEROPS publishes, per peptidase, how many observed cleavages carried each amino
acid at each of the eight substrate positions P4-P4'. Those counts become
log-odds position weight matrices.

The background matters more than usual here. De novo designs are famously
Ala/Glu/Lys-rich, so scoring them against a generic proteome background would
report enrichment that is an artifact of design composition rather than protease
preference. The background is therefore computed from the design corpus itself.

Schechter-Berger numbering: the scissile bond sits between P1 and P1'.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field

import httpx

from .paths import MEROPS_CACHE, ensure_caches

log = logging.getLogger(__name__)

PEPSUM_URL = "https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id={mid}"

POSITIONS = ("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
N_NONPRIME = 4  # P4..P1 precede the scissile bond

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# MEROPS identifiers for the endolysosomal proteases relevant to MHC class II
# antigen processing. `expect` is checked against the page so a wrong identifier
# cannot silently supply the wrong specificity model.
PEPTIDASES: dict[str, dict[str, str]] = {
    "cathepsin_S": {"merops": "C01.034", "expect": "cathepsin S", "family": "papain (Cys)"},
    "cathepsin_L": {"merops": "C01.032", "expect": "cathepsin L", "family": "papain (Cys)"},
    "cathepsin_B": {"merops": "C01.060", "expect": "cathepsin B", "family": "papain (Cys)"},
    "cathepsin_K": {"merops": "C01.036", "expect": "cathepsin K", "family": "papain (Cys)"},
    "cathepsin_D": {"merops": "A01.009", "expect": "cathepsin D", "family": "pepsin (Asp)"},
    "legumain": {"merops": "C13.004", "expect": "legumain", "family": "legumain (Cys)"},
}


@dataclass
class SpecificityMatrix:
    """Observed cleavage counts and the log-odds model derived from them."""

    name: str
    merops_id: str
    family: str
    counts: dict[str, dict[str, int]]  # position -> one-letter aa -> count
    n_cleavages: int
    cleavage_pattern: str = ""
    pwm: dict[str, dict[str, float]] = field(default_factory=dict)

    def column_total(self, pos: str) -> int:
        return sum(self.counts.get(pos, {}).values())

    def build_pwm(self, background: dict[str, float], pseudocount: float = 20.0) -> None:
        """Log2 odds of each residue at each position versus the background.

        The pseudocount is spread proportionally to the background, which keeps
        residues MEROPS never observed at a position from producing -inf while
        still letting strong depletions (Pro at cathepsin S P2, for instance)
        register as clearly negative.
        """
        self.pwm = {}
        for pos in POSITIONS:
            col = self.counts.get(pos, {})
            total = sum(col.values())
            if total == 0:
                self.pwm[pos] = dict.fromkeys(STANDARD_AA, 0.0)
                continue
            row: dict[str, float] = {}
            for aa in STANDARD_AA:
                bg = background.get(aa, 1.0 / 20.0)
                observed = col.get(aa, 0) + pseudocount * bg
                freq = observed / (total + pseudocount)
                row[aa] = math.log2(freq / bg) if bg > 0 else 0.0
            self.pwm[pos] = row

    def score(self, window: str) -> float:
        """Score an eight-residue P4-P4' window. Unknown residues contribute 0."""
        if len(window) != len(POSITIONS):
            raise ValueError(f"expected {len(POSITIONS)}-residue window, got {len(window)}")
        return sum(
            self.pwm.get(pos, {}).get(aa, 0.0) for pos, aa in zip(POSITIONS, window)
        )

    def max_score(self) -> float:
        return sum(max(self.pwm[p].values()) for p in POSITIONS if self.pwm.get(p))


def fetch_pepsum(merops_id: str, *, force: bool = False) -> str:
    ensure_caches()
    dest = MEROPS_CACHE / f"{merops_id}.html"
    if dest.exists() and dest.stat().st_size > 5000 and not force:
        return dest.read_text(encoding="utf-8", errors="replace")
    log.info("fetching MEROPS %s", merops_id)
    with httpx.Client(timeout=60.0, follow_redirects=True) as c:
        r = c.get(PEPSUM_URL.format(mid=merops_id))
        r.raise_for_status()
        dest.write_text(r.text, encoding="utf-8")
    return r.text


_MATRIX_TABLE = re.compile(r'<table[^>]*summary="matrix"[^>]*>(.*?)</table>', re.S | re.I)
_ROW = re.compile(r"<tr>(.*?)</tr>", re.S | re.I)
_CELL = re.compile(r"<t[dh][^>]*>(.*?)</t[dh]>", re.S | re.I)
_TAGS = re.compile(r"<[^>]+>")
_PATTERN = re.compile(
    r"Cleavage pattern.*?<td[^>]*>(.*?\(based on\s+([\d,]+)\s+cleavages\))\s*</td>", re.S | re.I
)


def _text(html: str) -> str:
    return _TAGS.sub("", html).replace("&nbsp;", " ").strip()


def parse_matrix(html: str, *, name: str, merops_id: str, family: str) -> SpecificityMatrix:
    """Pull the P4-P4' count table out of a MEROPS pepsum page."""
    m = _MATRIX_TABLE.search(html)
    if not m:
        raise ValueError(f"no specificity matrix table found for {merops_id}")

    counts: dict[str, dict[str, int]] = {p: {} for p in POSITIONS}
    header_seen = False
    for row_html in _ROW.findall(m.group(1)):
        cells = [_text(c) for c in _CELL.findall(row_html)]
        if not cells:
            continue
        if not header_seen:
            # The header row names the eight positions.
            if any(c.strip() == "P1'" for c in cells):
                header_seen = True
            continue
        aa3 = cells[0].strip()
        aa1 = THREE_TO_ONE.get(aa3)
        if aa1 is None:
            continue
        values = cells[1 : 1 + len(POSITIONS)]
        if len(values) < len(POSITIONS):
            continue
        for pos, raw in zip(POSITIONS, values):
            digits = raw.replace(",", "").strip()
            if digits.lstrip("-").isdigit():
                counts[pos][aa1] = int(digits)

    found = sum(len(v) for v in counts.values())
    if found < len(POSITIONS) * 15:
        raise ValueError(f"{merops_id}: only parsed {found} matrix cells, page layout changed?")

    pattern, n_cleavages = "", 0
    pm = _PATTERN.search(html)
    if pm:
        pattern = _text(pm.group(1))
        n_cleavages = int(pm.group(2).replace(",", ""))

    return SpecificityMatrix(
        name=name,
        merops_id=merops_id,
        family=family,
        counts=counts,
        n_cleavages=n_cleavages,
        cleavage_pattern=pattern,
    )


def verify_identity(html: str, expect: str, merops_id: str) -> None:
    """Guard against a wrong MEROPS identifier silently supplying a wrong model."""
    head = _text(html[:8000]).lower()
    if expect.lower() not in head:
        raise ValueError(
            f"MEROPS {merops_id} does not look like {expect!r}; refusing to use this matrix"
        )


def background_from_sequences(sequences: list[str]) -> dict[str, float]:
    """Amino-acid background frequencies from the corpus being analysed."""
    counts = dict.fromkeys(STANDARD_AA, 0)
    for seq in sequences:
        for aa in seq:
            if aa in counts:
                counts[aa] += 1
    total = sum(counts.values())
    if total == 0:
        return {aa: 1.0 / 20 for aa in STANDARD_AA}
    # Floor rare residues so a corpus lacking Trp cannot yield a divide-by-zero.
    return {aa: max(counts[aa] / total, 1e-4) for aa in STANDARD_AA}


def load_matrices(
    background: dict[str, float],
    *,
    names: list[str] | None = None,
    force: bool = False,
) -> dict[str, SpecificityMatrix]:
    """Fetch, verify, parse and calibrate every protease model."""
    out: dict[str, SpecificityMatrix] = {}
    for name, meta in PEPTIDASES.items():
        if names and name not in names:
            continue
        html = fetch_pepsum(meta["merops"], force=force)
        verify_identity(html, meta["expect"], meta["merops"])
        mat = parse_matrix(html, name=name, merops_id=meta["merops"], family=meta["family"])
        mat.build_pwm(background)
        out[name] = mat
        log.info(
            "%-14s %s  %6d cleavages  pattern=%s",
            name,
            meta["merops"],
            mat.n_cleavages,
            mat.cleavage_pattern or "-",
        )
    return out


def window_at(sequence: str, cut_index: int) -> str | None:
    """The P4-P4' window for a bond between ``cut_index-1`` and ``cut_index``.

    Returns None near the termini where a full eight-residue window does not fit,
    since a partial window is not comparable to a full one under the same matrix.
    """
    start = cut_index - N_NONPRIME
    end = cut_index + (len(POSITIONS) - N_NONPRIME)
    if start < 0 or end > len(sequence):
        return None
    return sequence[start:end]
