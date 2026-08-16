"""Dataset construction: labeled natural peptides (IEDB) and unlabeled de novo windows.

Every prediction unit is an isolated 15-mer. MHC class II assays are run on
synthetic peptides, so embedding each window in isolation (rather than slicing it
out of a parent-protein embedding) keeps training and inference on the same footing.
"""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from re_agent.immuno.config import AA_SET, DATA_RAW, PATHS, WINDOW, ensure_dirs

IEDB_URL = "https://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip"
IEDB_ZIP = DATA_RAW / "tcell_full_v3.zip"
IEDB_CSV = DATA_RAW / "tcell_full_v3.csv"

# Column offsets in the two-row-header IEDB export (data starts at row 3).
IEDB_COLS = {
    10: "object_type",
    11: "peptide",
    19: "source_molecule",
    23: "source_organism",
    43: "host",
    122: "qualitative",
    141: "mhc_allele",
    145: "mhc_class",
}

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
RFDIFF_FASTA = (
    "https://raw.githubusercontent.com/adaptyvbio/rfdiff_il7ra/main/"
    "original_data/IL7Ra_binders_sequences.fasta"
)
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

MIN_PEPTIDE = 13
MAX_PEPTIDE = 21

# Poly-His purification tags and their usual linkers are construct artifacts, not
# part of the designed protein, so they must not leak into the de novo pool.
_TAG = re.compile(r"(H{6,}|^MGSSHHHHHHSSGLVPRGSH|^GSHM|^MGS|GSGSGS+|LEHHHHHH$)")


def _http_json(url: str, payload: dict | None = None, timeout: int = 120) -> dict:
    data = json.dumps(payload).encode() if payload is not None else None
    headers = {"Content-Type": "application/json"} if data else {}
    req = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def clean_sequence(seq: str) -> str:
    """Uppercase, drop expression tags, and reject anything non-canonical."""
    if not isinstance(seq, str):
        return ""
    seq = seq.strip().upper()
    seq = _TAG.sub("", seq)
    return seq if seq and set(seq) <= AA_SET else ""


def windows(seq: str, size: int = WINDOW, stride: int = 1) -> list[tuple[int, str]]:
    """Tile a sequence into (start_index, subsequence) windows.

    Sequences shorter than `size` are returned whole; the model pads and masks them.
    """
    if len(seq) <= size:
        return [(0, seq)] if seq else []
    return [(i, seq[i : i + size]) for i in range(0, len(seq) - size + 1, stride)]


# --------------------------------------------------------------------------- #
# Labeled natural peptides (IEDB class II T-cell assays)
# --------------------------------------------------------------------------- #


def download_iedb(force: bool = False) -> Path:
    ensure_dirs()
    if IEDB_CSV.exists() and not force:
        return IEDB_CSV
    if not IEDB_ZIP.exists() or force:
        print(f"downloading IEDB T-cell export -> {IEDB_ZIP}")
        urllib.request.urlretrieve(IEDB_URL, IEDB_ZIP)
    with zipfile.ZipFile(IEDB_ZIP) as zf:
        zf.extractall(DATA_RAW)
    return IEDB_CSV


def parse_iedb(csv_path: Path | None = None, chunksize: int = 200_000) -> pd.DataFrame:
    """Collapse class II human T-cell assay rows into one labeled row per peptide."""
    csv_path = csv_path or download_iedb()
    idx = sorted(IEDB_COLS)
    names = [IEDB_COLS[i] for i in idx]

    frames = []
    for chunk in pd.read_csv(
        csv_path,
        skiprows=2,
        header=None,
        usecols=idx,
        names=names,
        chunksize=chunksize,
        low_memory=False,
        dtype=str,
    ):
        keep = (
            (chunk["object_type"] == "Linear peptide")
            & (chunk["mhc_class"] == "II")
            & chunk["host"].fillna("").str.contains("Homo sapiens", regex=False)
        )
        sub = chunk.loc[keep, ["peptide", "source_molecule", "source_organism", "qualitative"]]
        if len(sub):
            frames.append(sub)

    df = pd.concat(frames, ignore_index=True)
    df["peptide"] = df["peptide"].map(clean_sequence)
    df = df[df["peptide"].str.len().between(MIN_PEPTIDE, MAX_PEPTIDE)]

    qual = df["qualitative"].fillna("")
    df = df[qual.str.startswith(("Positive", "Negative"))].copy()
    df["is_pos"] = df["qualitative"].str.startswith("Positive").astype(int)

    def first_value(series: pd.Series) -> str:
        clean = series.dropna()
        return clean.iloc[0] if len(clean) else ""

    grouped = df.groupby("peptide", sort=False)
    out = grouped.agg(
        n_pos=("is_pos", "sum"),
        n_assay=("is_pos", "size"),
        source_molecule=("source_molecule", first_value),
        source_organism=("source_organism", first_value),
    ).reset_index()
    out["n_neg"] = out["n_assay"] - out["n_pos"]
    # A peptide counts as immunogenic if it ever drove a response; a single positive
    # assay is real evidence, whereas a negative only says "not in this donor set".
    out["label"] = (out["n_pos"] > 0).astype(int)
    return out


def build_labeled(force: bool = False) -> pd.DataFrame:
    """Peptide table -> 15-mer window table with group keys for leakage-free splits."""
    ensure_dirs()
    if PATHS.labeled.exists() and not force:
        return pd.read_parquet(PATHS.labeled)

    peptides = parse_iedb()
    rows = []
    for rec in peptides.itertuples(index=False):
        for start, win in windows(rec.peptide, stride=1):
            rows.append(
                {
                    "seq": win,
                    "parent": rec.peptide,
                    "start": start,
                    "label": rec.label,
                    "n_pos": rec.n_pos,
                    "n_neg": rec.n_neg,
                    "source_molecule": rec.source_molecule or rec.peptide,
                    "source_organism": rec.source_organism or "unknown",
                }
            )
    df = pd.DataFrame(rows)
    # Group by source protein so overlapping peptides from one antigen cannot be
    # split across train and test.
    df["group"] = df["source_molecule"].where(df["source_molecule"].str.len() > 0, df["parent"])
    df.to_parquet(PATHS.labeled, index=False)
    print(f"labeled windows: {len(df)} ({df.label.mean():.1%} positive) -> {PATHS.labeled}")
    return df


# --------------------------------------------------------------------------- #
# Unlabeled de novo sequences
# --------------------------------------------------------------------------- #


def fetch_rcsb_denovo_ids(limit: int = 5000) -> list[str]:
    """PDB entities classified DE NOVO PROTEIN: genuinely designed, not natural."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "struct_keywords.pdbx_keywords",
                "operator": "contains_phrase",
                "value": "DE NOVO PROTEIN",
            },
        },
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    data = _http_json(RCSB_SEARCH, payload)
    return [hit["identifier"] for hit in data.get("result_set", [])]


def fetch_rcsb_sequences(entity_ids: list[str], batch: int = 200) -> dict[str, str]:
    out: dict[str, str] = {}
    for i in range(0, len(entity_ids), batch):
        chunk = entity_ids[i : i + batch]
        ids = ",".join(f'"{e}"' for e in chunk)
        query = (
            f"{{polymer_entities(entity_ids:[{ids}])"
            "{rcsb_id entity_poly{pdbx_seq_one_letter_code_can type}}}"
        )
        try:
            data = _http_json(RCSB_GRAPHQL, {"query": query})
        except OSError as exc:
            print(f"  rcsb batch {i} failed: {exc}")
            continue
        for ent in data.get("data", {}).get("polymer_entities") or []:
            poly = ent.get("entity_poly") or {}
            if poly.get("type") != "polypeptide(L)":
                continue
            seq = clean_sequence(poly.get("pdbx_seq_one_letter_code_can", ""))
            if len(seq) >= WINDOW:
                out[ent["rcsb_id"]] = seq
    return out


def load_rfdiffusion_binders() -> dict[str, str]:
    """95 RFdiffusion/ProteinMPNN IL-7Ra minibinders — the demo-facing design set."""
    path = DATA_RAW / "IL7Ra_binders_sequences.fasta"
    if not path.exists():
        urllib.request.urlretrieve(RFDIFF_FASTA, path)
    return read_fasta(path)


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name and buf:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name and buf:
        seqs[name] = "".join(buf)
    return {k: clean_sequence(v) for k, v in seqs.items() if clean_sequence(v)}


def build_unlabeled(
    stride: int = 2, max_windows: int = 120_000, force: bool = False
) -> pd.DataFrame:
    ensure_dirs()
    if PATHS.unlabeled.exists() and not force:
        return pd.read_parquet(PATHS.unlabeled)

    proteins: dict[str, tuple[str, str]] = {}
    for name, seq in load_rfdiffusion_binders().items():
        proteins[f"rfdiff:{name}"] = (seq, "rfdiffusion")
    print("querying RCSB for de novo designed entities...")
    ids = fetch_rcsb_denovo_ids()
    print(f"  {len(ids)} entities; fetching sequences")
    for rid, seq in fetch_rcsb_sequences(ids).items():
        proteins[f"pdb:{rid}"] = (seq, "pdb_de_novo")

    rows = []
    for pid, (seq, src) in proteins.items():
        for start, win in windows(seq, stride=stride):
            rows.append({"seq": win, "parent": pid, "start": start, "source": src})
    df = pd.DataFrame(rows).drop_duplicates(subset="seq").reset_index(drop=True)
    if len(df) > max_windows:
        df = df.sample(max_windows, random_state=0).reset_index(drop=True)
    df.to_parquet(PATHS.unlabeled, index=False)
    print(f"unlabeled de novo windows: {len(df)} from {len(proteins)} designs -> {PATHS.unlabeled}")
    return df


# --------------------------------------------------------------------------- #
# Natural reference cohort (percentile scale)
# --------------------------------------------------------------------------- #


def fetch_uniprot(query: str, limit: int) -> dict[str, str]:
    url = (
        f"{UNIPROT_STREAM}?format=fasta&size={min(limit, 500)}"
        f"&query={urllib.parse.quote(query)}"
    )
    with urllib.request.urlopen(url, timeout=180) as resp:
        text = resp.read().decode()
    seqs, name, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name and buf:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name and buf:
        seqs[name] = "".join(buf)
    return {k: clean_sequence(v) for k, v in seqs.items() if clean_sequence(v)}


def build_reference(
    stride: int = 2, per_panel: int = 120, force: bool = False, seed: int = 0
) -> pd.DataFrame:
    """A mixed natural panel: human self-proteins plus microbial foreign proteins.

    Spanning both ends of the tolerance spectrum keeps the risk percentile
    interpretable instead of saturating every de novo binder at the top. Whole
    proteins are kept intact — a protein-level percentile needs complete tilings.
    """
    ensure_dirs()
    if PATHS.reference.exists() and not force:
        return pd.read_parquet(PATHS.reference)

    panels = {
        "human": "reviewed:true AND organism_id:9606 AND length:[80 TO 400]",
        "bacterial": "reviewed:true AND organism_id:83333 AND length:[80 TO 400]",
        "pathogen": "reviewed:true AND organism_id:83332 AND length:[80 TO 400]",
    }
    rng = np.random.default_rng(seed)
    rows = []
    for panel, query in panels.items():
        try:
            seqs = fetch_uniprot(query, per_panel)
        except OSError as exc:
            print(f"  uniprot {panel} failed: {exc}")
            continue
        names = sorted(seqs)
        if len(names) > per_panel:
            names = [names[i] for i in rng.choice(len(names), per_panel, replace=False)]
        print(f"  reference panel {panel}: {len(names)} proteins")
        for name in names:
            for start, win in windows(seqs[name], stride=stride):
                rows.append({"seq": win, "parent": name, "start": start, "source": panel})

    df = pd.DataFrame(rows).reset_index(drop=True)
    df.to_parquet(PATHS.reference, index=False)
    print(
        f"reference windows: {len(df)} from {df.parent.nunique()} proteins -> {PATHS.reference}"
    )
    return df


# --------------------------------------------------------------------------- #
# Splits
# --------------------------------------------------------------------------- #


def group_split(
    groups: pd.Series, fracs: tuple[float, float, float] = (0.7, 0.15, 0.15), seed: int = 0
) -> pd.Series:
    """Assign whole groups to train/val/test so no source protein spans two splits."""
    rng = np.random.default_rng(seed)
    uniq = np.asarray(groups.unique(), dtype=object)
    rng.shuffle(uniq)
    n = len(uniq)
    n_train = int(fracs[0] * n)
    n_val = int(fracs[1] * n)
    assignment = {}
    for i, g in enumerate(uniq):
        if i < n_train:
            assignment[g] = "train"
        elif i < n_train + n_val:
            assignment[g] = "val"
        else:
            assignment[g] = "test"
    return groups.map(assignment)


def main() -> None:
    build_labeled()
    build_unlabeled()
    build_reference()


if __name__ == "__main__":
    main()
