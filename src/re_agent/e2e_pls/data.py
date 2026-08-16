"""Track 1: dataset acquisition, tiling, clustering, and sampling.

Sources and licenses:
- Natural reference proteins: UniProt reviewed (Swiss-Prot) human entries,
  fetched via the UniProt REST API. UniProt data is CC BY 4.0.
- De novo designed proteins: RCSB PDB polymer entities matching a
  "de novo design" full-text search restricted to synthetic-construct
  source organism, fetched via the RCSB search + GraphQL data APIs. PDB
  structure/sequence data is public domain.
- DS613 (613 measured TAP-binding 9-mers, Diez-Rivero et al. 2010): no
  authoritative redistributable original was found. It is mirrored, with
  no declared license, at github.com/ChenWeiCCZU/CLTAP. Per the plan, this
  module downloads it locally only (gitignored data/raw/) and never
  commits the raw file -- only a provenance record goes in the dataset
  manifest (see dataset_card.py).
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

from re_agent.e2e_pls import schema

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
DS613_URL = "https://raw.githubusercontent.com/ChenWeiCCZU/CLTAP/main/Dataset/regression_DS613.csv"

DEFAULT_PEPTIDE_LEN = 9
DEFAULT_FLANK_LEN = 4
DEFAULT_HLA = "HLA-A*02:01"

STUDY_DE_NOVO = "rcsb_de_novo_v1"
STUDY_DS613 = "ds613"

LICENSE_NATURAL = "CC-BY-4.0 (UniProt)"
LICENSE_DE_NOVO = "CC0-1.0 (RCSB PDB, public domain)"
LICENSE_DS613 = "unverified-mirror -- see dataset_manifest.json; not redistributed"

# source_domain -> study_id, for every domain UniProt supplies (all share LICENSE_NATURAL)
STUDY_BY_DOMAIN = {
    "natural_human": "uniprot_human_reference_v1",
    "natural_viral": "uniprot_viral_reference_v1",
    "natural_bacterial": "uniprot_bacterial_reference_v1",
    "de_novo": STUDY_DE_NOVO,
}


@dataclass
class SourceProvenance:
    name: str
    url: str
    retrieved_at: str
    sha256: str
    license: str
    n_records: int
    notes: str = ""

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def now_iso() -> str:
    return datetime.now(UTC).isoformat()


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def _next_cursor(link_header: str | None) -> str | None:
    if not link_header:
        return None
    match = re.search(r"cursor=([^&>]+)", link_header)
    return match.group(1) if match else None


def _chunks(seq: list, size: int):
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


def _strip_expression_tags(seq: str) -> str:
    """Strip N/C-terminal His-tag purification handles (e.g. "MGSSHHHHHH...",
    "...LEHHHHHH") that are common in RCSB expression constructs but aren't
    part of the designed protein's biological sequence.
    """
    if m := re.search(r"H{4,}", seq[:25]):
        seq = seq[m.end() :]
    if m := re.search(r"H{4,}", seq[-25:]):
        cut = max(0, len(seq) - 25) + m.start()
        seq = seq[:cut]
    return seq


# --------------------------------------------------------------------------
# Source acquisition
# --------------------------------------------------------------------------

# Well-characterized human-pathogenic species, by UniProt organism_id. Swiss-Prot's
# manually-reviewed coverage is thin for most of these (a few dozen entries total),
# so viral/bacterial fetches also pull TrEMBL (unreviewed but real, submitted)
# sequences -- unlike the human fetch, which stays reviewed:true since human has
# ample reviewed coverage (20k+ entries).
VIRAL_ORGANISM_IDS = (
    11320,  # Influenza A virus
    11676,  # HIV-1
    2697049,  # SARS-CoV-2
    10359,  # Human cytomegalovirus
    10376,  # Epstein-Barr virus
    10407,  # Hepatitis B virus
    11103,  # Hepatitis C virus
    11234,  # Measles virus
    333760,  # Human papillomavirus type 16
    12637,  # Dengue virus
)
BACTERIAL_ORGANISM_IDS = (
    1773,  # Mycobacterium tuberculosis
    1280,  # Staphylococcus aureus
    562,  # Escherichia coli
    1639,  # Listeria monocytogenes
    28901,  # Salmonella enterica
    1314,  # Streptococcus pyogenes
    287,  # Pseudomonas aeruginosa
    487,  # Neisseria meningitidis
)


def _organism_clause(organism_ids: tuple[int, ...]) -> str:
    return "(" + " OR ".join(f"organism_id:{oid}" for oid in organism_ids) + ")"


def fetch_uniprot_proteins(
    n_proteins: int,
    organism_clause: str,
    min_length: int = 150,
    max_length: int = 500,
    reviewed: bool = True,
    seed: int = 0,
) -> dict[str, str]:
    """Sequences from UniProt matching `organism_clause` (a raw query fragment,
    e.g. "organism_id:9606" or an OR-group of several organism_ids).
    """
    reviewed_clause = " AND reviewed:true" if reviewed else ""
    query = f"{organism_clause}{reviewed_clause} AND length:[{min_length} TO {max_length}]"
    sequences: dict[str, str] = {}
    cursor: str | None = None
    with httpx.Client(timeout=30) as client:
        while len(sequences) < n_proteins * 2:  # overfetch a bit for shuffled sampling
            params = {"query": query, "format": "fasta", "size": 500}
            if cursor:
                params["cursor"] = cursor
            resp = client.get(UNIPROT_SEARCH_URL, params=params)
            resp.raise_for_status()
            records = _parse_fasta(resp.text)
            if not records:
                break
            for record_id, seq in records:
                accession = record_id.split("|")[1] if "|" in record_id else record_id
                sequences[f"uniprot_{accession}"] = seq
            cursor = _next_cursor(resp.headers.get("Link"))
            if not cursor:
                break
    rng = np.random.default_rng(seed)
    keys = list(sequences)
    rng.shuffle(keys)
    return {k: sequences[k] for k in keys[:n_proteins]}


def fetch_uniprot_human(n_proteins: int, seed: int = 0, **kwargs) -> dict[str, str]:
    return fetch_uniprot_proteins(
        n_proteins, "organism_id:9606", reviewed=True, seed=seed, **kwargs
    )


def fetch_uniprot_viral(n_proteins: int, seed: int = 0, **kwargs) -> dict[str, str]:
    return fetch_uniprot_proteins(
        n_proteins, _organism_clause(VIRAL_ORGANISM_IDS), reviewed=False, seed=seed, **kwargs
    )


def fetch_uniprot_bacterial(n_proteins: int, seed: int = 0, **kwargs) -> dict[str, str]:
    return fetch_uniprot_proteins(
        n_proteins, _organism_clause(BACTERIAL_ORGANISM_IDS), reviewed=False, seed=seed, **kwargs
    )


def fetch_rcsb_de_novo(
    n_proteins: int, min_length: int = 50, max_length: int = 300, seed: int = 0
) -> dict[str, str]:
    """Polymer entities annotated as de novo designs in RCSB PDB."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": "de novo design"},
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                        "operator": "exact_match",
                        "value": "synthetic construct",
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": min(5000, n_proteins * 6)}},
    }
    sequences: dict[str, str] = {}
    with httpx.Client(timeout=30) as client:
        resp = client.post(RCSB_SEARCH_URL, json=query)
        resp.raise_for_status()
        ids = [h["identifier"] for h in resp.json().get("result_set", [])]

        rng = np.random.default_rng(seed)
        rng.shuffle(ids)

        for batch in _chunks(ids, 50):
            if len(sequences) >= n_proteins:
                break
            gql = {
                "query": (
                    "query($ids: [String!]!) { polymer_entities(entity_ids: $ids) "
                    "{ rcsb_id entity_poly { pdbx_seq_one_letter_code_can } } }"
                ),
                "variables": {"ids": batch},
            }
            resp = client.post(RCSB_GRAPHQL_URL, json=gql)
            resp.raise_for_status()
            for entity in resp.json().get("data", {}).get("polymer_entities") or []:
                if entity is None:
                    continue
                raw_seq = entity["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "")
                seq = _strip_expression_tags(raw_seq)
                if (
                    len(seq) >= DEFAULT_PEPTIDE_LEN + 2 * DEFAULT_FLANK_LEN
                    and set(seq) <= schema.CANONICAL_RESIDUES
                ):
                    sequences[f"rcsb_{entity['rcsb_id']}"] = seq
    keys = list(sequences)[:n_proteins]
    return {k: sequences[k] for k in keys}


def fetch_ds613(dest_dir: str | Path) -> tuple[pd.DataFrame, SourceProvenance]:
    """Download DS613 into gitignored `dest_dir`; never committed to git."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_path = dest_dir / "ds613.csv"

    with httpx.Client(timeout=30) as client:
        resp = client.get(DS613_URL)
        resp.raise_for_status()
    dest_path.write_bytes(resp.content)
    sha256 = hashlib.sha256(resp.content).hexdigest()

    df = pd.read_csv(dest_path)
    lengths = df["Sequence"].str.len().unique()
    if len(df) != 613 or list(lengths) != [9]:
        raise ValueError(
            f"DS613 fetch doesn't match expected shape: {len(df)} rows, peptide lengths {lengths}"
        )

    provenance = SourceProvenance(
        name="DS613",
        url=DS613_URL,
        retrieved_at=now_iso(),
        sha256=sha256,
        license=LICENSE_DS613,
        n_records=len(df),
        notes="Diez-Rivero et al. 2010 TAP-binding peptides, mirrored from an unlicensed repo.",
    )
    return df, provenance


# --------------------------------------------------------------------------
# Tiling and candidate pool assembly
# --------------------------------------------------------------------------


def tile_protein(
    parent_id: str,
    sequence: str,
    source_domain: str,
    peptide_len: int = DEFAULT_PEPTIDE_LEN,
    flank_len: int = DEFAULT_FLANK_LEN,
) -> list[dict]:
    """9-mer windows with flanks. `start` runs from 1, not 0: a window at
    the protein's very first position has no upstream residue, so no valid
    N-terminal cleavage-site label exists for it (see label.py) -- simplest
    to exclude it here rather than carry a NaN special case downstream.
    """
    n = len(sequence)
    parent_hash = hashlib.sha256(sequence.encode()).hexdigest()[:16]
    rows = []
    for start in range(1, n - peptide_len + 1):
        end = start + peptide_len
        peptide = sequence[start:end]
        n_flank = sequence[max(0, start - flank_len) : start]
        c_flank = sequence[end : end + flank_len]
        if not set(peptide) <= schema.CANONICAL_RESIDUES:
            continue
        n_flank_ok = set(n_flank) <= schema.CANONICAL_RESIDUES
        c_flank_ok = set(c_flank) <= schema.CANONICAL_RESIDUES
        if not (n_flank_ok and c_flank_ok):
            continue
        rows.append(
            {
                "peptide": peptide,
                "length": peptide_len,
                "parent_sequence_id": parent_id,
                "parent_sequence_hash": parent_hash,
                "start": start,
                "end": end,
                "n_flank": n_flank,
                "c_flank": c_flank,
                "source_domain": source_domain,
            }
        )
    return rows


def build_candidate_pool(sources: dict[str, dict[str, str]]) -> pd.DataFrame:
    """`sources`: source_domain (e.g. "natural_human") -> {parent_id: sequence}."""
    rows: list[dict] = []
    for source_domain, seqs in sources.items():
        for pid, seq in seqs.items():
            rows.extend(tile_protein(pid, seq, source_domain))
    df = pd.DataFrame(rows)
    return df.drop_duplicates(subset="peptide", keep="first").reset_index(drop=True)


def quantile_stratified_sample(
    df: pd.DataFrame, target_n: int, score_col: str, n_bins: int = 10, seed: int = 0
) -> pd.DataFrame:
    """Sample ~evenly across `score_col` quantile bins so the corpus spans
    binders and non-binders rather than whatever a random protein sample
    happens to contain (mostly non-binders).
    """
    if len(df) <= target_n:
        return df.reset_index(drop=True)
    rng = np.random.default_rng(seed)
    bins = pd.qcut(df[score_col], q=n_bins, duplicates="drop")
    per_bin = max(1, target_n // bins.cat.categories.size)

    sampled_idx: list = []
    for _, group in df.groupby(bins, observed=True):
        n = min(len(group), per_bin)
        sampled_idx.extend(rng.choice(group.index, size=n, replace=False))

    if len(sampled_idx) < target_n:
        remaining_idx = df.index.difference(sampled_idx)
        extra_n = min(target_n - len(sampled_idx), len(remaining_idx))
        if extra_n > 0:
            sampled_idx.extend(rng.choice(remaining_idx, size=extra_n, replace=False))

    return df.loc[sampled_idx].sample(frac=1, random_state=seed).reset_index(drop=True)


# --------------------------------------------------------------------------
# Clustering and splits
# --------------------------------------------------------------------------


def _assign_splits(
    groups: np.ndarray, seed: int, ratios: tuple[float, float, float]
) -> dict[str, str]:
    unique = np.unique(groups)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(unique))
    n_train = int(len(unique) * ratios[0])
    n_val = int(len(unique) * ratios[1])
    split_of = np.empty(len(unique), dtype=object)
    split_of[order[:n_train]] = "train"
    split_of[order[n_train : n_train + n_val]] = "val"
    split_of[order[n_train + n_val :]] = "test"
    return dict(zip(unique, split_of, strict=True))


def assign_clusters_and_splits(
    df: pd.DataFrame, seed: int = 0, ratios: tuple[float, float, float] = (0.8, 0.1, 0.1)
) -> pd.DataFrame:
    """Whole parent-protein clusters go to a single split, so overlapping
    tiled windows from the same protein never leak across train/val/test.
    """
    df = df.copy()
    df["protein_cluster_id"] = df["parent_sequence_id"]
    df["peptide_cluster_id"] = df["parent_sequence_id"] + "_nbhd" + (df["start"] // 5).astype(str)
    split_map = _assign_splits(df["protein_cluster_id"].unique(), seed, ratios)
    df["split"] = df["protein_cluster_id"].map(split_map)
    return df


def build_ds613_rows(
    ds613_df: pd.DataFrame, seed: int = 0, ratios: tuple[float, float, float] = (0.8, 0.1, 0.1)
) -> pd.DataFrame:
    """DS613 peptides are isolated combinatorial-library measurements, not
    tiled from a shared parent protein, so each gets its own singleton
    cluster and a simple seeded random split (no overlap-window leakage is
    possible between them either way).
    """
    ds613_df = ds613_df.reset_index(drop=True)
    rows = []
    now = now_iso()
    for i, r in ds613_df.iterrows():
        peptide = r["Sequence"]
        rows.append(
            {
                "row_id": f"ds613:{i}",
                "peptide": peptide,
                "length": len(peptide),
                "parent_sequence_id": "ds613_synthetic_library",
                "parent_sequence_hash": "ds613",
                "start": 0,
                "end": len(peptide),
                "n_flank": "",
                "c_flank": "",
                "source_domain": "natural_human",
                "cleave_n_prob": np.nan,  # no real parent-protein flank context for Pepsickle
                "cleave_c_prob": np.nan,
                "tap_log_ic50_relative": float(r["log(IC50_relative)"]),
                "mhc_affinity_nm": np.nan,  # filled by label.py (peptide-only, doesn't need flanks)
                "mhc_percentile": np.nan,
                "hla_allele": DEFAULT_HLA,
                "label_origin": "measured",
                "label_model_version": "",
                "source_uri": DS613_URL,
                "license": LICENSE_DS613,
                "generated_at": now,
                "study_id": STUDY_DS613,
                "protein_cluster_id": f"ds613_{i}",
                "peptide_cluster_id": f"ds613_{i}",
                "encoder_model_id": schema.ENCODER_MODEL_ID,
                "pooling_recipe": "mean_9mer",
                "embedding_cache_key": schema.embedding_cache_key(peptide, "", ""),
                "has_structure_track": False,
            }
        )
    out = pd.DataFrame(rows)
    split_map = _assign_splits(out["row_id"].to_numpy(), seed, ratios)
    out["split"] = out["row_id"].map(split_map)
    return out


def finalize_candidate_columns(df: pd.DataFrame, hla_allele: str = DEFAULT_HLA) -> pd.DataFrame:
    """Fill the identity/provenance/ESM3 columns that don't depend on labels,
    for rows produced by `build_candidate_pool` + `assign_clusters_and_splits`.
    Label columns (cleave_*, tap_*, mhc_*) are left for label.py to fill in.
    """
    df = df.copy()
    now = now_iso()
    df["row_id"] = (
        df["parent_sequence_id"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    )
    df["hla_allele"] = hla_allele
    df["label_origin"] = "teacher"
    df["label_model_version"] = ""
    is_natural = df["source_domain"].str.startswith("natural_")
    df["license"] = np.where(is_natural, LICENSE_NATURAL, LICENSE_DE_NOVO)
    df["study_id"] = df["source_domain"].map(STUDY_BY_DOMAIN).fillna(STUDY_DE_NOVO)
    df["source_uri"] = np.where(
        is_natural,
        "https://www.uniprot.org/uniprotkb/"
        + df["parent_sequence_id"].str.removeprefix("uniprot_"),
        "https://www.rcsb.org/structure/"
        + df["parent_sequence_id"].str.removeprefix("rcsb_").str.split("_").str[0],
    )
    df["generated_at"] = now
    df["encoder_model_id"] = schema.ENCODER_MODEL_ID
    df["pooling_recipe"] = "mean_9mer"
    df["embedding_cache_key"] = [
        schema.embedding_cache_key(p, nf or "", cf or "")
        for p, nf, cf in zip(df["peptide"], df["n_flank"], df["c_flank"], strict=True)
    ]
    df["has_structure_track"] = False
    for col in (
        "cleave_n_prob",
        "cleave_c_prob",
        "tap_log_ic50_relative",
        "mhc_affinity_nm",
        "mhc_percentile",
    ):
        if col not in df.columns:
            df[col] = np.nan
    return df
