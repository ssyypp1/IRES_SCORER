"""I/O helpers for IRES_MLDV: read candidate IRES CSVs and write scored CSVs."""

from __future__ import annotations

import os
from typing import List, Sequence, Tuple

import pandas as pd

# Columns produced by :func:`load_to_csv` (one row per valid splice site).
SCORE_COLUMNS: List[str] = [
    "ires_source",
    "ires_seq",
    "intergrated-sequence",
    "ires_loop_seq",
    "splice-site",
    "loop_index",
    "pair-prob",
    "pair-mean-score",
    "pair-cnt",
    "loop_prob",
    "stem_prob",
    "t_pre_at_score",
    "t_behind_at_score",
    "ana_stem_prob",
    "mean_score",
]


def read_ires_csv(path: str) -> List[Tuple[str, str]]:
    """Read a CSV of candidate IRES sequences into ``(name, sequence)`` tuples.
    Accepts columns ``name,sequence`` or ``ires_name,seq`` (case-insensitive), or
    a single sequence column (names auto-generated)."""
    df = pd.read_csv(path)
    cols = {str(c).lower(): c for c in df.columns}

    name_col = None
    for cand in ("name", "ires_name", "id"):
        if cand in cols:
            name_col = cols[cand]
            break
    seq_col = None
    for cand in ("sequence", "seq"):
        if cand in cols:
            seq_col = cols[cand]
            break
    if seq_col is None and len(df.columns) == 1:
        seq_col = df.columns[0]

    if seq_col is None:
        raise ValueError(
            f"Could not find a sequence column in {path}; "
            f"expected one of 'sequence'/'seq' (got columns: {list(df.columns)})."
        )

    sequences = df[seq_col].astype(str).tolist()
    if name_col is not None:
        names = df[name_col].astype(str).tolist()
    else:
        names = [f"seq_{i:04d}" for i in range(len(sequences))]
    return list(zip(names, sequences))


def load_to_csv(res: Sequence[tuple], output_path: str = "mldv.csv") -> None:
    """Write splice-site info tuples to a CSV sorted by ``mean_score`` descending.
    Each tuple must have 15 fields matching :data:`SCORE_COLUMNS` (prepend the
    sequence name to the 14-field ``find_splice_loop`` info tuple)."""
    sorted_list = sorted(res, key=lambda x: -x[-1])
    rows = []
    for item in sorted_list:
        rows.append({col: item[idx] for idx, col in enumerate(SCORE_COLUMNS)})
    df = pd.DataFrame(rows, columns=SCORE_COLUMNS)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    df.to_csv(output_path, index=False)


__all__ = ["SCORE_COLUMNS", "read_ires_csv", "load_to_csv"]
