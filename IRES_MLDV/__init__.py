"""IRES_MLDV: a biophysics-based surrogate scorer for IRES activity.

Public API
----------
Scoring (see :mod:`IRES_MLDV.scorer`):
* :func:`find_splice_loop`     -- parallel scorer (per-sequence score + best site)
* :func:`mutation`             -- rewrite the 3' flanking acceptor context
* :data:`ANA_PRE/ANA_END/ORF/ORF_ADDED/IMPORT_SEQ` -- construct constants

I/O (see :mod:`IRES_MLDV.io_utils`):
* :func:`read_ires_csv` / :func:`load_to_csv` -- CSV helpers
"""

from __future__ import annotations

from .construct import (
    ANA_PRE,
    ANA_END,
    ORF,
    ORF_ADDED,
    IMPORT_SEQ,
    mutation,
)
from .scorer import (
    find_Hairpin_Loop,
    get_pair_prob,
    find_max_list,
    get_stem_prob_score,
    guess_loop,
    guess_pair,
    guess_ana_loop_num,
    compute_st_score,
    find_splice_loop,
)
from .io_utils import SCORE_COLUMNS, read_ires_csv, load_to_csv

__version__ = "1.0.0"

__all__ = [
    "ANA_PRE",
    "ANA_END",
    "ORF",
    "ORF_ADDED",
    "IMPORT_SEQ",
    "mutation",
    "find_Hairpin_Loop",
    "get_pair_prob",
    "find_max_list",
    "get_stem_prob_score",
    "guess_loop",
    "guess_pair",
    "guess_ana_loop_num",
    "compute_st_score",
    "find_splice_loop",
    "SCORE_COLUMNS",
    "read_ires_csv",
    "load_to_csv",
    "__version__",
]
