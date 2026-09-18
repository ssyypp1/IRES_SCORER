"""Sequence construction constants and helpers for IRES_MLDV.

Fixed flanking / reporter sequences used to build the integrated construct
scored by :mod:`IRES_MLDV.scorer`, plus the :func:`mutation` helper that
rewrites the 5' end of the 3' flank based on a candidate splice site. Defaults
reproduce the EGFP reporter construct; override them in scorer calls to use a
different reporter / vector backbone.

Construct layout::

    ana_pre + IRES[j:] + orf + orf_added + IRES[:j] + mutation_ana_end
"""

from __future__ import annotations

from typing import Dict

# 5' flanking / upstream backbone (ana_pre, 122 nt).
ANA_PRE: str = (
    "TGACTTACAACTAATCGGAAGGTGCAGAGACTCGACGGGAGCTACCCTAACGTCAAGACGAGGGTAAAGA"
    "GAGAGTCCAATTCTCAAAGCCAATAGGCAGTAGCGAAAGCTGCAAGAGAATG"
)

# 3' flanking / downstream backbone (ana_end, 129 nt). The first 11 nt
# (AAATAATTGAG) are rewritten by :func:`mutation` per splice site.
ANA_END: str = (
    "AAATAATTGAGCCTTAAAGAAGAAATTCTTTAAGTGGATGCTCTCAAACTCAGGGAAACCTAAATCTAGT"
    "TATAGACAAGGCAATCCTGAGCCAAGCCGAAGTAGTAATTAGTAAGTCAACAATAGATG"
)

# Reporter ORF: EGFP coding sequence (720 nt).
ORF: str = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAA"
    "ACGGCCACAAGTTCAGCGTGTCTGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTT"
    "CATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAG"
    "TGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACG"
    "TCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGG"
    "CGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCAC"
    "AAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGG"
    "CGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC"
    "CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAA"
    "GACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGC"
    "ATGGACGAGCTGTACAAGTAA"
)

# poly-A-like spacer after the ORF (50 nt).
ORF_ADDED: str = "AAAAACAAAAAACAAAAAAAACAAAAAAAAAACCAAAAAAACAAAACACA"

# Conserved import sequence located at the 5' end of :data:`ANA_END`.
IMPORT_SEQ: str = "AAATAATTGAG"

# Watson-Crick complement table used by :func:`mutation`.
_COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def mutation(loop_seq: str, splice_site: int, ana_end: str = ANA_END) -> str:
    """Rewrite the first 11 nt of ``ana_end`` (the conserved AAATAATTGAG motif)
    as the complement of the loop context around ``splice_site``, so the 3'
    flank forms a compatible acceptor context. The rest of ``ana_end`` is
    unchanged. ``splice_site`` is the offset of the T splice signal within
    ``loop_seq``. Returns the rewritten flank (same length as ``ana_end``)."""
    important_seq = list(IMPORT_SEQ)
    important_seq[-1] = _COMPLEMENT[loop_seq[splice_site - 1 - 2]]  # 2nd base before T
    important_seq[-2] = _COMPLEMENT[loop_seq[splice_site - 1 - 1]]  # 1st base before T
    important_seq[-4] = _COMPLEMENT[loop_seq[splice_site]]          # 1st base after T
    important_seq[-5] = _COMPLEMENT[loop_seq[splice_site + 1]]      # 2nd base after T
    important_seq[0] = _COMPLEMENT[important_seq[-4]]
    important_seq[1] = _COMPLEMENT[important_seq[-5]]
    new_head = "".join(important_seq)
    return new_head + ana_end[len(new_head):]


__all__ = [
    "ANA_PRE",
    "ANA_END",
    "ORF",
    "ORF_ADDED",
    "IMPORT_SEQ",
    "mutation",
]
