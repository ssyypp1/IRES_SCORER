#!/usr/bin/env python3
"""CLI: score candidate IRES sequences from a CSV file.

Example
-------
    python run_score.py --input data/example_ires.csv
    python run_score.py --input data/example_ires.csv --output result/scores.csv --verbose

Input CSV accepts columns ``name,sequence`` or ``ires_name,seq``
(case-insensitive). Output (default ``result/scores.csv``) has one row per
valid splice site, sorted by ``mean_score`` descending, with columns::

    ires_source, ires_seq, intergrated-sequence, ires_loop_seq, splice-site,
    loop_index, pair-prob, pair-mean-score, pair-cnt, loop_prob, stem_prob,
    t_pre_at_score, t_behind_at_score, ana_stem_prob, mean_score

Uses :func:`find_splice_loop` (one row per sequence, the best valid splice site).
"""

from __future__ import annotations

import argparse
import sys

from IRES_MLDV.scorer import find_splice_loop
from IRES_MLDV.io_utils import read_ires_csv, load_to_csv


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Score candidate IRES sequences with the IRES_MLDV surrogate.",
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input CSV of candidate IRES sequences (columns: name,sequence or ires_name,seq).",
    )
    parser.add_argument(
        "--output", "-o", default="result/scores.csv",
        help="Output CSV path (one row per valid splice site, sorted by mean_score desc; "
             "default result/scores.csv).",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print the six sub-scores for every valid splice site.",
    )
    args = parser.parse_args(argv)

    data = read_ires_csv(args.input)
    print(f"[IRES_MLDV] Loaded {len(data)} candidate IRES sequences from {args.input}")

    names = [n for n, _ in data]
    seqs = [s for _, s in data]
    ans, info = find_splice_loop(seqs, verbose=args.verbose)

    # Prepend name to each 14-field info tuple → 15-field rows for load_to_csv.
    rows = [(names[k],) + tuple(info[k]) for k in range(len(data)) if info[k]]
    load_to_csv(rows, args.output)

    n_valid = sum(1 for a in ans if a > 0)
    print(f"[IRES_MLDV] {n_valid}/{len(data)} sequences had >=1 valid splice site")
    print(f"[IRES_MLDV] {len(rows)} rows written to {args.output}")
    if ans:
        best_idx = max(range(len(ans)), key=lambda k: ans[k])
        print(f"[IRES_MLDV] Best sequence: {names[best_idx]}  score={ans[best_idx]:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
