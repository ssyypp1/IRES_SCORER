# IRES_MLDV

A biophysics-based **surrogate scorer for IRES (Internal Ribosome Entry Site)
activity**. Given a candidate IRES sequence, it predicts how well the IRES can
drive downstream ORF expression by combining [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)
secondary-structure prediction with a six-component heuristic score built around
the hairpin-loop / splice-site context of the construct.

---

## Installation

```bash
pip install -r requirements.txt
```

Requirements: `ViennaRNA>=2.4`, `numpy>=1.20`, `pandas>=1.3`. Python 3.8+.

> If `pip install ViennaRNA` is unavailable, install it from your system package
> manager or <https://www.tbi.univie.ac.at/RNA/#download>.

---

## Quick start (CLI)

```bash
cd IRES_MLDV
python run_score.py --input data/example_ires.csv
```

Writes `result/scores.csv` (one row per valid splice site, sorted by
`mean_score` descending), with columns:

```
ires_source, ires_seq, intergrated-sequence, ires_loop_seq, splice-site,
loop_index, pair-prob, pair-mean-score, pair-cnt, loop_prob, stem_prob,
t_pre_at_score, t_behind_at_score, ana_stem_prob, mean_score
```

Use `--output` to change the path, `--verbose` to print the six sub-scores for
every valid splice site.

---

## Python API

```python
from IRES_MLDV import find_splice_loop, load_to_csv

seqs = ["TTAAAACAGC...", "CCCCTCTCCC..."]
ans, info = find_splice_loop(seqs, workers=10)
# ans[k]  = score for seqs[k]: max(ana_stable) over valid splice sites (0 if none)
# info[k] = 14-field tuple describing the best site for seqs[k] (or [])
```

Override the construct components to use a different reporter / ribozyme
backbone: `ana_pre` (upstream ribozyme half), `ana_end` (downstream ribozyme
half, carries the internal guide sequence), `orf` (reporter ORF), `orf_added`
(spacer after the ORF).

---

## Scoring principle

This scorer evaluates a **permuted intron–exon (PIE) circularization
construct**: a group-I self-splicing ribozyme is split into two halves
(`ana_pre`, `ana_end`) and placed around the candidate IRES so that, after
*trans*-splicing, the two ends of the insert are joined into a circle. For
every candidate IRES `sequence`, the scorer:

1. **Folds** `sequence` with ViennaRNA (MFE) and parses hairpin loops.
2. **Keeps** loops whose loop region is unpaired (mean pair probability `<= 0.1`)
   and whose flanking stem is well-paired (stem-prob score `>= 0.7`).
3. **Enumerates candidate splice sites** inside each kept loop — a `T` signal
   followed by enough downstream bases (the catalytic splice-site rule).
4. **Rewrites the internal guide sequence (IGS)** at the head of `ana_end` to
   pair with the sequence upstream of the candidate site, then builds the PIE
   construct for splice site `j`:

   ```
   all_seq = ana_pre + sequence[j:] + orf + orf_added + sequence[:j] + mutation_ana_end
   ```

5. **Validates** that the ribozyme halves still fold into an active
   conformation (`>= 4` hairpin loops at the 5' flank, `>= 2` at the 3' flank).
6. **Computes six sub-scores** and averages them into `mean_score`:

   | sub-score | meaning |
   | --- | --- |
   | `pair-mean-score` | mean IGS-target pairing probability at the splice junction |
   | `loop_prob` (`1-p`) | loop openness (1 − mean loop pair probability) |
   | `stem_prob` | stem stability of the hairpin loop |
   | `t_pre_at_score` | AT-content score just upstream of the splice site |
   | `t_behind_at_score` | AT fraction of the rewritten IGS |
   | `ana_stem_prob` | normalized stability of the ribozyme-half stems |

The per-sequence score is `max(ana_stable)` over valid splice sites
(`mean_score` is recorded in the info tuple but not used as the score).

---

## Construct layout (default = EGFP reporter)

```
ana_pre + IRES[j:] + orf + orf_added + IRES[:j] + mutation_ana_end
  122nt     ...    EGFP    50nt       ...        129nt
```

* `ana_pre` — upstream ribozyme half (122 nt)
* `ana_end` — downstream ribozyme half (129 nt); its first 11 nt are the
  internal guide sequence (IGS), rewritten per splice site. The
  `guess_ana_loop_num` check probes the last **120** nt of this half.
* `orf` — EGFP coding sequence (720 nt)
* `orf_added` — poly-A-like spacer after the ORF (50 nt)

---

## Repository layout

```
IRES_MLDV/
├── README.md
├── requirements.txt
├── run_score.py            # CLI: score candidate IRES sequences
├── data/
│   ├── example_ires.csv    # example viral IRES sequences
│   └── algo.tex / algo.pdf # algorithm pseudocode (paper figure)
├── result/                 # scoring outputs (gitignored)
└── IRES_MLDV/               # importable package
    ├── __init__.py
    ├── construct.py         # construct constants + IGS rewriting
    ├── scorer.py            # scoring functions (ViennaRNA-based)
    └── io_utils.py          # read_ires_csv / load_to_csv
```
---

