# IRES_MLDV

Given a candidate IRES (Internal Ribosome Entry Site) sequence, this tool
searches for hairpin stem-loops that satisfy a set of structural criteria and,
within each qualifying loop, identifies a suitable catalytic splice site for
ribozyme-mediated circularization, using [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) secondary-structure prediction.


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
`mean_score` descending). Use `--output` to change the path, `--verbose` to
print the six sub-scores for every valid splice site.

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
    ├── scorer.py             # scoring functions (ViennaRNA-based)
    └── io_utils.py           # read_ires_csv / load_to_csv
```
