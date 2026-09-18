"""IRES surrogate scoring functions.

Biophysics-based scorer that estimates the ability of a candidate IRES sequence
to drive downstream ORF expression, using ViennaRNA for secondary-structure
prediction (MFE + partition function).
"""

from __future__ import annotations

import concurrent.futures
import statistics
from collections import Counter
from typing import List, Sequence, Tuple

import numpy as np
import RNA

from .construct import ANA_PRE, ANA_END, ORF, ORF_ADDED, mutation


def find_Hairpin_Loop(seq: str, structure: str) -> Tuple[List[str], List[List[int]]]:
    """Parse hairpin loops from a dot-bracket secondary structure.

    Returns the loop sub-sequences (including flanking stem base-pairs) and
    their 1-based ``[start, end]`` indices into ``seq``.
    """
    i, pre, ans, index = 0, 0, [], []
    flag = False
    while i < len(structure):
        if structure[i] == "(":
            pre += 1
            flag = False
        elif structure[i] == ")":
            pre -= 1
            flag = False
        else:
            s = ""
            flag = False
            if pre != 0 and i - 1 >= 0 and structure[i - 1] == "(":
                flag = True
                s += seq[i - 1]
                original_index = i - 1
                while i < len(structure) and structure[i] == ".":
                    s += seq[i]
                    i += 1
            if i < len(structure) and structure[i] == ")" and len(s) > 1:
                s += seq[i]
                ans.append(s)
                index.append([original_index + 1, i - 1])
        if not flag:
            i += 1
    return ans, index


def get_pair_prob(
    sequence: str, sitea: int = -1, siteb: int = -1
) -> Tuple[List[float], List[float], int]:
    """Per-base maximum pairing probability via ViennaRNA partition function.

    When ``sitea``/``siteb`` are set, the stem-probability list is restricted to
    bases outside ``[sitea, siteb]`` (the flanking region of an integrated
    construct). Returns ``(per-base max prob, sorted stem probs (top 30), stem
    count)``.
    """
    RNA.cvar.temperature = 37
    fc = RNA.fold_compound(sequence)
    fc.pf()

    plist = fc.plist_from_probs(cutoff=0)

    length = len(sequence)

    matrix = [[0.0 for _ in range(length)] for _ in range(length)]

    for pair in plist:
        matrix[pair.i - 1][pair.j - 1] = pair.p
        matrix[pair.j - 1][pair.i - 1] = pair.p
    res = []
    my_set = set()
    stem_prob = []
    for i, row in enumerate(matrix):
        if max(row) in my_set and sitea != -1 and (i < sitea or i > siteb):
            stem_prob.append(max(row))
        elif sitea != -1 and (i < sitea or i > siteb):
            my_set.add(max(row))
        res.append(max(row))
    stem_num = len(stem_prob)
    stem_prob = sorted(stem_prob)[:30]
    return res, stem_prob, stem_num


def find_max_list(index: Sequence[Sequence[int]]) -> List[int]:
    """Return the ``[start, end]`` pair with the largest span."""
    ans, res = 0, []
    for i in range(len(index)):
        if index[i][1] - index[i][0] > ans:
            res = index[i]
            ans = index[i][1] - index[i][0]
    return res


def get_stem_prob_score(all_pair_prob: Sequence[float], stem_pair_prob: Sequence[float]) -> float:
    """Stem stability score: reward paired bases (prob appears twice), penalize unpaired."""
    number_counts = Counter(all_pair_prob)
    res = 0
    for _, num in enumerate(stem_pair_prob):
        if number_counts[num] == 2:
            res += num
        else:
            res += 1 - num
    return res / len(stem_pair_prob)


def guess_loop(seq: str) -> List[Tuple[List[int], float, float]]:
    """Locate candidate functional hairpin loops in ``seq``.

    A loop is kept if its loop region is unpaired (mean pair prob ``<= 0.1``)
    and its flanking stem is well-paired (stem prob score ``>= 0.7``). Returns
    ``(loop_index, mean_loop_pair_prob, stem_prob_score)`` tuples.
    """
    res_loop = []
    RNA.cvar.temperature = 37
    fc = RNA.fold_compound(seq)
    structure, value = fc.mfe()
    _, index = find_Hairpin_Loop(seq, structure)
    all_pair_prob, _, _ = get_pair_prob(seq)
    for _, max_loop in enumerate(index):
        loop_pair_prob, stem_pair_prob = [], []
        for i in range(max(max_loop[0], 0), min(len(seq) - 1, max_loop[1] + 1)):
            loop_pair_prob.append(all_pair_prob[i])
        for i in range(max(max_loop[0] - 8, 0), min(len(seq) - 1, max_loop[0])):
            stem_pair_prob.append(all_pair_prob[i])
        for i in range(max(max_loop[1] + 1, 0), min(len(seq) - 1, max_loop[1] + 9)):
            stem_pair_prob.append(all_pair_prob[i])
        if statistics.mean(loop_pair_prob) <= 0.1 and get_stem_prob_score(all_pair_prob, stem_pair_prob) >= 0.7:
            res_loop.append(
                (max_loop, statistics.mean(loop_pair_prob), get_stem_prob_score(all_pair_prob, stem_pair_prob))
            )
    return res_loop


def guess_pair(all_pair_prob: Sequence[float]) -> List[float]:
    """Pairing probabilities at 5 fixed positions of the 3' acceptor region."""
    tmp = []
    for i in range(-129 - 3, -127):
        tmp.append(all_pair_prob[i])
    return tmp


def guess_ana_loop_num(seq: str, ana_pre: str = ANA_PRE, ana_end: str = ANA_END) -> bool:
    """Check the construct forms ``>= 4`` hairpin loops at the 5' flank and
    ``>= 2`` at the 3' flank (last 120 nt)."""
    RNA.cvar.temperature = 37
    fc = RNA.fold_compound(seq)
    structure, value = fc.mfe()
    _, index1 = find_Hairpin_Loop(seq[: len(ana_pre)], structure[: len(ana_pre)])
    _, index2 = find_Hairpin_Loop(seq[-120:], structure[-120:])
    if len(index1) >= 4 and len(index2) >= 2:
        return True
    return False


def compute_st_score(t_pre: str, t_behenid: str) -> Tuple[float, float]:
    """AT-content scores around the splice site. Returns ``(pre_dict[AT%], AT%
    of the rewritten 3' flank head)`` where ``pre_dict = {1.0: 0.5, 0.5: 1.0,
    0: 0.3}``."""
    cnt_pre, cnt_behenid = 0, 0
    for i in range(len(t_pre)):
        if t_pre[i] == "A" or t_pre[i] == "T":
            cnt_pre += 1
    for i in range(len(t_behenid)):
        if t_behenid[i] == "A" or t_behenid[i] == "T":
            cnt_behenid += 1
    pre_percent, behenid_percent = cnt_pre / len(t_pre), cnt_behenid / len(t_behenid)
    pre_dict = {1.0: 0.5, 0.5: 1.0, 0: 0.3}
    return pre_dict[pre_percent], behenid_percent


def _pair_count(tmp: Sequence[float]) -> int:
    """Count positions with pair probability ``>= 0.7``."""
    cnt = 0
    for i in range(len(tmp)):
        if tmp[i] >= 0.7:
            cnt += 1
    return cnt


def find_splice_loop(
    data: Sequence[str],
    ana_pre: str = ANA_PRE,
    ana_end: str = ANA_END,
    orf: str = ORF,
    orf_added: str = ORF_ADDED,
    workers: int = 10,
    verbose: bool = False,
) -> Tuple[List[float], List[list]]:
    """Score a list of candidate IRES sequences (parallel).

    For each sequence returns ``(score, best_site_info)`` where ``score =
    max(ana_stable)`` over valid splice sites (``0`` if none) and
    ``best_site_info`` is the 14-field tuple of the best site (or ``[]``):

    ``(sequence, all_seq, loop_seq, splice_site, loop_index, pair_prob_list,
      mean(pair_prob), pair_count, 1-loop_prob, stem_prob, t_pre_at_score,
      t_behind_at_score, ana_stable, mean_score)``

    A splice site is valid when the preceding base is ``T``, ``splice_site >= 4``,
    ``len(loop_seq) - splice_site >= 3``, the construct passes
    :func:`guess_ana_loop_num`, and :func:`guess_pair` yields ``>= 3`` positions
    with pair prob ``>= 0.7``. ``verbose`` prints the six sub-scores per site.
    """
    RNA.cvar.temperature = 37

    def score_one(sequence: str):
        tmp_ans, res = [], []
        index = guess_loop(sequence)
        for _, t in enumerate(index):
            indexs, loop_prob, stem_prob = t
            for j in range(indexs[0] + 1, indexs[1]):
                splice_site = j - indexs[0]
                loop_seq = sequence[indexs[0]:indexs[1] + 1]
                if not (
                    splice_site >= 1
                    and loop_seq[splice_site - 1] == "T"
                    and splice_site >= 4
                    and len(loop_seq) - splice_site >= 3
                ):
                    continue
                mutation_ana_end = mutation(loop_seq, splice_site, ana_end)
                all_seq = (
                    ana_pre + sequence[j:] + orf + orf_added + sequence[:j] + mutation_ana_end
                )
                t_pre_at_score, t_behind_at_score = compute_st_score(
                    sequence[j - 3:j - 1], mutation_ana_end[:2]
                )
                if not guess_ana_loop_num(all_seq, ana_pre, ana_end):
                    continue
                all_pair_prob, ana_stable, stem_num = get_pair_prob(
                    all_seq, len(ana_pre), len(all_seq) - len(ana_end)
                )
                ana_stable = statistics.mean(ana_stable) / (100 - stem_num) * 15
                tmp = guess_pair(all_pair_prob)
                if _pair_count(tmp) >= 3:
                    mean_score = statistics.mean(
                        [
                            statistics.mean(tmp),
                            1 - loop_prob,
                            stem_prob,
                            t_pre_at_score,
                            t_behind_at_score,
                            ana_stable,
                        ]
                    )
                    if verbose:
                        print(f"pair-prob:{tmp}")
                        print(f"pair-mean-score:{statistics.mean(tmp)}")
                        print(f"pair-cnt:{_pair_count(tmp)}")
                        print(f"loop_prob:{1 - loop_prob}")
                        print(f"stem_prob:{stem_prob}")
                        print(f"t_pre_at_score:{t_pre_at_score}")
                        print(f"t_behind_at_score:{t_behind_at_score}")
                        print(f"ana_stem_prob:{ana_stable}")
                        print(f"mean_score:{mean_score}")
                        print("-------")
                    res.append(
                        (
                            sequence,
                            all_seq,
                            sequence[indexs[0]:indexs[1] + 1],
                            j - indexs[0],
                            indexs,
                            tmp,
                            statistics.mean(tmp),
                            _pair_count(tmp),
                            1 - loop_prob,
                            stem_prob,
                            t_pre_at_score,
                            t_behind_at_score,
                            ana_stable,
                            mean_score,
                        )
                    )
                    tmp_ans.append(res[-1][-2])
        if tmp_ans:
            best = int(np.argmax(tmp_ans))
            return max(tmp_ans), res[best]
        return 0, []

    ans: List[float] = [0.0] * len(data)
    info: List[list] = [[] for _ in range(len(data))]
    if workers and workers > 1 and len(data) > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(score_one, seq): k for k, seq in enumerate(data)}
            for fut in concurrent.futures.as_completed(futures):
                k = futures[fut]
                score_k, info_k = fut.result()
                ans[k] = score_k
                info[k] = info_k
    else:
        for k, seq in enumerate(data):
            score_k, info_k = score_one(seq)
            ans[k] = score_k
            info[k] = info_k
    return ans, info


__all__ = [
    "find_Hairpin_Loop",
    "get_pair_prob",
    "find_max_list",
    "get_stem_prob_score",
    "guess_loop",
    "guess_pair",
    "guess_ana_loop_num",
    "compute_st_score",
    "find_splice_loop",
]
