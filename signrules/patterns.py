"""
Sign patterns used across the package.

* ``ELEVEN``: the eleven four-species ecologies conjectured impossible in Meng et al. (2025) and proved
  impossible in Çelik–Isoldi–Portakal–Zucal (2026); names S-a … S-e are Figure 2 (a)–(e) of that paper,
  A-a … A-f Figure 3 (a)–(f).
* ``parity_pattern``: the benefactor-permutation family (rule R8; every species of a cycle is obligate and
  benefits only from the next one, the fixed species are free and benefit from everyone).
* ``TREES``: presets whose search tree is worth watching, with the rules to run them under (``TREE_RULES``).
* ``PATTERNS``: five- and six-species presets of the live search (``signrules run --preset``);
  ``presets()`` is both tables, ``resolve_preset`` finds one by its short id (``S-d``) or a prefix.
* ``all_patterns``, ``canonical_classes``, ``trivial``: enumeration of the patterns on n species and of
  their classes up to relabelling; ``encode`` / ``decode`` / ``canonical`` / ``key``: the integer code of a
  pattern used by the extension screen (``signrules.screen``) and in the data files.
* ``random_pattern``: a uniformly random pattern.

Conventions: species are 0-based, b_ii > 0, B[i, j] < 0 means species j benefits species i (so predation
"i eats j" is B[i, j] < 0, B[j, i] > 0); a_i < 0 means species i cannot grow alone (it is obligate).
"""
import itertools
import random as _random
from functools import lru_cache
from itertools import permutations
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

Pattern = Tuple[List[int], np.ndarray]


# ----------------------------------------------------------------------------- four species
def sym(pairs_neg, n=4):
    """Symmetric pattern: listed unordered pairs mutualistic (-), all others competitive (+)."""
    B = np.ones((n, n), dtype=int)
    for i, j in pairs_neg:
        B[i, j] = B[j, i] = -1
    return B


def cycle_pattern(n, cycle, extra=None):
    """All pairs compete except the directed cycle i -> i+ (i eats i+), plus `extra` edits."""
    B = np.ones((n, n), dtype=int)
    for i, ip in zip(cycle, cycle[1:] + cycle[:1]):
        B[i, ip], B[ip, i] = -1, 1
    for (i, j), v in (extra or {}).items():
        B[i, j] = v
    return B


ELEVEN: Dict[str, Tuple[List[int], object]] = {
    # ---- symmetric (Figure 2) ----
    "S-a obligate mutualism": ([-1, -1, -1, -1], sym([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])),
    "S-b {0}|{1,2,3}, triple obligate": ([1, -1, -1, -1], sym([(1, 2), (1, 3), (2, 3)])),
    "S-c {0,1}|{2,3}, {2,3} obligate, a01=(+,+)": ([1, 1, -1, -1], sym([(0, 1), (2, 3)])),
    "S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)": ([1, -1, -1, -1], sym([(0, 1), (2, 3)])),
    "S-e {0,1}|{2,3}, both obligate": ([-1, -1, -1, -1], sym([(0, 1), (2, 3)])),
    # ---- asymmetric (Figure 3) ----
    "A-a Lemma 3.8": ([1, 1, -1, -1], [[1, 1, 1, 1], [-1, 1, -1, -1], [1, 1, 1, -1], [1, 1, -1, 1]]),
    "A-b Lemma 3.9": ([1, 1, -1, -1], [[1, -1, -1, -1], [-1, 1, -1, -1], [1, 1, 1, -1], [1, 1, -1, 1]]),
    "A-c Lemma 3.10": ([1, -1, -1, -1], [[1, 1, 1, 1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1]]),
    "A-d obligate 4-cycle (Thm 3.4)": ([-1, -1, -1, -1], cycle_pattern(4, [0, 1, 2, 3])),
    "A-e competitor + 3-cycle (Thm 3.5)": ([1, -1, -1, -1], cycle_pattern(4, [1, 2, 3])),
    "A-f super-predator + 3-cycle (Thm 3.6)": ([1, -1, -1, -1],
                                               cycle_pattern(4, [1, 2, 3], {(0, 1): -1, (0, 2): -1, (0, 3): -1})),
}


# ----------------------------------------------------------------------------- the parity family (R8)
def parity_pattern(n, cycles, fixed):
    """Sign pattern of the parity family; cycles are lists of species, i benefits from the next one."""
    assert sorted(sum(cycles, []) + list(fixed)) == list(range(n))
    B = np.ones((n, n), dtype=int)
    a = np.ones(n, dtype=int)
    for i in fixed:
        B[i, [j for j in range(n) if j != i]] = -1
    for cyc in cycles:
        for i, nxt in zip(cyc, cyc[1:] + cyc[:1]):
            B[i, nxt] = -1
            a[i] = -1
    return list(a), B


def cycle_types(n):
    """All (cycles, fixed) shapes up to relabeling: partitions of a subset of [n] into parts ≥ 2."""
    def parts(m, max_part):
        if m == 0:
            yield []
            return
        for p in range(min(m, max_part), 1, -1):
            for rest in parts(m - p, p):
                yield [p] + rest
    for size in range(2, n + 1):
        for shape in parts(size, size):
            cycles, start = [], 0
            for length in shape:
                cycles.append(list(range(start, start + length)))
                start += length
            yield cycles, list(range(size, n))


# ----------------------------------------------------------------------------- presets (n = 5, 6)
def cycle(n, cyc, a_sign=-1):
    """All pairs compete except i -> next (i eats next); species on the cycle get growth-rate sign a_sign."""
    B = np.ones((n, n), dtype=int)
    a = np.ones(n, dtype=int)
    for i, nxt in zip(cyc, cyc[1:] + cyc[:1]):
        B[i, nxt], B[nxt, i] = -1, 1
        a[i] = a_sign
    return list(a), B


def _with(pattern, a_edits=None, B_edits=None):
    a, B = pattern
    a, B = list(a), np.array(B, dtype=int)
    for i, v in (a_edits or {}).items():
        a[i] = v
    for (i, j), v in (B_edits or {}).items():
        B[i, j] = v
    return a, B


def _n6(a, rows):
    return list(a), np.array(rows, dtype=int)


def _pairs_plus_benefactor(a):
    """Two mutualistic pairs {0,1}, {2,3} and species 4 benefiting from everybody (R8, one wrong growth sign)."""
    B = np.ones((5, 5), dtype=int)
    for i, j in ((0, 1), (1, 0), (2, 3), (3, 2)):
        B[i, j] = -1
    B[4, :4] = -1
    return list(a), B


# Presets of the live search. "parity_k2" and "r8_*" are killed by R8 alone (the real eigenvalue of
# benefactor-permutation patterns) and open under R1-R7; the next three are five-species patterns next to
# proved-impossible ones ("kill #k" is the id in signrules/data/n5/killed_patterns.json); "n6_hard_*" are
# six-species patterns that are possible, whose certificates of possibility (tropical witnesses, in the
# repository feasible_ecologies) were the hardest to find.
PATTERNS: Dict[str, Pattern] = {
    # even-k parity family with a fixed point: two obligate pairs plus one universal beneficiary (killed by R8 only)
    "parity_k2": parity_pattern(5, [[0, 1], [2, 3]], [4]),
    # the obligate 5-cycle (kill #21) with one species allowed to grow on its own
    "5cycle_one_growing": _with(cycle(5, [0, 1, 2, 3, 4]), a_edits={0: 1}),
    # competitor + obligate 4-cycle (kill #20) where the competitor is also eaten by species 1
    "competitor_4cycle_eaten": _with(cycle(5, [1, 2, 3, 4]), B_edits={(1, 0): -1}),
    # obligate 3-cycle + pair (kill #17) with one cross-benefit between the two groups
    "3cycle_pair_bridge": _with(cycle(5, [0, 1, 2]), B_edits={(3, 4): -1, (4, 3): -1, (0, 3): -1}),
    # killed by R8 alone (n = 5)
    "r8_pairs_benefactor_+---+": _pairs_plus_benefactor([1, -1, -1, -1, 1]),
    "r8_pairs_benefactor_-----": _pairs_plus_benefactor([-1, -1, -1, -1, -1]),
    # killed by R8 alone (n = 6; codes of signrules/data/n6/extension_screen.jsonl)
    "r8_n6_58455985": _n6([1, -1, -1, -1, 1, 1], [[1, -1, 1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, 1, 1, -1, 1, 1],
                                                 [1, 1, -1, 1, 1, 1], [-1, -1, -1, -1, 1, -1], [-1, -1, -1, -1, -1, 1]]),
    "r8_n6_58455952": _n6([-1, -1, -1, -1, 1, -1], [[1, -1, 1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, 1, 1, -1, 1, 1],
                                                   [1, 1, -1, 1, 1, 1], [-1, -1, -1, -1, 1, -1], [-1, -1, -1, -1, -1, 1]]),
    "r8_n6_1593702248": _n6([-1, -1, -1, 1, -1, 1], [[1, 1, -1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, -1, 1, 1, 1, 1],
                                                    [1, 1, 1, 1, -1, 1], [1, 1, 1, -1, 1, 1], [-1, -1, -1, -1, -1, 1]]),
    "r8_n6_1593702208": _n6([-1, -1, -1, -1, -1, -1], [[1, 1, -1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, -1, 1, 1, 1, 1],
                                                      [1, 1, 1, 1, -1, 1], [1, 1, 1, -1, 1, 1], [-1, -1, -1, -1, -1, 1]]),
    "r8_n6_1593702241": _n6([1, -1, -1, -1, -1, 1], [[1, 1, -1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, -1, 1, 1, 1, 1],
                                                    [1, 1, 1, 1, -1, 1], [1, 1, 1, -1, 1, 1], [-1, -1, -1, -1, -1, 1]]),
    # possible six-species patterns whose tropical witnesses (feasible_ecologies) were the hardest to find
    "n6_hard_1065088896": _n6([-1, -1, -1, -1, -1, -1], [[1, -1, 1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, 1, 1, -1, 1, 1],
                                                        [1, 1, -1, 1, 1, 1], [1, 1, 1, 1, 1, -1], [-1, -1, -1, -1, -1, 1]]),
    "n6_hard_1316222848": _n6([-1, -1, -1, -1, -1, -1], [[1, -1, 1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, 1, 1, -1, -1, 1],
                                                        [1, 1, -1, 1, -1, 1], [1, 1, -1, -1, 1, 1], [-1, -1, -1, -1, -1, 1]]),
    "n6_hard_3741182912": _n6([-1, -1, -1, -1, -1, -1], [[1, 1, 1, 1, 1, -1], [1, 1, -1, 1, 1, 1], [1, -1, 1, 1, 1, 1],
                                                        [1, 1, 1, 1, -1, 1], [1, 1, 1, -1, 1, 1], [1, -1, -1, -1, -1, 1]]),
    "n6_hard_15518791488": _n6([-1, -1, -1, -1, -1, -1], [[1, 1, -1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, -1, 1, 1, 1, 1],
                                                         [1, 1, 1, 1, -1, -1], [1, 1, 1, -1, 1, -1], [1, 1, 1, -1, -1, 1]]),
    "n6_hard_16072439616": _n6([-1, -1, -1, -1, -1, -1], [[1, 1, -1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, -1, 1, 1, 1, 1],
                                                         [1, 1, 1, 1, 1, -1], [1, 1, 1, 1, 1, -1], [1, 1, 1, -1, -1, 1]]),
    "n6_hard_33277474544": _n6([-1, -1, -1, -1, 1, 1], [[1, 1, 1, -1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, -1, 1, 1, 1, 1],
                                                       [1, 1, -1, 1, 1, 1], [1, 1, 1, 1, 1, -1], [1, 1, 1, 1, -1, 1]]),
    "n6_hard_1132197808": _n6([-1, -1, -1, -1, 1, 1], [[1, -1, 1, 1, 1, 1], [-1, 1, 1, 1, 1, 1], [1, 1, 1, -1, 1, 1],
                                                      [1, 1, -1, 1, 1, 1], [-1, -1, -1, -1, 1, 1], [-1, -1, -1, -1, -1, 1]]),
}


# Presets whose search tree is worth watching: the root propagation does not decide them under the rules given
# with each (TREE_RULES; all rules when absent), and the search ends in a tree of completions and contradictions.  Under all rules every
# pattern on n <= 4 species that survives the root ends almost only in completions (at most 7 contradictions in
# 1,734 nodes, and no search kills); the contradictions show best under the Annals algorithm R1, R2n.
TREES: Dict[str, Pattern] = {
    "tree_S-d_annals": ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"],
    "tree_S-a_search_kill": ELEVEN["S-a obligate mutualism"],
    "tree_hub3": ([1, 1, 1], np.array([[1, 1, -1], [1, 1, -1], [-1, -1, 1]])),
    "tree_chain4": ([1, -1, -1, -1], np.array([[1, 1, -1, -1], [-1, 1, 1, 1], [1, 1, 1, -1], [1, 1, -1, 1]])),
    "tree_preyed_pair4": ([1, 1, -1, -1], np.array([[1, 1, -1, -1], [-1, 1, -1, -1], [1, 1, 1, -1], [1, 1, -1, 1]])),
    "tree_mutualists4": ([1, 1, 1, -1], np.array([[1, 1, 1, -1], [1, 1, -1, 1], [-1, -1, 1, 1], [-1, 1, 1, 1]])),
}
TREE_RULES: Dict[str, Tuple[str, ...]] = {          # the others run under all rules, the default
    "tree_S-d_annals": ("R1", "R2n"), "tree_S-a_search_kill": ("R1", "R2n"),
}
TREE_NOTES: Dict[str, str] = {
    "tree_S-d_annals": "S-d under the Annals algorithm: 64 completions, 131 nodes, 4 contradictions",
    "tree_S-a_search_kill": "S-a under R1, R2n: survives the root, every branch closes (5 contradictions)",
    "tree_hub3": "competitors 0, 1, each a mutualist of 2: 30 completions, 60 nodes, 1 contradiction",
    "tree_chain4": "obligate pair {2,3}, free 0 eats both, obligate 1 eats 0: 416 completions, 2 contradictions",
    "tree_preyed_pair4": "obligate pair {2,3} eaten by the free pair {0,1}: 416 completions, 1 contradiction",
    "tree_mutualists4": "free pair {1,2}, obligate 3 with 0, 2 eats 0: 864 completions, 7 contradictions",
}


def preset_rules(name: str) -> Optional[Tuple[str, ...]]:
    """The rules a preset is meant to be run with (the ``TREES`` presets), or None."""
    return TREE_RULES.get(name)


_SIGN = {"+": 1, "1": 1, "+1": 1, "-": -1, "−": -1, "-1": -1}


def parse_signs(text: str) -> List[int]:
    """``"1,-1,-1"``, ``"+,-,-"`` or the compact ``"+--"`` (``−`` is read as ``-``) -> [1, -1, -1]."""
    text = text.strip()
    if "," in text or text in _SIGN:
        parts = [t.strip() for t in text.split(",")]
    else:
        parts = list(text)
    try:
        return [_SIGN[t] for t in parts]
    except KeyError as e:
        raise ValueError(f"{text!r}: signs are written 1,-1,... or +,-,... or +-...; {e.args[0]!r} is none of these") from None


def parse_pattern(a_str: str, B_str: str) -> Pattern:
    """``"-1,-1,1"`` and ``"1,-1,1/1,1,-1/-1,1,1"`` (rows separated by ``/``; every row may also be written
    ``"+-+"``) -> (a, B).  Raises ValueError when B is not n×n for n = len(a) or a diagonal entry is not +."""
    a = parse_signs(a_str)
    rows = [parse_signs(row) for row in B_str.split("/")] if B_str.strip() else []
    n = len(a)
    if len(rows) != n or any(len(r) != n for r in rows):
        raise ValueError(f"a has {n} signs, so B needs {n} rows of {n} signs; "
                         f"got {len(rows)} rows of lengths {[len(r) for r in rows]}")
    B = np.array(rows, dtype=int)
    bad = [i for i in range(n) if B[i, i] <= 0]
    if bad:
        raise ValueError(f"the diagonal of B is positive (self-regulation); entries {', '.join(f'b_{i}{i}' for i in bad)} are not")
    return a, B


def short_id(name: str) -> str:
    """The short name of a preset: ``"S-d"`` for ``"S-d {0,1}|{2,3}, ..."``, the name itself for the others."""
    return name.split()[0]


def presets() -> Dict[str, Pattern]:
    """Every named pattern: the search-tree presets (``TREES``), the eleven four-species patterns (``ELEVEN``),
    then the five- and six-species presets (``PATTERNS``), as plain lists."""
    out = {}
    for table in (TREES, ELEVEN, PATTERNS):
        for k, (a, B) in table.items():
            out[k] = (list(map(int, a)), [list(map(int, r)) for r in B])
    return out


def resolve_preset(name: str) -> str:
    """The full preset name for ``name``: the full name, a short id (``S-d``, ``A-e``) or an unambiguous
    prefix of either, ignoring case.  Raises ValueError, naming the candidates, otherwise."""
    names = list(presets())
    low = name.lower()
    for test in (lambda k: k == name, lambda k: short_id(k).lower() == low, lambda k: k.lower().startswith(low)):
        hits = [k for k in names if test(k)]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise ValueError(f"preset {name!r} is ambiguous: " + ", ".join(short_id(k) for k in hits))
    import difflib
    close = difflib.get_close_matches(name, [short_id(k) for k in names], n=5, cutoff=0.4)
    hint = f"; did you mean {', '.join(close)}?" if close else ""
    raise ValueError(f"no preset {name!r}{hint} ('signrules presets' lists them)")


# ----------------------------------------------------------------------------- random and exhaustive
def random_pattern(n: int, seed: Optional[int] = None):
    """A uniformly random sign pattern on n species: signs of a and of the off-diagonal of B; b_ii > 0."""
    rng = _random.Random(seed)
    a = [rng.choice((-1, 1)) for _ in range(n)]
    B = [[1 if i == j else rng.choice((-1, 1)) for j in range(n)] for i in range(n)]
    return a, B


def all_patterns(n) -> Iterator[Pattern]:
    """Every labelled pattern on n species: 2^n growth signs times 2^(n(n-1)) off-diagonal signs."""
    for a in itertools.product([1, -1], repeat=n):
        for off_signs in itertools.product([1, -1], repeat=n * (n - 1)):
            B = np.ones((n, n), dtype=int)
            k = 0
            for i in range(n):
                for j in range(n):
                    if i != j:
                        B[i, j] = off_signs[k]
                        k += 1
            yield list(a), B


def trivial(a, B) -> bool:
    """A species with a < 0 and no benefactor cannot have (Bx)_i = a_i < 0 with x > 0."""
    n = len(a)
    return any(a[i] < 0 and all(B[i][j] > 0 for j in range(n) if j != i) for i in range(n))


def canonical_classes(n) -> List[Pattern]:
    """One representative of every class of labelled patterns up to relabelling (brute force over S_n)."""
    ps = perms(n)
    seen, reps = set(), []
    for a, B in all_patterns(n):
        k = min((tuple(a[i] for i in p), tuple(int(B[i, j]) for i in p for j in p)) for p in ps)
        if k not in seen:
            seen.add(k)
            reps.append((a, B))
    return reps


# ----------------------------------------------------------------------------- integer codes
@lru_cache(maxsize=None)
def off(n):
    """The off-diagonal slots (i, j), i != j, in row-major order."""
    return [(i, j) for i in range(n) for j in range(n) if i != j]


@lru_cache(maxsize=None)
def perms(n):
    """All permutations of range(n), as tuples."""
    return list(permutations(range(n)))


def key(a, B):
    """Canonical form of a pattern under relabelling (brute force over S_n)."""
    n = len(a)
    return min((tuple(int(a[i]) for i in p), tuple(int(B[i][j]) for i in p for j in p)) for p in perms(n))


def encode(a, B) -> int:
    """Pattern -> int: n bits of a then n(n-1) bits of off-diagonal B (bit = 1 for +)."""
    n = len(a)
    code = 0
    for i in range(n):
        code |= int(a[i] > 0) << i
    for k, (i, j) in enumerate(off(n)):
        code |= int(B[i][j] > 0) << (n + k)
    return code


def decode(code: int, n: int) -> Pattern:
    """Inverse of ``encode``."""
    a = [1 if (code >> i) & 1 else -1 for i in range(n)]
    B = np.eye(n, dtype=int)
    for k, (i, j) in enumerate(off(n)):
        B[i, j] = 1 if (code >> (n + k)) & 1 else -1
    return a, B


def canonical(codes: Sequence[int], n: int) -> np.ndarray:
    """Minimal code over the n! relabellings, vectorised over an array of codes."""
    codes = np.asarray(codes, dtype=np.int64)
    m = n + len(off(n))
    bits = ((codes[:, None] >> np.arange(m)) & 1).astype(np.int64)
    slot = {ij: k for k, ij in enumerate(off(n))}
    weights = np.int64(1) << np.arange(m, dtype=np.int64)
    best = codes.copy()
    for p in perms(n):
        # relabel i -> p[i]: new entry (p[i], p[j]) holds old (i, j); new a_{p[i]} = old a_i
        inv = np.argsort(p)
        idx = [inv[i] for i in range(n)] + [n + slot[(inv[i], inv[j])] for (i, j) in off(n)]
        np.minimum(best, bits[:, idx] @ weights, out=best)
    return best


__all__ = ["Pattern", "ELEVEN", "PATTERNS", "sym", "cycle_pattern", "parity_pattern", "cycle_types", "cycle",
           "TREES", "TREE_RULES", "TREE_NOTES", "preset_rules",
           "parse_signs", "parse_pattern", "short_id", "presets", "resolve_preset", "random_pattern", "all_patterns", "trivial", "canonical_classes",
           "off", "perms", "key", "encode", "decode", "canonical"]
