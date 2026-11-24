"""
R8 (real eigenvalue): the benefactor-permutation patterns (scripts/families.py, cases (c) and (d)).
The statement is checked numerically, and R8 is checked never to kill an n = 3 pattern that
R1-R7 leave open (all of those are certified possible: signrules.data.census(4)).
"""
import itertools

import numpy as np
import pytest

from signrules.patterns import parity_pattern

from signrules import Ecology, complete
from signrules import labels as L

PARITY_K2 = parity_pattern(5, [[0, 1], [2, 3]], [4])


@pytest.mark.parametrize("engine", ["trigger", "sweep"])
def test_parity_k2_is_killed_by_r8(engine):
    eco = Ecology(*PARITY_K2)
    assert complete(eco, use=[r for r in Ecology.ALL_RULES if r != "R8"], engine=engine, max_nodes=0).root_conflict is None
    res = complete(eco, engine=engine, max_completions=1)
    assert res.root_conflict.startswith("R8 real eigenvalue")
    assert L.rule(res.root_conflict, 5).rule == "R8"
    assert L.conflict(eco, res.root_conflict)[1].startswith("R8")


def _sample(a, B, rng):
    """C = B diag(x*) with the signs of B and row sums with the signs of a (log-uniform magnitudes)."""
    C = B * np.exp(rng.uniform(-3, 3, size=B.shape))
    for i in range(len(a)):
        pos, neg = C[i][C[i] > 0].sum(), -C[i][C[i] < 0].sum()
        f = rng.uniform(1.01, 5)
        if a[i] < 0 and pos >= neg:
            C[i][C[i] < 0] *= f * pos / neg
        elif a[i] > 0 and neg >= pos:
            C[i][C[i] < 0] *= pos / (f * neg)
    assert (np.sign(C.sum(1)) == a).all()
    return C


@pytest.mark.parametrize("n,cycles,fixed", [(5, [[0, 1], [2, 3]], [4]), (6, [[0, 1], [2, 3]], [4, 5]),
                                            (4, [[0, 1]], [2, 3]), (6, [[0, 1, 2], [3, 4]], [5])])
def test_parity_patterns_have_a_negative_real_eigenvalue(n, cycles, fixed):
    a, B = parity_pattern(n, cycles, fixed)
    a, B = np.array(a), np.array(B)
    assert complete(Ecology(a, B), use=("R8",), max_nodes=0).root_conflict is not None
    rng = np.random.default_rng(n)
    for _ in range(300):
        ev = np.linalg.eigvals(_sample(a, B, rng))           # spectrum of B X; the Jacobian is similar to -B X
        assert any(abs(e.imag) < 1e-9 * abs(e) and e.real < 0 for e in ev)


def test_r8_kills_nothing_that_r1_to_r7_leave_open_at_n3():
    rest = [r for r in Ecology.ALL_RULES if r != "R8"]
    fired = 0
    for a in itertools.product((1, -1), repeat=3):
        for off in itertools.product((1, -1), repeat=6):
            B = [[1, off[0], off[1]], [off[2], 1, off[3]], [off[4], off[5], 1]]
            eco = Ecology(list(a), B)
            if complete(eco, use=("R8",), max_nodes=0).root_conflict:
                fired += 1
                assert not complete(eco, use=rest, max_completions=1, keep_derivations=False).completions
    assert fired > 0


# ---------------------------------------------------------------- case (d): one wrong growth sign

def _flip(pattern, i):
    a, B = pattern
    a = list(a)
    a[i] = -a[i]
    return a, B


# the two five-species and the five six-species patterns that R1-R7 and the witnesses leave open
ONE_WRONG = [_flip(parity_pattern(5, [[0, 1], [2, 3]], [4]), 0),           # (+,-,-,-,+)
             _flip(parity_pattern(5, [[0, 1], [2, 3]], [4]), 4),           # (-,-,-,-,-)
             _flip(parity_pattern(6, [[0, 1], [2, 3]], [4, 5]), 0),        # (+,-,-,-,+,+)
             _flip(parity_pattern(6, [[0, 1], [2, 3]], [4, 5]), 5),        # (-,-,-,-,+,-)
             _flip(parity_pattern(6, [[0, 1, 2], [3, 4]], [5]), 5),        # (-,-,-,-,-,-)
             _flip(parity_pattern(6, [[0, 1, 2], [3, 4]], [5]), 0),        # (+,-,-,-,-,+)
             _flip(parity_pattern(6, [[0, 1, 2], [3, 4]], [5]), 3)]        # (-,-,-,+,-,+)


@pytest.mark.parametrize("engine", ["trigger", "sweep"])
@pytest.mark.parametrize("pattern", ONE_WRONG, ids=lambda p: "".join("+" if x > 0 else "-" for x in p[0]))
def test_one_wrong_sign_even_cycles_is_killed_by_r8(pattern, engine):
    eco = Ecology(*pattern)
    res = complete(eco, engine=engine, max_completions=1)
    assert res.root_conflict is not None and res.root_conflict.startswith("R8 real eigenvalue")
    assert "even number of cycles" in res.root_conflict
    assert L.rule(res.root_conflict, len(pattern[0])).rule == "R8"


@pytest.mark.parametrize("pattern", ONE_WRONG, ids=lambda p: "".join("+" if x > 0 else "-" for x in p[0]))
def test_one_wrong_sign_patterns_have_a_negative_real_eigenvalue(pattern):
    a, B = np.array(pattern[0]), np.array(pattern[1])
    rng = np.random.default_rng(len(a))
    for _ in range(300):
        C = _sample(a, B, rng)
        ev = np.linalg.eigvals(C)
        assert (ev.real < 0).any()                                   # never positive stable
        if np.linalg.det(C) > 0:                                     # case (d): a negative real eigenvalue
            assert any(abs(e.imag) < 1e-9 * abs(e) and e.real < 0 for e in ev)


def test_one_wrong_sign_needs_even_cycles():
    """One wrong growth sign with k odd (obligate 3-cycle, one fixed point): R8 stays silent, and the
    pattern (all four species obligate, a 3-cycle under a species benefiting from everyone) is possible: it has a
    certificate of possibility in signrules.data.census(4)."""
    a, B = _flip(parity_pattern(4, [[0, 1, 2]], [3]), 3)              # (-,-,-,-), k = 1
    assert complete(Ecology(a, B), use=("R8",), max_nodes=0).root_conflict is None
    a, B = _flip(parity_pattern(5, [[0, 1], [2, 3]], [4]), 0)
    a[4] = -a[4]                                                       # two wrong signs, k = 2
    assert complete(Ecology(a, B), use=("R8",), max_nodes=0).root_conflict is None


def _twisted(a, B):
    """Cheap prefilter for the hypothesis of R8: a sign vector s and a permutation with s_i B_ij > 0
    exactly at j = sigma(i), and at most one i with s_i a_i < 0."""
    n = len(a)
    for s in itertools.product((1, -1), repeat=n):
        sigma = []
        for i in range(n):
            pos = [j for j in range(n) if s[i] * B[i][j] > 0]
            if len(pos) != 1:
                break
            sigma.append(pos[0])
        if len(sigma) == n and len(set(sigma)) == n and sum(s[i] * a[i] < 0 for i in range(n)) <= 1:
            return True
    return False


def test_r8_with_one_wrong_sign_kills_nothing_open_at_n4():
    """Every n = 4 class that R1-R7 leave open has a certificate of possibility (signrules.data.census(4)), so R8 must not fire on any of
    them (exhaustive over the 65,536 four-species patterns)."""
    rest = [r for r in Ecology.ALL_RULES if r != "R8"]
    fired = checked = 0
    for a in itertools.product((1, -1), repeat=4):
        for off in itertools.product((1, -1), repeat=12):
            B = [[1, off[0], off[1], off[2]], [off[3], 1, off[4], off[5]],
                 [off[6], off[7], 1, off[8]], [off[9], off[10], off[11], 1]]
            if not _twisted(a, B):
                continue
            eco = Ecology(list(a), B)
            res = complete(eco, use=("R8",), max_nodes=0)
            if res.root_conflict:
                fired += 1
                if "even number of cycles" in res.root_conflict:
                    checked += 1
                assert not complete(eco, use=rest, max_completions=1, keep_derivations=False).completions
    assert fired == 35 and checked > 0
