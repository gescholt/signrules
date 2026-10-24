"""
R4 is subsumed one branching step deep, exhaustively at n = 3, 4 (bead Python-jr9).

The propagators are monotone, so it suffices to start from the minimal table on which R4 applies (its
premises) and check the closure under the other rules against R4's conclusions: every principal and
almost-principal minor of the twisted block in the M-matrix case, and in the not-M case, for every
ordering of the block and every position, the negative leading minor once the others are positive.
scripts/r4_subsumption.py enumerates the configurations (block, twist, case, sign of det B_PP, growth
signs); up to relabelling P and Q are initial segments.  Under the 3-term relations with R1 and R7
alone every configuration passes at the fixpoint except four whole-matrix ones at n = 4 (P empty,
Q = [n]), which pass "by probe": the sign R4 forces is not derived, but its opposite conflicts at the
fixpoint, so one branching on that minor settles it.

R4 is not subsumed at the fixpoint level: at n = 5 (not-M case, P empty, |Q| = 4) there are tables,
consistent under the other rules, on which R4 forces a 3 x 3 leading minor of the twisted block
negative when the 1 x 1, 2 x 2 and 4 x 4 ones are positive; the opposite sign is refuted one step
deeper (an R2 Cramer conflict).  The census of whole search trees (test_subsumption.py) shows the DFS,
which decides small minors first, never reaches such a table.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

import r4_subsumption as R  # noqa: E402
from signrules import Ecology, Conflict  # noqa: E402
from signrules.engine import Engine  # noqa: E402
from signrules.plucker_index import UNKNOWN  # noqa: E402

OTHERS = tuple(r for r in Ecology.ALL_RULES if r != "R4")


@pytest.mark.parametrize("n", [3, 4])
def test_r4_conclusions_follow_from_its_premises(n):
    use = ("R1", "R7")
    counts, probed = {}, []
    for cfg in R.configs(n, wlog=True):
        r = R.check_config(cfg + (use,))
        counts[r["status"]] = counts.get(r["status"], 0) + 1
        if r["status"] == "ok by probe":
            n_, P, Q, perm, case, dPP, a = cfg
            probed.append(dict(P=P, Q=Q, sigma=perm, case=case, dPP=dPP, a=a))
    assert counts["ok"] > 30 and "R4 forces more" not in counts and "R4 conflicts on the premises" not in counts
    for f in probed:                                     # whole-matrix cases, not-M: all leading minors
        assert f["P"] == () and f["Q"] == tuple(range(n)) and f["case"] == "notM", f
    assert len(probed) == (0 if n == 3 else 4)          # n = 4: a Z-matrix with a < 0; three involutions with a > 0


def test_r4_forces_a_sign_the_others_only_refute_by_probing():
    """The n = 5 table: R4 fires at the fixpoint of the other rules; the opposite sign conflicts one step deeper."""
    a = [-1, 1, 1, 1, 1]
    B = np.array([[1, -1, -1, -1, 1], [1, 1, -1, 1, -1], [1, 1, 1, -1, -1], [1, -1, 1, 1, -1], [1, 1, 1, 1, 1]])
    eco = Ecology(a, B)
    P, Q, perm = (), (0, 1, 2, 3), (0, 2, 3, 1)
    sigma = dict(zip(Q, perm))
    chi0, s, A = R.premises(eco, P, Q, sigma, "notM", 1)
    assert all(s[i] * A[i] < 0 for i in Q)              # the not-M case of the twist
    leading = [eco.block_minor(None, P, s, sigma, tuple(range(m)), tuple(range(m))) for m in range(1, 5)]
    chi = eco.initial_table()
    for S, v in enumerate(chi0):
        if v != UNKNOWN:
            assert chi[S] in (UNKNOWN, v)
            chi[S] = v
    for m in (1, 2, 4):                                 # leading minors 1, 2, 4 of Z positive
        coef, key = leading[m - 1]
        S, c = eco._minor_key[key]
        chi[S] = c * coef
    Engine(eco, use=OTHERS).fixpoint(chi, [])
    coef3, key3 = leading[2]
    S3, c3 = eco._minor_key[key3]
    assert chi[S3] == UNKNOWN                            # the others leave the 3 x 3 leading minor open
    with_r4 = list(chi)
    assert eco.z_block(with_r4, [], P, Q)                # R4 fires and sets it negative
    assert c3 * coef3 * with_r4[S3] == -1
    assert R.refuted(eco, OTHERS, Engine(eco, use=OTHERS), chi, S3, c3 * coef3)
    assert R.verdict(eco, chi0, P, Q, sigma, "notM", 1, s, OTHERS) == "ok by probe"
    wrong = list(chi)
    wrong[S3] = c3 * coef3
    with pytest.raises(Conflict):
        Engine(eco, use=OTHERS).fixpoint(wrong, [])
