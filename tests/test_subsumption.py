"""
R4 and R6 are subsumed (bead Python-7d3): on a table closed under the other rules they force nothing.

R6(a) can only fire on a 2×2 block (a k×k matrix without zero entries has permutation terms of both
signs for k ≥ 3), and there Sylvester's identity is a 3-term Grassmann–Plücker relation; R6(b) is the
R7 circuit of the row of the Schur system, whose x_R coefficients R5 sets first.  So R6 is silent on
tables closed under GP, R5 and R7 alone.  R4's M-matrix and not-M-matrix conclusions are reached by
R2, R5 and R7 (R3 and R8 on the whole matrix): the tests build tables on which R4 applies on a chosen
block and twist and check that the other rules leave it nothing to add.  docs/SIGNRULES.md has the
argument and the census (all patterns at n = 3, 4; random samples at n = 5, 6).  At the level of
propagation R4 is subsumed only up to one branching (test_r4_subsumption.py has the n = 5 table on
which it fires after the other rules, and the exhaustive check from its premises).
"""
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

from signrules import Ecology, Conflict  # noqa: E402
from signrules.engine import Engine  # noqa: E402
from signrules.patterns import random_pattern  # noqa: E402
from signrules.plucker_index import UNKNOWN  # noqa: E402

ALL = Ecology.ALL_RULES


def random_table(eco, rng, frac):
    """The initial table with a random fraction of the open coordinates set at random."""
    chi = eco.initial_table()
    for S in rng.sample(range(len(chi)), rng.randrange(0, max(1, int(frac * len(chi))))):
        if chi[S] == UNKNOWN:
            chi[S] = rng.choice((1, -1))
    return chi


def fires(eco, rule, chi, *args):
    """Does the rule force a sign or raise a conflict on this table?"""
    try:
        return getattr(eco, rule)(chi, [], *args)
    except Conflict:
        return True


def closed(eco, use, chi):
    """The table closed under GP and the rules in ``use``, or None if it is inconsistent."""
    try:
        Engine(eco, use=use).fixpoint(chi, [])
    except Conflict:
        return None
    return chi


@pytest.mark.parametrize("n, trials", [(3, 400), (4, 200)])
def test_r6_is_silent_on_tables_closed_under_gp_r5_r7(n, trials):
    rng = random.Random(100 + n)
    applicable = silent = 0
    for _ in range(trials):
        a, B = random_pattern(n, rng.randrange(10 ** 9))
        eco = Ecology(a, B)
        chi = random_table(eco, rng, 0.4)
        gp = closed(eco, (), list(chi))
        if gp is None:
            continue
        applicable += fires(eco, "rule_schur", gp)
        table = closed(eco, ("R5", "R7"), chi)
        if table is not None:
            assert not fires(eco, "rule_schur", table), (a, B)
            silent += 1
    assert applicable > 10 and silent > trials // 4      # the test is not vacuous


def z_block_table(eco, chi, P, Q, sigma, case, rng):
    """
    Set det B_PP, the entries of the Schur block D = (B/B_PP)_QQ and the right-hand side A_Q so
    that S D P_σ is a Z-matrix with positive diagonal and S A_Q > 0 (case "M") or < 0 ("notM");
    None when that clashes with a sign already in the table.
    """
    dPP_key, entries, rhs = eco._block_data(P, Q)

    def put(S, v):
        if chi[S] == UNKNOWN:
            chi[S] = v
        return chi[S] == v

    if dPP_key is None:
        dPP = 1
    else:
        S, c = dPP_key
        dPP = c * chi[S] if chi[S] != UNKNOWN else rng.choice((1, -1))
        if not put(S, c * dPP):
            return None
    A = {}
    for i, (mi, terms) in zip(Q, rhs):
        if terms is None:
            return None
        t_i = next((par if key is None else par * c * chi[S]
                    for par, S, c, key in terms if key is None or key == (P, P)), None)
        if t_i is None or not all(S is None or put(S, par * c * t_i) for par, S, c, key in terms):
            return None
        A[i] = mi * t_i * dPP
    s = {i: (A[i] if case == "M" else -A[i]) for i in Q}
    for i, j, S, c, key, bit in entries:
        if not put(S, c * s[i] * (1 if sigma[i] == j else -1) * dPP):
            return None
    return s


@pytest.mark.parametrize("n, trials", [(3, 1500), (4, 1500)])
def test_r4_is_silent_on_tables_closed_under_the_other_rules(n, trials):
    rng = random.Random(200 + n)
    others = tuple(r for r in ALL if r != "R4")
    shapes = set()
    for _ in range(trials):
        a, B = random_pattern(n, rng.randrange(10 ** 9))
        eco = Ecology(a, B)
        P = tuple(sorted(rng.sample(range(n), rng.randrange(0, n - 1))))
        rest = [i for i in range(n) if i not in P]
        Q = tuple(sorted(rng.sample(rest, rng.randrange(2, len(rest) + 1))))
        perm = list(Q)
        rng.shuffle(perm)
        sigma = dict(zip(Q, perm))
        case = rng.choice(("M", "notM"))
        chi = eco.initial_table()
        if z_block_table(eco, chi, P, Q, sigma, case, rng) is None:
            continue
        gp = closed(eco, (), list(chi))
        if gp is None or not fires(eco, "z_block", gp, P, Q):
            continue                                    # R4 has nothing to say even before the others
        shapes.add((len(P), len(Q), case))
        table = closed(eco, others, chi)
        if table is not None:
            assert not fires(eco, "z_block", table, P, Q), (a, B, P, Q, sigma, case)
    assert len(shapes) >= 4 and any(p > 0 for p, _, _ in shapes)      # the test is not vacuous


def test_n3_census_whole_trees_r4_r6_idle():
    """All 512 three-species patterns, every completion: R4 and R6 never set a sign nor conflict."""
    import rule_ablation
    recs = [rule_ablation._census_job((3, name, a, B, 10 ** 6, 10 ** 9))
            for name, a, B in rule_ablation.census_patterns(3, 0, 0)]
    assert len(recs) == 512 and not any(r["stopped"] for r in recs)
    assert all(not any(r["fired"].values()) for r in recs)
    status = {s: sum(r["status"] == s for r in recs) for s in ("killed_root", "open")}
    assert status == {"killed_root": 178, "open": 334}
    runs = sum(r["rule_stats"][rule][w]["runs"] for r in recs for rule in ("R4", "R6") for w in ("root", "below"))
    assert runs > 10_000                                 # they did run


def test_idle_auto_names_the_subsumed_rules_and_leaves_the_search_unchanged():
    from signrules import complete
    from signrules.search import subsumed_rules
    assert subsumed_rules(ALL) == ("R4", "R6")
    assert subsumed_rules(("R1", "R2n", "R5", "R6")) == ()            # R7 is absent: R6 is not subsumed there
    assert subsumed_rules(("R1", "R5", "R6", "R7")) == ("R6",)
    assert subsumed_rules(("R1", "R2", "R3", "R4", "R5", "R7", "R8")) == ("R4",)
    with pytest.raises(ValueError):
        complete(Ecology(*random_pattern(3, 1)), idle="all")
    kw = dict(keep_derivations=False, max_completions=5, max_nodes=300, stats=True)
    ran = 0
    for seed in (176, 3, 11):
        eco = Ecology(*random_pattern(4, seed))
        auto, full = complete(eco, **kw), complete(eco, idle=(), **kw)
        assert (auto.nodes, auto.conflicts, auto.completions, auto.root_log, auto.root_conflict) == \
               (full.nodes, full.conflicts, full.completions, full.root_log, full.root_conflict)
        assert all(auto.rule_stats[r][w]["runs"] == 0 for r in ("R4", "R6") for w in ("root", "below"))
        ran += sum(full.rule_stats[r][w]["runs"] for r in ("R4", "R6") for w in ("root", "below"))
    assert ran > 0                                       # some pattern survives the cheap rules and reaches them
