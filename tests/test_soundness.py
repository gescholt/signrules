"""
Soundness: the chirotope of a genuine feasible-stable Lotka-Volterra system must
never be rejected by the sign rules, and the index conventions must be right.
"""
import numpy as np
import pytest

from signrules import PluckerIndex, Ecology, Conflict, build_relation_ids, propagate, violates


def random_feasible_stable(n, rng, balanced=None, tries=2000):
    """Random (a, B) with b_ii > 0, x* > 0 and Hurwitz Jacobian; optional signature twist."""
    for _ in range(tries):
        x = rng.uniform(0.2, 2.0, n)
        B = rng.normal(size=(n, n))
        np.fill_diagonal(B, rng.uniform(0.5, 2.0, n))
        if balanced is not None:
            s = np.array(balanced)
            off = -np.abs(B) * np.outer(s, s)
            np.fill_diagonal(off, np.diag(B))
            B = off
        a = B @ x
        J = -np.diag(x) @ B
        if np.all(np.linalg.eigvals(J).real < 0) and np.all(np.abs(a) > 1e-3):
            return a, B
    raise RuntimeError("no instance")


@pytest.mark.parametrize("n", [2, 3, 4])
def test_index_conventions(n):
    PluckerIndex(n).check_conventions(seed=n)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_initial_table_matches_numeric(n):
    rng = np.random.default_rng(1)
    a, B = random_feasible_stable(n, rng)
    eco = Ecology(np.sign(a).astype(int), np.sign(B).astype(int))
    chi = eco.initial_table()
    num = eco.idx.numeric_chi(a, B)
    assert all(c in (0, v) for c, v in zip(chi, num))


@pytest.mark.parametrize("n,balanced", [(3, None), (4, None), (4, [1, -1, 1, -1]), (4, [1, 1, -1, -1]), (3, [1, 1, 1])])
def test_rules_accept_true_chirotope(n, balanced):
    rng = np.random.default_rng(7)
    rels = build_relation_ids(n)
    for trial in range(15):
        a, B = random_feasible_stable(n, rng, balanced)
        eco = Ecology(np.sign(a).astype(int), np.sign(B).astype(int))
        num = eco.idx.numeric_chi(a, B)
        if any(v == 0 for v in num):
            continue
        assert not violates(num, rels)
        chi = list(num)
        try:
            eco.apply_rules(chi, [])           # must not raise on a realized table
        except Conflict as c:
            pytest.fail(f"true system rejected: {c}")
        # partial table: rules must never force a sign different from the truth
        part = eco.initial_table()
        ok, _ = propagate(part, rels)
        assert ok
        eco.apply_rules(part, [])
        assert all(c in (0, v) for c, v in zip(part, num))


@pytest.mark.parametrize("n", [3, 4])
def test_tope_rule_subsumes_cramer(n):
    """R7 (circuit orthogonality) determines at least every sign that R2 (subsystem Cramer) does."""
    rng = np.random.default_rng(11)
    rels = build_relation_ids(n)
    for _ in range(40):
        Bs = rng.choice([-1, 1], size=(n, n)); np.fill_diagonal(Bs, 1)
        a = rng.choice([-1, 1], size=n)
        tables = {}
        for use in (("R2",), ("R7",)):
            eco = Ecology(a, Bs); chi = eco.initial_table(); propagate(chi, rels)
            try:
                eco.apply_rules(chi, [], use=use)
                tables[use] = chi
            except Conflict:
                tables[use] = None
        r2, r7 = tables[("R2",)], tables[("R7",)]
        if r7 is None:
            continue
        assert r2 is not None
        assert all(v7 == v2 or v2 == 0 for v7, v2 in zip(r7, r2))
