"""
complete(root_only=...) and complete(stats=True): the ablation switches used by scripts/rule_ablation.py.

root_only keeps the propagators registered (same watcher counts, same branching order) and stops
waking them below the root; stats counts runs, signs and conflicts per rule.  Neither may change
the search when off, and root_only may never reject a realized sign table.
"""
import json

import numpy as np
from signrules.data import DATA_DIR
import pytest

from signrules.patterns import ELEVEN

from signrules import Ecology, complete, Conflict
from signrules.engine import Engine
from signrules.patterns import random_pattern
from test_flat_table import GOLDEN_N5, stats
from test_soundness import random_feasible_stable


def n5_pattern(code):
    rows = [json.loads(line) for line in open(DATA_DIR / "n5/extension_screen.jsonl")]
    rec = next(r for r in rows if r["code"] == code)
    return rec["a"], np.array(rec["B"])


def test_engine_keeps_order_and_watchers_under_root_only():
    a, B = n5_pattern(1555808)
    full, part = Engine(Ecology(a, B)), Engine(Ecology(a, B), root_only=("R4", "R6"))
    assert full.order == part.order and full.watchers == part.watchers and full.level == part.level
    assert set(part.rule_of) == {"GP", *Ecology.ALL_RULES}
    assert all(bool(part.root_only[pid]) == (part.rule_of[pid] in ("R4", "R6")) for pid in range(len(part.run)))
    assert not any(full.root_only)
    with pytest.raises(ValueError):
        Engine(Ecology(a, B), use=("R1", "R2n"), root_only=("R6",))


@pytest.mark.parametrize("code", sorted(GOLDEN_N5))
def test_stats_and_empty_root_only_leave_the_search_unchanged(code):
    a, B = n5_pattern(code)
    kw = dict(max_completions=1, max_nodes=20000, keep_derivations=False)
    plain = complete(Ecology(a, B), **kw)
    assert stats(plain) + (plain.root_unknowns,) == GOLDEN_N5[code]
    assert plain.rule_stats is None
    r = complete(Ecology(a, B), root_only=(), stats=True, **kw)
    assert stats(r) == stats(plain) and r.completions == plain.completions and r.root_log == plain.root_log
    assert set(r.rule_stats) <= set(Ecology.ALL_RULES)
    for counts in r.rule_stats.values():
        for where in ("root", "below"):
            c = counts[where]
            assert c["runs"] >= c["conflicts"] and c["runs"] >= 0 and c["signs"] >= 0
    assert sum(c["root"]["signs"] for c in r.rule_stats.values()) == len(r.root_log)


@pytest.mark.parametrize("code", sorted(GOLDEN_N5))
def test_root_only_r4_r6_keeps_the_root_and_runs_them_only_there(code):
    a, B = n5_pattern(code)
    kw = dict(max_completions=1, max_nodes=20000, keep_derivations=False, stats=True)
    full = complete(Ecology(a, B), **kw)
    part = complete(Ecology(a, B), root_only=("R4", "R6"), **kw)
    assert part.root_log == full.root_log and part.root_unknowns == full.root_unknowns
    for rule in ("R4", "R6"):
        assert part.rule_stats[rule]["root"] == full.rule_stats[rule]["root"]
        assert part.rule_stats[rule]["below"] == dict(runs=0, signs=0, conflicts=0)
    # weaker pruning below the root can only enlarge the tree
    assert part.nodes >= full.nodes and len(part.completions) == len(full.completions) == 1


@pytest.mark.parametrize("name", [n for n, (a, B) in ELEVEN.items()])
def test_root_only_keeps_the_eleven_root_kills(name):
    a, B = ELEVEN[name]
    full = complete(Ecology(a, B), keep_derivations=False)
    part = complete(Ecology(a, B), root_only=("R4", "R6"), keep_derivations=False)
    assert full.root_conflict is not None
    assert part.root_conflict == full.root_conflict and part.root_log == full.root_log


def test_root_only_and_stats_need_the_trigger_engine():
    a, B = random_pattern(3, 0)
    with pytest.raises(ValueError):
        complete(Ecology(a, B), engine="sweep", root_only=("R6",))
    with pytest.raises(ValueError):
        complete(Ecology(a, B), engine="sweep", stats=True)


@pytest.mark.parametrize("n", [3, 4])
def test_root_only_never_rejects_a_true_chirotope(n):
    """Along the branch that follows a realized sign table, propagation with root_only never
    conflicts and never forces a sign different from the truth (the walk of one dfs path)."""
    rng = np.random.default_rng(7)
    for trial in range(15):
        a, B = random_feasible_stable(n, rng)
        eco = Ecology(np.sign(a).astype(int), np.sign(B).astype(int))
        num = eco.idx.numeric_chi(a, B)
        if any(v == 0 for v in num):
            continue
        engine = Engine(eco, root_only=("R4", "R6"))
        chi = eco.initial_table()
        try:
            engine.fixpoint(chi, [])
            while True:
                assert all(c in (0, v) for c, v in zip(chi, num))
                open_ids = [S for S in engine.order if chi[S] == 0]
                if not open_ids:
                    break
                S = open_ids[0]
                chi[S] = num[S]
                engine.fixpoint(chi, [], None, [S])
        except Conflict as c:
            pytest.fail(f"true system rejected: {c}")
        assert chi == list(num)


def test_root_only_completions_contain_the_full_ones_on_random_patterns():
    for seed in range(8):
        a, B = random_pattern(4, seed)
        kw = dict(branching="first", max_nodes=200, keep_derivations=False)
        full = complete(Ecology(a, B), **kw)
        part = complete(Ecology(a, B), root_only=("R4", "R6"), **kw)
        assert (full.root_conflict is None) == (part.root_conflict is None), seed
        if full.stopped or part.stopped:
            continue
        assert set(map(tuple, full.completions)) <= set(map(tuple, part.completions)), seed
