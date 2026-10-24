"""
The flat sign table (a list indexed by basis id, UNKNOWN = 0) and the search statistics it must
reproduce.  The golden numbers were recorded with the earlier dict-based table (feasible_ecologies
history, 2026-10) before the change: same propagators, same priority levels, same branching order, hence the same trees.
"""
import json

import numpy as np
from signrules.data import DATA_DIR
import pytest

from signrules.patterns import ELEVEN

from signrules import (Ecology, complete, UNKNOWN, index, build_relations,
                           build_relation_ids, propagate, violates)

# name -> ((root kill, nodes, conflicts, completions) under R1,R2n, the same under all rules)
GOLDEN_N4 = {
    'S-a obligate mutualism': ((False, 4, 5, 0), (True, 0, 0, 0)),
    'S-b {0}|{1,2,3}, triple obligate': ((False, 4, 5, 0), (True, 0, 0, 0)),
    'S-c {0,1}|{2,3}, {2,3} obligate, a01=(+,+)': ((False, 4, 5, 0), (True, 0, 0, 0)),
    'S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)': ((False, 131, 4, 64), (True, 0, 0, 0)),
    'S-e {0,1}|{2,3}, both obligate': ((False, 131, 4, 64), (True, 0, 0, 0)),
    'A-a Lemma 3.8': ((False, 3993, 2, 1996), (True, 0, 0, 0)),
    'A-b Lemma 3.9': ((False, 31025, 2, 15512), (True, 0, 0, 0)),
    'A-c Lemma 3.10': ((False, 3060, 5, 1528), (True, 0, 0, 0)),
    'A-d obligate 4-cycle (Thm 3.4)': ((False, 31021, 6, 15508), (True, 0, 0, 0)),
    'A-e competitor + 3-cycle (Thm 3.5)': ((False, 3996, 5, 1996), (True, 0, 0, 0)),
    'A-f super-predator + 3-cycle (Thm 3.6)': ((False, 31027, 4, 15512), (True, 0, 0, 0)),
}
# code -> (root kill, nodes, conflicts, completions, open signs after root propagation), all rules,
# first completion, node limit 20000 (as scripts/extension_screen.py)
GOLDEN_N5 = {
    1162680: (False, 51, 0, 1, 125),
    1539539: (False, 67, 0, 1, 141),
    1555808: (False, 81, 6, 1, 110),
    1031600: (False, 59, 0, 1, 103),
    3636688: (False, 69, 0, 1, 125),
}
# the same for five open n = 6 patterns of signrules/data/n6/extension_screen.jsonl (recorded with
# scripts/rule_ablation.py, variant "all"; the full 237-pattern record is signrules/data/rule_ablation/all.jsonl)
GOLDEN_N6 = {
    15518791496: (False, 188, 29, 1, 723),
    1132328688: (False, 123, 0, 1, 588),
    16906974688: (False, 165, 4, 1, 601),
    33252308816: (False, 147, 2, 1, 507),
    58455936: (False, 162, 5, 1, 657),
}


def stats(r):
    return r.root_conflict is not None, r.nodes, r.conflicts, len(r.completions)


@pytest.mark.parametrize("name", [n for n in GOLDEN_N4 if GOLDEN_N4[n][0][1] < 5000])
def test_golden_n4(name):
    a, B = ELEVEN[name]
    weak, full = GOLDEN_N4[name]
    assert stats(complete(Ecology(a, B), use=("R1", "R2n"), keep_derivations=False)) == weak
    assert stats(complete(Ecology(a, B), keep_derivations=False)) == full


@pytest.mark.parametrize("n,code", [(5, c) for c in sorted(GOLDEN_N5)] + [(6, c) for c in sorted(GOLDEN_N6)])
def test_golden_screen_search(n, code):
    golden = GOLDEN_N5 if n == 5 else GOLDEN_N6
    rows = [json.loads(line) for line in open(DATA_DIR / f"n{n}/extension_screen.jsonl")]
    rec = next(r for r in rows if r["code"] == code)
    r = complete(Ecology(rec["a"], np.array(rec["B"])), max_completions=1, max_nodes=20000, keep_derivations=False)
    assert stats(r) + (r.root_unknowns,) == golden[code]
    table = r.completions[0]
    assert len(table) == r.completions[0].count(1) + table.count(-1) == index(n).m
    eco = Ecology(rec["a"], np.array(rec["B"]))
    assert eco.minor(table, tuple(range(n)), tuple(range(n))) == 1                # R1: det B > 0


@pytest.mark.parametrize("n", [3, 4])
def test_table_conversions(n):
    idx = index(n)
    assert idx.m == len(idx.bases) and all(idx.id[S] == i for i, S in enumerate(idx.bases))
    eco = Ecology([1, -1, 1, -1][:n], np.where(np.eye(n) > 0, 1, -1))
    chi = eco.initial_table()
    assert len(chi) == idx.m and set(chi) <= {-1, 0, 1}
    d = idx.to_dict(chi)
    assert list(d) == idx.bases and idx.from_dict(d) == chi
    assert UNKNOWN == 0 and idx.from_dict({}) == [0] * idx.m
    # the id relations are the basis relations, term by term
    for rel, irel in zip(build_relations(n), build_relation_ids(n)):
        assert [(idx.id[k1], idx.id[k2], c) for k1, k2, c in rel] == list(irel)


def test_sweep_matches_trigger_on_the_eleven():
    for name, (a, B) in ELEVEN.items():
        if GOLDEN_N4[name][0][1] >= 5000:
            continue
        runs = [complete(Ecology(a, B), use=("R1", "R2n"), engine=e, branching="first", keep_derivations=False)
                for e in ("sweep", "trigger")]
        assert runs[0].root_unknowns == runs[1].root_unknowns
        assert sorted(map(tuple, runs[0].completions)) == sorted(map(tuple, runs[1].completions))
        for t in runs[0].completions:
            assert not violates(t, build_relation_ids(4))
        # the sweep propagation logs bases, not ids
        chi, log = Ecology(a, B).initial_table(), []
        propagate(chi, build_relation_ids(4), log, index(4).bases)
        assert all(isinstance(S, tuple) and all(isinstance(p, tuple) for p in prem) for S, _, prem in log)
