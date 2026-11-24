"""
Regression tests for the search and the derivation graphs: the documented results of
docs/SIGNRULES.md §4 (the eleven n = 4 kills, the R1,R2n ablation counts) and the
Graphviz renderers.  Everything here runs in well under a second per case.
"""

import pytest

from signrules.patterns import ELEVEN

from signrules import Ecology, complete, certificate_dot, completion_dot, UNKNOWN


@pytest.mark.parametrize("name", list(ELEVEN))
def test_eleven_killed_at_root(name):
    a, B = ELEVEN[name]
    res = complete(Ecology(a, B, name))
    assert res.root_conflict is not None, name
    assert res.root_conflict_premises
    assert "digraph" in certificate_dot(Ecology(a, B), res)


@pytest.mark.parametrize("name,count", [
    ("S-a obligate mutualism", 0),
    ("S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)", 64),
    ("S-e {0,1}|{2,3}, both obligate", 64),
])
def test_ablation_counts_r1_r2n(name, count):
    """Weak-completion counts of the algorithm of the Annals paper (R1 + R2n), as quoted in Çelik–Isoldi–Portakal–Zucal."""
    a, B = ELEVEN[name]
    res = complete(Ecology(a, B), use=("R1", "R2n"))
    assert res.root_conflict is None
    assert len(res.completions) == count
    assert len(res.derivations) == count


@pytest.mark.parametrize("use", [("R1", "R2n"), ("R1", "R2n", "R5", "R6"), Ecology.ALL_RULES])
@pytest.mark.parametrize("name", ["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)", "A-c Lemma 3.10"])
def test_trigger_engine_matches_sweep(name, use):
    """Same root fixed point and the same completions for both engines and both branching orders."""
    a, B = ELEVEN[name]
    runs = [complete(Ecology(a, B), use=use, engine=e, branching=b, keep_derivations=False)
            for e, b in (("sweep", "first"), ("trigger", "first"), ("trigger", "constrained"), ("trigger", "conflict"))]
    assert len({r.root_conflict is None for r in runs}) == 1
    assert len({r.root_unknowns for r in runs}) == 1
    tables = [sorted(map(tuple, r.completions)) for r in runs]
    assert tables[0] == tables[1] == tables[2] == tables[3]


def test_trigger_engine_matches_sweep_on_random_patterns():
    from signrules.patterns import random_pattern
    for seed in range(12):
        a, B = random_pattern(4, seed)
        eco = Ecology(a, B)
        r0 = complete(eco, engine="sweep", branching="first", max_nodes=60, keep_derivations=False)
        r1 = complete(eco, engine="trigger", branching="first", max_nodes=60, keep_derivations=False)
        assert (r0.root_conflict is None) == (r1.root_conflict is None), seed
        assert r0.root_unknowns == r1.root_unknowns
        assert r0.completions == r1.completions          # same branching order: same tree, same leaves


def test_completion_derivation_accounts_for_every_sign():
    a, B = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]
    eco = Ecology(a, B)
    res = complete(eco, use=("R1", "R2n"), max_completions=3)
    bases = eco.idx.bases
    initial = eco.idx.to_dict(eco.initial_table())
    for table, d in zip(res.completions, res.derivations):
        assert all(v != UNKNOWN for v in table)
        # every sign of the completion is either initial, a branch choice, or forced by a logged step
        explained = {S for S, v in initial.items() if v != UNKNOWN} | {S for S, _ in d.branches}
        explained |= {S for S, *_ in d.gp_log}
        explained |= {bases[eco._minor_key[(tuple(P), tuple(Q))][0]] for _, P, Q, _, _ in d.log}
        assert explained == set(bases)
        for S, sign in d.branches:
            assert table[eco.idx.id[S]] == sign
    dot = completion_dot(eco, res, index=0)
    d0 = res.derivations[0]
    assert "root propagation" in dot and "GP · 3-term relation" in dot
    assert all(f"branch {k}" in dot for k in range(1, len(d0.branches) + 1))       # every branch choice is a node
    assert dot.count("subgraph \"cluster_band") <= len(d0.branches) + 1               # empty branches share a band
    full = completion_dot(eco, res, index=0, edges="all")
    assert dot.count(" -> ") <= full.count(" -> ")                                    # band-local edges by default
    flat = completion_dot(eco, res, index=1, structured=False)
    assert "GP · 3-term relation" in flat and "p_" in flat
    stages = list(res.derivations[0].stages())
    assert stages[0][0] is None and [b for b, *_ in stages[1:]] == res.derivations[0].branches
    assert sum(len(l) for _, l, _ in stages) == len(res.derivations[0].log)
    assert sum(len(g) for _, _, g in stages) == len(res.derivations[0].gp_log)


def test_completion_dot_requires_a_completion():
    a, B = ELEVEN["S-a obligate mutualism"]
    eco = Ecology(a, B)
    with pytest.raises(ValueError):
        completion_dot(eco, complete(eco))
    with pytest.raises(ValueError):
        certificate_dot(eco, complete(eco, use=("R1", "R2n")))


def test_search_tree_records_every_node():
    from signrules import search_tree_text, search_tree_dot
    a, B = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]
    eco = Ecology(a, B)
    res = complete(eco, use=("R1", "R2n"), record_tree=True)
    status = [nd["status"] for nd in res.tree]
    assert status.count("completion") == len(res.completions) == 64
    assert status.count("conflict") == res.conflicts
    assert len(res.tree) == res.nodes + res.conflicts                 # visited nodes + pruned children
    assert all(nd["parent"] is None or res.tree[nd["parent"]]["branch_on"] == nd["basis"] for nd in res.tree)
    text = search_tree_text(eco, res, max_depth=2)
    assert text.startswith("root:") and "branch on" in text
    assert "digraph" in search_tree_dot(eco, res)
    # completion limit: the unsearched subtrees are marked, not mislabelled as completions
    res2 = complete(eco, use=("R1", "R2n"), record_tree=True, max_completions=3)
    assert [nd["status"] for nd in res2.tree].count("completion") == 3
    assert any(nd["status"] == "stopped" for nd in res2.tree)


def test_tree_layout_and_matplotlib_plot(tmp_path):
    from signrules.treeplot import tree_graph, tree_layout, save_search_tree_png
    a, B = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]
    eco = Ecology(a, B)
    res = complete(eco, use=("R1", "R2n"), record_tree=True)
    G, pos = tree_graph(res), tree_layout(res)
    assert set(G) == set(pos) and len(G) == len(res.tree)
    leaves = sorted(pos[i][0] for i in G if G.out_degree(i) == 0)
    assert leaves == [k + 0.5 for k in range(len(leaves))]   # leaves evenly spaced in search order
    for i in G:
        kids = list(G.successors(i))
        if kids:
            assert abs(pos[i][0] - sum(pos[k][0] for k in kids) / len(kids)) < 1e-9   # parent centred
            assert all(pos[k][1] == pos[i][1] - 1 for k in kids)
    out = tmp_path / "tree.png"
    save_search_tree_png(eco, res, str(out))
    assert out.stat().st_size > 1000


@pytest.mark.parametrize("probe", ["root", "all"])
@pytest.mark.parametrize("name", ["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)", "S-a obligate mutualism"])
def test_probing_keeps_the_completions(name, probe):
    """Failed-literal probing is sound: the same completions and never more nodes (S-a, a search kill under
    R1, R2n, becomes a root kill)."""
    a, B = ELEVEN[name]
    r0 = complete(Ecology(a, B), use=("R1", "R2n"), keep_derivations=False)
    r1 = complete(Ecology(a, B), use=("R1", "R2n"), keep_derivations=False, probe=probe)
    assert sorted(map(tuple, r0.completions)) == sorted(map(tuple, r1.completions))
    assert r1.nodes <= r0.nodes
    if name.startswith("S-a"):
        assert r0.root_conflict is None and r0.nodes > 0 and r1.root_conflict.startswith("failed-literal")


def test_seed_changes_the_order_not_the_set():
    # used downstream: random coarse cells via complete(decide=..., max_completions=1, seed=k)
    a = [1, -1, 1, -1]
    B = [[1, -1, 1, 1], [1, 1, -1, 1], [-1, 1, 1, -1], [1, -1, 1, 1]]
    eco = Ecology(a, B)
    ref = sorted(map(tuple, complete(eco, keep_derivations=False).completions))
    firsts = set()
    for seed in range(4):
        r = complete(eco, keep_derivations=False, seed=seed)
        assert sorted(map(tuple, r.completions)) == ref
        one = complete(eco, keep_derivations=False, seed=seed, max_completions=1)
        assert len(one.completions) == 1 and tuple(one.completions[0]) in ref
        firsts.add(tuple(one.completions[0]))
    assert len(firsts) > 1
    with pytest.raises(ValueError):
        complete(eco, seed=0, symmetry=True)
