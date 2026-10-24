"""
Relabelling symmetry: the stabilizer acts on sign tables, the rules commute with it, and
complete(symmetry=True) keeps the verdicts and finds every completion up to the stabilizer.
"""
import sys
from pathlib import Path

import pytest

from signrules.patterns import ELEVEN  # noqa: E402

from signrules import Ecology, complete, completion_dot  # noqa: E402
from signrules.engine import Engine  # noqa: E402
from signrules.patterns import random_pattern  # noqa: E402
from signrules.symmetry import (pattern_stabilizer, stabilizer_actions, expand_orbits,  # noqa: E402
                                    orbit_representatives)


def key(t):
    return tuple(t)


def test_stabilizer_orders():
    assert len(pattern_stabilizer(*ELEVEN["S-a obligate mutualism"])) == 24
    assert len(pattern_stabilizer(*ELEVEN["A-a Lemma 3.8"])) == 2


@pytest.mark.parametrize("name", ["S-a obligate mutualism", "A-a Lemma 3.8", "A-e competitor + 3-cycle (Thm 3.5)"])
def test_root_tables_are_invariant(name):
    eco = Ecology(*ELEVEN[name])
    G = stabilizer_actions(eco)
    assert G
    chi = eco.initial_table()
    assert all(g.fixes(chi) for g in G)
    Engine(eco, ("R1", "R2n")).fixpoint(chi, [])     # the root fixed point commutes with G
    assert all(g.fixes(chi) for g in G)


@pytest.mark.parametrize("name,use", [
    ("A-a Lemma 3.8", ("R1", "R2n")),
    ("A-c Lemma 3.10", ("R1", "R2n")),
    ("A-e competitor + 3-cycle (Thm 3.5)", ("R1", "R2n")),
    ("S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)", ("R1", "R2n")),
])
def test_one_completion_per_orbit(name, use):
    eco = Ecology(*ELEVEN[name])
    full = complete(eco, use=use, keep_derivations=False)
    sym = complete(eco, use=use, symmetry=True)
    assert sym.symmetry
    assert {key(t) for t in expand_orbits(sym.completions, sym.symmetry)} == {key(t) for t in full.completions}
    assert len(sym.completions) == len(orbit_representatives(full.completions, sym.symmetry))
    assert sym.nodes < full.nodes
    ks = [i for i, d in enumerate(sym.derivations) if any(st[0].startswith("SYM") for st in d.log)]
    assert ks or sym.symmetry_cut
    for k in ks[:1]:
        assert "digraph" in completion_dot(eco, sym, k)  # SYM steps render like any other rule


def test_branching_kill_is_kept():
    """An impossibility that needs branching stays impossible, with a smaller tree."""
    a, B = random_pattern(4, 176)
    eco = Ecology(a, B)
    use = ("R1", "R2n", "R5", "R6")
    full = complete(eco, use=use, keep_derivations=False)
    sym = complete(eco, use=use, keep_derivations=False, symmetry=True)
    assert full.root_conflict is None and not full.completions
    assert not sym.completions and sym.nodes < full.nodes


@pytest.mark.parametrize("seed", [21, 42, 181])
def test_random_patterns(seed):
    eco = Ecology(*random_pattern(4, seed))
    full = complete(eco, use=("R1", "R2n"), keep_derivations=False)
    sym = complete(eco, use=("R1", "R2n"), keep_derivations=False, symmetry=True)
    assert {key(t) for t in expand_orbits(sym.completions, sym.symmetry)} == {key(t) for t in full.completions}
