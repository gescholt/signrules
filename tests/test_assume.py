"""
complete(assume=...) and the exact sign-correlation test of scripts/sign_correlations_exact.py: starting
from a table with extra signs, and deciding whether two Plücker coordinates are correlated on every
completion without enumerating them.
"""
import itertools

import numpy as np
import pytest

from sign_correlations_exact import analyse, prove_pair

from signrules import Ecology, complete

HUB = ([-1, -1, -1], [[1, 1, -1], [1, 1, -1], [-1, -1, 1]])     # the obligate hub (species 2) with competing satellites
F, G = ((0, 2), (0, 1)), ((1, 2), (0, 1))


def hub_ids():
    eco = Ecology(*HUB)
    (Sf, cf), (Sg, cg) = eco._minor_key[F], eco._minor_key[G]
    return eco, Sf, cf, Sg, cg


def test_assumed_signs_hold_in_every_completion():
    eco, Sf, cf, Sg, cg = hub_ids()
    r = complete(eco, assume={Sf: 1}, keep_derivations=False)
    assert r.assumed == {Sf: 1} and r.completions and all(t[Sf] == 1 for t in r.completions)
    full = complete(eco, keep_derivations=False).completions
    assert sorted(map(tuple, r.completions)) == sorted(tuple(t) for t in full if t[Sf] == 1)


def test_assumption_against_the_initial_table_is_a_root_conflict():
    eco = Ecology(*HUB)
    chi = eco.initial_table()
    S = next(i for i, v in enumerate(chi) if v != 0)
    r = complete(eco, assume={S: -chi[S]}, keep_derivations=False)
    assert r.root_conflict and "assumption" in r.root_conflict and r.nodes == 0
    with pytest.raises(ValueError):
        complete(eco, assume={S: 2})


def test_hub_minors_are_correlated_on_every_completion():
    eco, Sf, cf, Sg, cg = hub_ids()
    # the two coordinates: only two of the four sign cells survive propagation
    r = complete(eco, decide=[Sf, Sg], keep_derivations=False)
    cells = {(t[Sf], t[Sg]) for t in r.completions}
    assert len(cells) == 2 and all(cf * s == cg * t for s, t in cells)      # sgn f = sgn g
    eps = cf * cg
    assert prove_pair(eco, Sf, Sg, eps, 10 ** 6)[:2] == ("proved", "propagation")
    assert prove_pair(eco, Sf, Sg, -eps, 10 ** 6)[0] == "refuted"             # the opposite relation fails
    # the forbidden combinations die at the root, the allowed ones have completions
    for s, t in itertools.product((1, -1), repeat=2):
        r = complete(eco, assume={Sf: s, Sg: t}, max_completions=1, keep_derivations=False)
        assert bool(r.completions) == (cf * s == cg * t) and (r.root_conflict is None) == bool(r.completions)


def test_analyse_matches_the_enumeration_on_n3():
    """On every n = 3 pattern the exact test and the enumerated completions agree on the correlation classes."""
    for a in itertools.product((1, -1), repeat=3):
        for off in itertools.product((1, -1), repeat=6):
            B = np.eye(3, dtype=int)
            B[~np.eye(3, dtype=bool)] = off
            rec = analyse(list(a), B, sample=None, max_nodes=10 ** 6)
            if rec["status"] != "open":
                continue
            assert rec["exhaustive"] and rec["undecided_pairs"] == 0 and rec["refuted_pairs"] == 0
            assert rec["disconnected"] == (list(a) == [-1, -1, -1] and rec["candidate_classes"] == 1)
            if rec["disconnected"]:
                cls = rec["classes"][0]
                assert cls["size"] == 2 and cls["plus"] == cls["minus"] == 6 and cls["status"] == "proved"


@pytest.mark.parametrize("sample", [1, 3])
def test_screen_then_prove_on_a_truncated_sample(sample):
    """With one or three sampled completions the hub pair is constant on the sample (so the freeness step
    must find it) and spurious correlations appear (so the proof step must refute them)."""
    a, B = HUB
    rec = analyse(a, B, sample=sample, max_nodes=10 ** 6)
    assert rec["sampled"] == sample and not rec["exhaustive"]
    assert rec["undecided_free"] == [] and rec["forced"] + rec["free"] == Ecology(a, B).idx.m
    assert rec["disconnected"] and rec["undecided_pairs"] == 0
    proved = [p for c in rec["classes"] for p in c["pairs"] if p["verdict"] == "proved"]
    assert [sorted(p["ids"]) for p in proved] == [sorted(hub_ids()[1::2])]
    full = analyse(a, B, sample=None, max_nodes=10 ** 6)
    assert rec["free"] == full["free"] and full["forced"] == rec["forced"]
