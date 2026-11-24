"""complete(assume=...): starting the search from a table with extra signs."""

import pytest

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
