"""Tests for signrules/replay.py: independent replay of impossibility certificates."""

import numpy as np
import pytest

from signrules import Ecology, complete, replay_kill, replay_certificate
from signrules.certificate import certificate_steps
from signrules.replay import Replayer


def _kill(a, B):
    eco = Ecology(a, B)
    res = complete(eco, max_nodes=0, max_completions=1, keep_derivations=False)
    assert res.root_conflict is not None
    steps, prem, _ = certificate_steps(eco, res, backward_slice=True)
    return eco, res, steps, prem


def test_initial_table_matches_the_search():
    rng = np.random.default_rng(3)
    for n in (3, 4):
        for _ in range(5):
            a = list(rng.choice([1, -1], size=n))
            B = rng.choice([1, -1], size=(n, n))
            np.fill_diagonal(B, 1)
            eco, rp = Ecology(a, B), Replayer(a, B)
            chi = eco.initial_table()
            for S, v in eco.idx.to_dict(chi).items():
                assert rp.chi(S) == v, (S, v, rp.chi(S))


def test_the_eleven_n4_kills_replay():
    from signrules.patterns import ELEVEN
    for name, (a, B) in ELEVEN.items():
        rep = replay_kill(Ecology(a, B, name))
        assert rep.ok, (name, rep.message)
        assert rep.steps >= 0 and rep.rules


def test_obligate_3cycle_and_balanced_patterns_replay():
    # the obligate 3-cycle (R2/R7 territory) and a balanced symmetric pattern (R3)
    assert replay_kill(Ecology([-1, -1, -1], [[1, -1, 1], [1, 1, -1], [-1, 1, 1]])).ok
    B = np.array([[1, -1, 1, 1], [-1, 1, 1, 1], [1, 1, 1, -1], [1, 1, -1, 1]])
    rep = replay_kill(Ecology([-1, -1, -1, -1], B))
    assert rep.ok, rep.message


def test_r8_parity_pattern_replays():
    from signrules.patterns import parity_pattern
    a, B = parity_pattern(5, [[0, 1], [2, 3]], [4])
    rep = replay_kill(Ecology(list(a), np.array(B)))
    assert rep.ok, rep.message
    assert "R8" in rep.rules


def test_tampered_certificates_fail():
    from signrules.patterns import ELEVEN
    for name, (a, B) in ELEVEN.items():
        eco, res, steps, prem = _kill(a, B)
        assert replay_certificate(a, B, steps, res.root_conflict, prem).ok
        for k in range(len(steps)):
            S, sign, why, p = steps[k]
            bad = list(steps)
            bad[k] = (S, -sign, why, p)
            assert not replay_certificate(a, B, bad, res.root_conflict, prem).ok, (name, k, why)
        # dropping all steps: the conflict is no longer supported (unless it reads the pattern alone)
        if steps:
            rep = replay_certificate(a, B, [], res.root_conflict, prem)
            assert not rep.ok or res.root_conflict.startswith("R8")


@pytest.mark.parametrize("n, count", [(5, 24), (6, 41)])
def test_every_minimal_impossible_pattern_replays_and_resists_tampering(n, count):
    """The certificates the classification rests on: each root kill of the n-species screen replays, and
    flipping the sign of any single step makes the replay fail."""
    from signrules.data import root_kills
    kills = root_kills(n)
    assert len(kills) == count
    for r in kills:
        a, B = r["a"], np.array(r["B"])
        eco, res, steps, prem = _kill(a, B)
        rep = replay_certificate(a, B, steps, res.root_conflict, prem)
        assert rep.ok, (r["code"], rep.message)
        for k, (S, sign, why, p) in enumerate(steps):
            bad = list(steps)
            bad[k] = (S, -sign, why, p)
            assert not replay_certificate(a, B, bad, res.root_conflict, prem).ok, (r["code"], k, why)


def test_a_possible_pattern_has_nothing_to_replay():
    rep = replay_kill(Ecology([-1, 1], [[1, -1], [1, 1]]))
    assert not rep.ok and "no kill" in rep.message


def test_a_wrong_conflict_is_rejected():
    a, B = [-1, -1, -1], np.array([[1, -1, 1], [1, 1, -1], [-1, 1, 1]])
    eco, res, steps, prem = _kill(a, B)
    bogus = "R8 real eigenvalue s=(1, 1, 1) σ=(1, 2, 0): nonsense"
    assert not replay_certificate(a, B, steps, bogus, []).ok
    assert not replay_certificate(a, B, steps, "R1 det B > 0: det B_(0, 1, 2),(0, 1, 2) must be > 0 but is < 0", []).ok


def test_malformed_steps_fail_the_replay_instead_of_crashing():
    from signrules.patterns import ELEVEN
    from signrules.replay import Replayer, ReplayError
    a, B = ELEVEN["S-c {0,1}|{2,3}, {2,3} obligate, a01=(+,+)"]
    eco, res, steps, prem = _kill(a, B)
    k, (S, sign, why, p) = next((k, st) for k, st in enumerate(steps) if st[2].startswith("R2"))
    bad = list(steps)
    bad[k] = (S, sign, why.replace("x*_2", "x*_9"), p)            # j not in J
    rep = replay_certificate(a, B, bad, res.root_conflict, prem)
    assert not rep.ok and "malformed" in rep.message
    rep = replay_certificate(a, B, [steps[0][:3]] + list(steps[1:]), res.root_conflict, prem)   # wrong arity
    assert not rep.ok and "malformed" in rep.message
    with pytest.raises(ReplayError):                                 # a GP step without its relation
        Replayer(a, B).verify_step((0, 4, 6, 7), 1, "Grassmann–Plücker 3-term relation", [])
