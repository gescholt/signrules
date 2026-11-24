"""The certificates of possibility: together with the root kills they decide every candidate of the screens."""
from fractions import Fraction

import numpy as np
import pytest

from signrules.data import census, possible, root_kills, screen_records
from signrules.possible import charpoly, check_matrix, check_witness, hurwitz_stable, witness_matrix


@pytest.mark.parametrize("n, kills, certified", [(5, 24, 59), (6, 41, 102)])
def test_every_candidate_is_killed_or_certified_possible(n, kills, certified):
    screen = {r["code"]: r for r in screen_records(n)}
    killed = {r["code"] for r in root_kills(n)}
    certs = possible(n)
    assert len(killed) == kills and len(certs) == certified and len(screen) == kills + certified
    assert {r["code"] for r in certs} == set(screen) - killed
    for r in certs:
        assert (r["a"], r["B"]) == (screen[r["code"]]["a"], screen[r["code"]]["B"])
        assert check_witness(r), r["code"]


def test_four_species_census_is_decided():
    rows = census(4)
    assert len(rows) == 3044
    assert sum(r["status"] == "killed_root" for r in rows) == 715
    certs = [r for r in rows if r["status"] == "possible"]
    assert len(certs) == 2329 and all(check_witness(r) for r in certs)


def test_a_tampered_witness_fails():
    r = possible(5)[0]
    bad = dict(r, B=[row[:] for row in r["B"]])
    bad["B"][0][1] = -bad["B"][0][1]                       # the matrix no longer has the signs of the pattern
    assert not check_matrix(r["a"], r["B"], witness_matrix(bad["B"], r["U"], r["T"]))
    assert not check_witness(dict(r, a=[-s for s in r["a"]]))   # row sums of the wrong signs


def test_charpoly_and_stability_against_numpy():
    rng = np.random.default_rng(0)
    for _ in range(40):
        n = int(rng.integers(2, 7))
        C = [[Fraction(int(rng.integers(-9, 10)), int(rng.integers(1, 6))) for _ in range(n)] for _ in range(n)]
        A = np.array(C, dtype=float)
        assert np.allclose([float(c) for c in charpoly(C)], np.poly(-A))
        ev = np.linalg.eigvals(A).real.min()
        if abs(ev) > 1e-6:
            assert hurwitz_stable(C) == (ev > 0)
