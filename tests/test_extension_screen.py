"""The extension screen (signrules.screen): canonical forms and the n=5 candidate list."""

import sys
import random
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

from signrules.screen import base_classes, candidates
from signrules.patterns import canonical, decode, encode, key


def test_canonical_is_relabelling_invariant():
    rng = random.Random(1)
    for n in (5, 6):
        for _ in range(10):
            code = rng.getrandbits(n * n)
            a, B = decode(code, n)
            p = rng.sample(range(n), n)
            a2 = [0] * n
            B2 = np.eye(n, dtype=int)
            for i in range(n):
                a2[p[i]] = a[i]
                for j in range(n):
                    B2[p[i], p[j]] = B[i, j]
            assert decode(encode(a, B), n)[0] == a
            assert canonical([code], n)[0] == canonical([encode(a2, B2)], n)[0]


def test_n5_candidates_contain_the_known_kills():
    base = base_classes(5, "census", workers=4)
    assert len(base) == 11                                   # the eleven non-trivial n=4 kills
    _, n_reps, cands = candidates(5, base)
    assert (n_reps, len(cands)) == (2142, 83)
    found = {key(*decode(c, 5)) for c in cands}
    known = {key(a, B) for a, B in base_classes(6, "known", workers=1)}
    assert len(known) == 24 and known <= found
