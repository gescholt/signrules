"""
Relabelling symmetry of a sign pattern, acting on Plücker sign tables.

A permutation π of the species with  a_{π(i)} = a_i  and  b_{π(i)π(j)} = b_ij  (signs) fixes the
pattern.  It maps a realisation (a, B) to (Πa, ΠBΠᵀ), which is feasible and stable iff (a, B) is,
so it maps realisable sign tables to realisable sign tables.  On minors,

    det B'_{π(P), π(Q)} = sgn(π|P) · sgn(π|Q) · det B_{P,Q},

with sgn(π|P) the sign of the permutation sorting (π(p) for p in P); the a-factor of the Plücker
coordinate is unchanged because π fixes sign a.  Every rule is stated for all index sets alike, so
the rules commute with this action, and so do the fixed points, the conflicts and the set of
completions: g·Comp(T) = Comp(g·T).

``Action`` stores one group element as the image id and sign of every basis id:
(g·chi)[image[S]] = sign[S] · chi[S].  Tables are flat lists indexed by basis id (plucker_index.py).
"""

from itertools import permutations
from typing import List, Sequence, Tuple

from .plucker_index import UNKNOWN, perm_sign


def pattern_stabilizer(a: Sequence[int], B) -> List[Tuple[int, ...]]:
    """All permutations π of the species fixing (sign a, sign B), the identity first."""
    n = len(a)
    return [p for p in permutations(range(n))
            if all(a[p[i]] == a[i] for i in range(n))
            and all(B[p[i]][p[j]] == B[i][j] for i in range(n) for j in range(n) if i != j)]


class Action:
    """The action of one relabelling π on the bases of Gr(n, 2n) of an ``Ecology``."""

    def __init__(self, eco, perm: Sequence[int]):
        self.perm = tuple(perm)
        m = eco.idx.m
        self.image: List[int] = [0] * m
        self.sign: List[int] = [0] * m
        for S, (_, P, Q, _) in enumerate(eco.idx.decomps):
            _, c = eco._minor_key[(P, Q)]               # sign det B_PQ = c · chi[S]
            P2 = [perm[i] for i in P]
            Q2 = [perm[j] for j in Q]
            S2, c2 = eco._minor_key[(tuple(sorted(P2)), tuple(sorted(Q2)))]
            self.image[S] = S2
            self.sign[S] = c * c2 * perm_sign(P2) * perm_sign(Q2)

    def apply(self, chi: List[int]) -> List[int]:
        """The table g·chi."""
        out = [UNKNOWN] * len(chi)
        image, sign = self.image, self.sign
        for S, v in enumerate(chi):
            out[image[S]] = sign[S] * v
        return out

    def fixes(self, chi: List[int]) -> bool:
        """g·chi == chi."""
        image, sign = self.image, self.sign
        for S, v in enumerate(chi):
            if chi[image[S]] != sign[S] * v:
                return False
        return True

    def literal(self, S: int, s: int) -> Tuple[int, int]:
        """Image of the literal chi[S] = s (S an id)."""
        return self.image[S], self.sign[S] * s


def stabilizer_actions(eco) -> List[Action]:
    """The non-identity elements of the pattern stabilizer of ``eco``, as actions on its bases."""
    return [Action(eco, p) for p in pattern_stabilizer(eco.a, eco.B)[1:]]


def orbit_representatives(tables: List[List[int]], actions: List[Action]) -> List[List[int]]:
    """One table per orbit (the first met), in input order; ``actions`` must form a group with the identity."""
    seen, reps = set(), []
    for t in tables:
        key = tuple(t)
        if key in seen:
            continue
        reps.append(t)
        seen.add(key)
        for g in actions:
            seen.add(tuple(g.apply(t)))
    return reps


def expand_orbits(tables: List[List[int]], actions: List[Action]) -> List[List[int]]:
    """The union of the orbits of ``tables``, without repetitions."""
    seen, out = set(), []
    for t in tables:
        for u in [t] + [g.apply(t) for g in actions]:
            key = tuple(u)
            if key not in seen:
                seen.add(key)
                out.append(u)
    return out
