"""
Bookkeeping between minors of B and Plücker coordinates of M = [diag(a) | B].

Every maximal minor of M is a minor of B dressed with a monomial in a:

    p_S = eps(S) * prod_{j in T_a} a_j * det B_{P,Q},

where S = T_a ∪ (n + Q) is a basis (sorted n-subset of range(2n)), P = [n] \\ T_a,
and eps(S) = (-1)^{sum(T_a) - |T_a|(|T_a|-1)/2} (Laplace expansion along the
diag(a) columns, 0-based indices).  Hence the Plücker coordinates of M are
*all* minors of B of all sizes, with |rows| = |cols|, up to known signs.

A *sign table* is a flat list of ints indexed by basis id (the position of the basis
in ``PluckerIndex.bases``): +1, -1, or UNKNOWN = 0.  Ids and bases translate through
``PluckerIndex.id`` and ``PluckerIndex.bases``; ``to_dict`` / ``from_dict`` convert a
table to and from the dict {basis: sign} used in displays and tests.
"""

from functools import lru_cache
from itertools import combinations, permutations
from typing import Dict, List, Tuple, Sequence, Iterable

import numpy as np

UNKNOWN = 0
Basis = Tuple[int, ...]
Table = List[int]


def perm_sign(seq: Sequence[int]) -> int:
    """Sign of the permutation sorting `seq` (assumes distinct entries)."""
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


@lru_cache(maxsize=None)
def signed_perms(k: int) -> Tuple[Tuple[int, Tuple[int, ...]], ...]:
    """(sign, perm) for every permutation of range(k), in itertools order."""
    return tuple((perm_sign(p), p) for p in permutations(range(k)))


def sgn(x: float, tol: float = 1e-9) -> int:
    return 0 if abs(x) < tol else (1 if x > 0 else -1)


class PluckerIndex:
    """
    Index of the bases of Gr(n, 2n) for the ecological matrix [diag(a) | B].

    Columns 0..n-1 of M are diag(a), columns n..2n-1 are B.  A basis S (a sorted
    n-subset of range(2n)) is decomposed as S = T_a ∪ (n + Q) with P = [n] \\ T_a,
    and  p_S = eps(S) · Π_{j∈T_a} a_j · det B_{P,Q}.  ``decompose`` and ``key``
    translate between bases and minors; ``initial_table`` reads off the signs
    determined by the entries alone.  Signs are +1, -1 or UNKNOWN (0).

    Depends on n only: use the cached ``index(n)`` rather than constructing one per pattern.
    """

    def __init__(self, n: int):
        self.n = n
        self.N = 2 * n
        self.bases: List[Basis] = [tuple(c) for c in combinations(range(self.N), n)]
        self.m = len(self.bases)
        self.id: Dict[Basis, int] = {S: i for i, S in enumerate(self.bases)}
        self.decomps: List[Tuple[Basis, Basis, Basis, int]] = []     # id -> (T_a, P, Q, eps)
        self._decomp: Dict[Basis, Tuple[Basis, Basis, Basis, int]] = {}
        self._minor_key: Dict[Tuple[Basis, Basis], Tuple[Basis, int]] = {}
        self._minor_id: Dict[Tuple[Basis, Basis], Tuple[int, int]] = {}
        for i, S in enumerate(self.bases):
            T_a = tuple(j for j in S if j < n)
            Q = tuple(j - n for j in S if j >= n)
            P = tuple(j for j in range(n) if j not in T_a)
            k = len(T_a)
            eps = -1 if (sum(T_a) - k * (k - 1) // 2) % 2 else 1
            self.decomps.append((T_a, P, Q, eps))
            self._decomp[S] = (T_a, P, Q, eps)
            self._minor_key[(P, Q)] = (S, eps)
            self._minor_id[(P, Q)] = (i, eps)

    # ----- translation between Plücker signs and minor signs -----------------

    def decompose(self, S: Basis):
        """Return (T_a, P, Q, eps) for basis S."""
        return self._decomp[S]

    def key(self, P: Iterable[int], Q: Iterable[int]) -> Tuple[Basis, int]:
        """Basis S and sign eps with p_S = eps * prod_{j not in P} a_j * det B_{P,Q}."""
        return self._minor_key[(tuple(sorted(P)), tuple(sorted(Q)))]

    def key_id(self, P: Iterable[int], Q: Iterable[int]) -> Tuple[int, int]:
        """Like ``key`` but with the basis id instead of the basis."""
        return self._minor_id[(tuple(sorted(P)), tuple(sorted(Q)))]

    def a_factor(self, a_signs: Sequence[int], P: Iterable[int]) -> int:
        s = 1
        Pset = set(P)
        for j in range(self.n):
            if j not in Pset:
                s *= a_signs[j]
        return s

    def minor_sign(self, chi: Table, a_signs: Sequence[int], P: Iterable[int], Q: Iterable[int]):
        """Sign of det B_{P,Q} (rows/cols sorted) from the chi table, or UNKNOWN."""
        i, eps = self.key_id(P, Q)
        return eps * self.a_factor(a_signs, P) * chi[i]

    # ----- tables ---------------------------------------------------------------

    def to_dict(self, chi: Table) -> Dict[Basis, int]:
        """The table as {basis: sign} (UNKNOWN entries included)."""
        return dict(zip(self.bases, chi))

    def from_dict(self, d: Dict[Basis, int]) -> Table:
        """The flat table of a {basis: sign} dict; missing bases are UNKNOWN."""
        chi = [UNKNOWN] * self.m
        for S, v in d.items():
            chi[self.id[S]] = int(v)
        return chi

    # ----- initial table from the ecological sign pattern --------------------

    def initial_table(self, a_signs: Sequence[int], B_signs: np.ndarray) -> Table:
        """
        Sign of every Plücker coordinate that is determined by the signs of the
        entries alone (all permutation terms of the minor agree), UNKNOWN otherwise.
        """
        B = [[int(v) for v in row] for row in np.asarray(B_signs)]
        chi: Table = [UNKNOWN] * self.m
        for i, (T_a, P, Q, eps) in enumerate(self.decomps):
            k = len(P)
            signs = set()
            for t, perm in signed_perms(k):
                for r, c in zip(P, perm):
                    t *= B[r][Q[c]]
                signs.add(t)
                if len(signs) == 2:
                    break
            if len(signs) == 1:
                chi[i] = eps * self.a_factor(a_signs, P) * signs.pop()
        return chi

    # ----- numeric sanity checks ----------------------------------------------

    def numeric_chi(self, a: np.ndarray, B: np.ndarray) -> Table:
        M = np.hstack([np.diag(a), B])
        return [sgn(np.linalg.det(M[:, list(S)])) for S in self.bases]

    def check_conventions(self, seed: int = 0) -> None:
        """Verify p_S = eps * prod a_j * det B_{P,Q} on a random instance."""
        rng = np.random.default_rng(seed)
        a = rng.normal(size=self.n)
        B = rng.normal(size=(self.n, self.n))
        M = np.hstack([np.diag(a), B])
        for S in self.bases:
            T_a, P, Q, eps = self._decomp[S]
            lhs = np.linalg.det(M[:, list(S)])
            sub = B[np.ix_(P, Q)] if P else np.ones((0, 0))
            rhs = eps * np.prod([a[j] for j in T_a]) * (np.linalg.det(sub) if P else 1.0)
            assert abs(lhs - rhs) < 1e-9 * (1 + abs(lhs)), (S, lhs, rhs)


@lru_cache(maxsize=None)
def index(n: int) -> PluckerIndex:
    """The shared, immutable index of Gr(n, 2n)."""
    return PluckerIndex(n)
