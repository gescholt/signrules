"""
Independent replay of impossibility certificates.

A kill at the root comes with a certificate (``certificate.certificate_steps``): the backward slice of
the propagation log, a list of steps (basis, sign, why, premises) ending in a conflict.  This module
re-derives every step from the sign pattern and the signs established before it, using only the
identity or inequality that the step names, and then checks that the conflict is a genuine one.  It
shares with the search the bookkeeping between bases of Gr(n, 2n) and minors of B (``PluckerIndex``)
and nothing else: the propagation engine, the rule drivers and the search are not called.  A
certificate that replays is a proof of impossibility that depends only on the soundness of the eight
rules and the 3-term relations, each applied to explicit minors, not on the correctness of the code
that found it.

What is checked for each kind of step (the first token of ``why`` names the rule)::

    GP    a 3-term Grassmann–Plücker relation (the five other bases are the premises): two terms
          known and equal, the third opposite
    R1    det B > 0; a coefficient of the characteristic polynomial: all other principal k-minors
          negative, so this one is positive
    R2    Cramer: x*_j det B_JJ = Σ_k ± w_k det B_{J∖k, J∖j}, with w_k read off (a, B)
    R3    the pattern is symmetric and balanced; principal minors > 0, almost-principal minors have
          the sign of b_ij
    R4    a Z-twist (s, σ) of the Schur block (P, Q), all entries and right-hand sides known; the
          M-matrix or not-an-M-matrix conclusion about one minor of S D P_σ
    R5    a Laplace expansion of det B_PQ along some row or column
    R6    a minor of the Schur block whose permutation monomials all agree; a row of D x*_Q = A_Q
    R7    a circuit of the oriented matroid of [I | Bᵀ], one-sided against the tope (+, sign a)
    R8    a Z-twist of B itself with the right growth signs (a conflict, checked from the pattern)

Steps are replayed in an order in which every premise is known (the logs are not topologically
ordered); a step whose premises never become known fails.  Usage::

    from signrules import Ecology, complete, replay_kill
    report = replay_kill(Ecology(a, B))          # runs the root propagation and replays its certificate
    assert report.ok, report.message
"""
from __future__ import annotations

import ast
import functools
import re
from dataclasses import dataclass, field
from itertools import combinations, permutations
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .plucker_index import UNKNOWN, perm_sign, index

Index = Tuple[int, ...]
MinorKey = Tuple[Index, Index]
GP_RULE = "Grassmann–Plücker 3-term relation"


class ReplayError(Exception):
    """A step or the conflict of a certificate does not follow; the message says which and why."""


def _malformed_is_replay_error(method):
    """A step that names bases, indices or rules that do not fit the pattern is a failed replay, not a crash:
    ValueError, IndexError, KeyError and TypeError raised while checking it become ReplayError."""
    @functools.wraps(method)
    def wrapped(self, *args, **kw):
        try:
            return method(self, *args, **kw)
        except ReplayError:
            raise
        except (ValueError, IndexError, KeyError, TypeError) as e:
            raise ReplayError(f"malformed step ({type(e).__name__}: {e})") from e
    return wrapped


@dataclass
class ReplayReport:
    """Outcome of a replay: ``ok``, a message, the number of steps verified and the steps per rule."""
    ok: bool
    message: str
    steps: int = 0                                  # steps verified
    rules: Dict[str, int] = field(default_factory=dict)   # rule -> number of steps (incl. the conflict)
    conflict: str = ""

    def __bool__(self):
        return self.ok


def _merge_sign(P: Sequence[int], I: Sequence[int]) -> int:
    """Sign of the permutation sorting the concatenation of the sorted disjoint lists P and I."""
    inv = sum(1 for p in P for i in I if i < p)
    return -1 if inv % 2 else 1


def _tuple(s: str) -> Index:
    v = ast.literal_eval(s)
    return tuple(v) if isinstance(v, tuple) else (v,)


class Replayer:
    """The state of a replay: the sign pattern and the minors whose signs have been established."""

    def __init__(self, a_signs: Sequence[int], B_signs):
        self.a = [int(x) for x in a_signs]
        self.B = np.array(B_signs, dtype=int)
        self.n = n = len(self.a)
        if self.B.shape != (n, n) or any(self.B[i, i] != 1 for i in range(n)):
            raise ValueError("B must be n×n with +1 on the diagonal")
        self.idx = index(n)
        self.species = tuple(range(n))
        # basis S <-> minor (P, Q): sign det B_PQ = c(S) · chi[S]
        self.key_of: Dict[Index, MinorKey] = {}
        self.c_of: Dict[Index, int] = {}
        self.basis_of: Dict[MinorKey, Index] = {}
        for S in self.idx.bases:
            T_a, P, Q, eps = self.idx.decompose(S)
            c = eps
            for j in T_a:
                c *= self.a[j]
            self.key_of[S], self.c_of[S], self.basis_of[(P, Q)] = (P, Q), c, S
        self.minor: Dict[MinorKey, int] = {}        # established signs of minors of B
        self.minor[((), ())] = 1
        self._initial()

    # ------------------------------------------------------------------ the initial table

    def _initial(self):
        """Signs of the minors all of whose Leibniz terms agree (recomputed, not read from the search)."""
        B = self.B
        for (P, Q) in self.basis_of:
            k = len(P)
            if k == 0:
                continue
            signs = set()
            for perm in permutations(range(k)):
                t = perm_sign(perm)
                for r, pc in zip(P, perm):
                    t *= int(B[r, Q[pc]])
                signs.add(t)
                if len(signs) > 1:
                    break
            if len(signs) == 1:
                self.minor[(P, Q)] = signs.pop()

    # ------------------------------------------------------------------ helpers

    def sign(self, key: MinorKey):
        return self.minor.get((tuple(key[0]), tuple(key[1])), UNKNOWN)

    def chi(self, S: Index):
        m = self.sign(self.key_of[S])
        return UNKNOWN if m == UNKNOWN else self.c_of[S] * m

    def w(self, k: int, J: Sequence[int]):
        """Sign of w_k = a_k − Σ_{r∉J} b_kr x*_r, when the terms agree."""
        outside = [r for r in self.species if r not in J]
        return self.a[k] if all(-self.B[k, r] == self.a[k] for r in outside) else UNKNOWN

    def set(self, key: MinorKey, sign: int, what: str):
        key = (tuple(key[0]), tuple(key[1]))
        cur = self.minor.get(key, UNKNOWN)
        if cur == UNKNOWN:
            self.minor[key] = sign
        elif cur != sign:
            raise ReplayError(f"{what}: det B_{key[0]},{key[1]} already has the opposite sign")

    # a signed sum  Σ coef_t · sign det B_{key_t}  = target  (coef UNKNOWN: opaque term; key None: constant)
    def _term_signs(self, terms):
        out = []
        for coef, key in terms:
            if coef == UNKNOWN:
                out.append(UNKNOWN)
            elif key is None:
                out.append(coef)
            else:
                m = self.sign(key)
                out.append(UNKNOWN if m == UNKNOWN else coef * m)
        return out

    def check_sum_step(self, terms, target, target_key, key, sign, what: str):
        """The step forces sign det B_key = sign from the identity Σ terms = target."""
        key = (tuple(key[0]), tuple(key[1]))
        signs = self._term_signs(terms)
        if target_key is not None and (tuple(target_key[0]), tuple(target_key[1])) == key:
            # the right-hand side is forced: every term known and of one sign
            if any(s == UNKNOWN for s in signs) or len(set(signs)) != 1:
                raise ReplayError(f"{what}: terms do not determine the right-hand side")
            if signs[0] != sign:
                raise ReplayError(f"{what}: right-hand side would be {signs[0]}, not {sign}")
            return
        if target == UNKNOWN:
            if target_key is None:
                raise ReplayError(f"{what}: right-hand side unknown")
            target = self.sign(target_key)
            if target == UNKNOWN:
                raise ReplayError(f"{what}: right-hand side unknown")
        pos = [t for t, (coef, k) in enumerate(terms)
               if k is not None and (tuple(k[0]), tuple(k[1])) == key]
        if len(pos) != 1:
            raise ReplayError(f"{what}: the forced minor is not a term of the identity")
        t0 = pos[0]
        coef = terms[t0][0]
        if coef == UNKNOWN:
            raise ReplayError(f"{what}: the coefficient of the forced term is unknown")
        others = [s for t, s in enumerate(signs) if t != t0]
        if any(s == UNKNOWN for s in others) or any(s != -target for s in others):
            raise ReplayError(f"{what}: the other terms are not all known and of sign {-target}")
        if coef * target != sign:
            raise ReplayError(f"{what}: identity forces {coef * target}, step claims {sign}")

    def check_sum_conflict(self, terms, target, target_key, what: str):
        """The identity ``sum of terms = target`` is violated: every term is known and has the sign opposite to the target."""
        signs = self._term_signs(terms)
        if target == UNKNOWN:
            target = self.sign(target_key) if target_key is not None else UNKNOWN
        if target == UNKNOWN:
            raise ReplayError(f"{what}: right-hand side unknown")
        if any(s == UNKNOWN for s in signs) or any(s != -target for s in signs):
            raise ReplayError(f"{what}: not every term is known and of sign {-target}")

    # ------------------------------------------------------------------ the identities, rule by rule

    def r1_terms(self, k: int):
        """R1: the principal minors of order k, whose sum is the k-th coefficient of the characteristic polynomial."""
        return [(+1, (I, I)) for I in combinations(self.species, k)]

    def r2_terms(self, J: Index, j: int):
        """R2: Cramer's rule for species j in the subsystem J, as signed terms (sign of a_k times cofactor, minor)."""
        pj = J.index(j)
        cols = tuple(c for c in J if c != j)
        terms = []
        for pk, k in enumerate(J):
            rows = tuple(r for r in J if r != k)
            cof = -1 if (pj + pk) % 2 else 1
            wk = self.w(k, J)
            terms.append((UNKNOWN if wk == UNKNOWN else cof * wk, (rows, cols)))
        return terms

    def r5_expansions(self, P: Index, Q: Index):
        """R5: the Laplace expansions of det B_{P,Q} along each row and each column, one list of terms each."""
        for along_rows in (True, False):
            for pr, r in enumerate(P if along_rows else Q):
                terms = []
                for pc, c in enumerate(Q if along_rows else P):
                    i, j = (r, c) if along_rows else (c, r)
                    rows = tuple(x for x in P if x != i)
                    cols = tuple(x for x in Q if x != j)
                    coef = (-1 if (pr + pc) % 2 else 1) * int(self.B[i, j])
                    terms.append((coef, (rows, cols)))
                yield terms

    def balanced_signature(self) -> Optional[List[int]]:
        B = self.B
        if not np.array_equal(B, B.T):
            return None
        s = [1] + [-int(B[0, j]) for j in range(1, self.n)]
        for i, j in permutations(self.species, 2):
            if s[i] * s[j] * B[i, j] != -1:
                return None
        return s

    def schur_block(self, P: Index, Q: Index):
        """(dPP, D, A): signs of the Schur block D = (B/B_PP)_QQ and of A_Q in D x*_Q = A_Q."""
        dPP = self.sign((P, P)) if P else 1
        if dPP == UNKNOWN:
            raise ReplayError(f"det B_{P},{P} unknown")
        D, A = {}, {}
        for i in Q:
            for j in Q:
                m = self.sign((tuple(sorted(P + (i,))), tuple(sorted(P + (j,)))))
                D[(i, j)] = UNKNOWN if m == UNKNOWN else _merge_sign(P, (i,)) * _merge_sign(P, (j,)) * m * dPP
        for i in Q:
            rows = tuple(sorted(P + (i,)))
            ts = set()
            for pos, r in enumerate(rows):
                wr = self.w(r, P + Q)
                sub = self.sign((tuple(x for x in rows if x != r), P)) if P else 1
                ts.add(UNKNOWN if (wr == UNKNOWN or sub == UNKNOWN) else (-1 if (pos + len(P)) % 2 else 1) * wr * sub)
            A[i] = _merge_sign(P, (i,)) * ts.pop() * dPP if (len(ts) == 1 and UNKNOWN not in ts) else UNKNOWN
        return dPP, D, A

    @staticmethod
    def is_z_twist(Q: Index, D, s: Dict[int, int], sigma: Dict[int, int]) -> bool:
        """s_i d_{i,σ(i)} > 0 and s_i d_ij < 0 for j ≠ σ(i); σ a permutation of Q."""
        if sorted(sigma[i] for i in Q) != sorted(Q):
            return False
        for i in Q:
            for j in Q:
                d = D[(i, j)]
                if d == UNKNOWN or (s[i] * d > 0) != (j == sigma[i]):
                    return False
        return True

    def twisted_minor(self, P: Index, s, sigma, I: Index, J: Index) -> Tuple[int, MinorKey]:
        """det (S D P_σ)_{I,J} = coef · sign det B_key · det B_PP."""
        cols = [sigma[c] for c in J]
        coef = perm_sign(cols) * _merge_sign(P, I) * _merge_sign(P, sorted(cols))
        for i in I:
            coef *= s[i]
        return coef, (tuple(sorted(P + tuple(I))), tuple(sorted(P + tuple(cols))))

    def circuit_terms(self, keys: Sequence[MinorKey]):
        """The circuit of the oriented matroid of N = [I | Bᵀ] spanned by the minors ``keys`` = (rows, cols),
        and its terms (coef, (Q, P)): T(e_k) C(e_k) = coef · sign det B_QP.  The element set of N is
        {0..n-1} (the identity) ∪ {n..2n-1} (the columns of Bᵀ); det B_{Q,P} is the basis ([n]∖P) ∪ (n+Q)."""
        n = self.n
        nbases = [set(x for x in range(n) if x not in cols) | set(n + r for r in rows) for rows, cols in keys]
        C = sorted(set().union(*nbases))
        if len(C) != n + 1 or any(len(b) != n for b in nbases):
            raise ReplayError("R7: the minors do not form one circuit")
        T = [1] * n + list(self.a)
        terms = []
        for k, e in enumerate(C):
            S = tuple(x for x in C if x != e)
            T1 = tuple(x for x in S if x < n)
            Q = tuple(x - n for x in S if x >= n)
            P = tuple(x for x in range(n) if x not in T1)
            eps = -1 if (sum(T1) - len(T1) * (len(T1) - 1) // 2) % 2 else 1
            terms.append((T[e] * (-1 if k % 2 else 1) * eps, (Q, P)))
        return terms

    # ------------------------------------------------------------------ parsing the rule strings

    _R4 = re.compile(r"R4 Z-block P=(\(.*?\)) Q=(\(.*?\)) s=(\(.*?\)) σ=(\(.*?\))(?:: (.*))?$")
    _R6 = re.compile(r"R6 Schur block P=(\(.*?\)) Q=(\(.*?\))(?:: (.*))?$")
    _R2 = re.compile(r"R2 Cramer: x\*_(\d+) det B_JJ, J=(\(.*?\))$")
    _R5 = re.compile(r"R5 Laplace expansion of det B_(\(.*?\)),(\(.*?\))$")
    _R1C = re.compile(r"R1 coefficient of λ\^(\d+)")
    _R8 = re.compile(r"R8 real eigenvalue s=(\(.*?\)) σ=(\(.*?\)):")
    _ORDER = re.compile(r"ordering (\(.*?\))")

    # ------------------------------------------------------------------ one step

    @_malformed_is_replay_error
    def verify_step(self, S: Index, chi_sign: int, why: str, premises: Sequence[Index]):
        """Check that sign chi[S] = chi_sign follows by the rule ``why`` from the established signs, and record it."""
        key = self.key_of[S]
        sign = self.c_of[S] * chi_sign                       # the claimed sign of det B_key
        P, Q = key
        what = why.split(":")[0]
        if why == GP_RULE:
            self._verify_gp_step(S, chi_sign, premises)
        elif why == "R1 det B > 0":
            if key != (self.species, self.species) or sign != 1:
                raise ReplayError("R1: det B > 0 misapplied")
        elif self._R1C.match(why):
            k = self.n - int(self._R1C.match(why).group(1))
            if P != Q or len(P) != k:
                raise ReplayError("R1: not a principal k-minor")
            self.check_sum_step(self.r1_terms(k), +1, None, key, sign, what)
        elif self._R2.match(why):
            m = self._R2.match(why)
            j, J = int(m.group(1)), _tuple(m.group(2))
            self.check_sum_step(self.r2_terms(J, j), UNKNOWN, (J, J), key, sign, what)
        elif why.startswith("R3 balanced"):
            if self.balanced_signature() is None:
                raise ReplayError("R3: the pattern is not symmetric and balanced")
            if P == Q:
                if sign != 1:
                    raise ReplayError("R3: principal minor must be positive")
            else:
                I = tuple(sorted(set(P) & set(Q)))
                i, j = (set(P) - set(I)), (set(Q) - set(I))
                if len(i) != 1 or len(j) != 1:
                    raise ReplayError("R3: not an almost-principal minor")
                (i,), (j,) = i, j
                if sign != int(self.B[i, j]) * _merge_sign(I, (i,)) * _merge_sign(I, (j,)):
                    raise ReplayError("R3: wrong sign for the almost-principal minor")
        elif self._R4.match(why):
            self._verify_r4(why, key, sign, conflict=False)
        elif self._R5.match(why):
            m = self._R5.match(why)
            P5, Q5 = _tuple(m.group(1)), _tuple(m.group(2))
            errors = []
            for terms in self.r5_expansions(P5, Q5):
                try:
                    self.check_sum_step(terms, UNKNOWN, (P5, Q5), key, sign, what)
                    break
                except ReplayError as e:
                    errors.append(str(e))
            else:
                raise ReplayError("R5: no Laplace expansion forces the sign: " + errors[0])
        elif self._R6.match(why):
            self._verify_r6(why, key, sign)
        elif why == "R7 tope: circuit orthogonality":
            terms = self.circuit_terms([self.key_of[p] for p in premises] + [key])
            signs = self._term_signs(terms)
            pos = [t for t, (_, k) in enumerate(terms) if k == key]
            if len(pos) != 1:
                raise ReplayError("R7: the forced minor is not in the circuit")
            others = {s for t, s in enumerate(signs) if t != pos[0]}
            if UNKNOWN in others or len(others) != 1:
                raise ReplayError("R7: the other terms are not known and equal")
            (common,) = others
            if terms[pos[0]][0] * sign != -common:
                raise ReplayError("R7: orthogonality forces the other sign")
        else:
            raise ReplayError(f"unsupported step: {why}")
        self.set(key, sign, what)

    def _verify_gp_step(self, S: Index, chi_sign: int, premises: Sequence[Index]):
        """Two terms of the 3-term relation known and equal; the term containing S has the other sign."""
        bases = list(dict.fromkeys(list(premises) + [S]))      # the six bases of the relation
        if len(bases) != 6:
            raise ReplayError("GP: a step names the five other bases of its relation")
        C = set(bases[0]).intersection(*map(set, bases[1:]))
        rest = sorted(set().union(*map(set, bases)) - C)
        if len(C) != self.n - 2 or len(rest) != 4:
            raise ReplayError("GP: the premises do not form a 3-term relation")
        a, b, c, d = rest
        Cl = sorted(C)

        def term(u, v, w, z, coef):
            k1, k2 = sorted([u, v] + Cl), sorted([w, z] + Cl)
            return tuple(k1), tuple(k2), coef * perm_sign([u, v] + Cl) * perm_sign([w, z] + Cl)

        rel = [term(a, b, c, d, 1), term(a, c, b, d, -1), term(a, d, b, c, 1)]
        mine = [t for t in rel if S in t[:2]]
        if len(mine) != 1:
            raise ReplayError("GP: the forced basis is not in the relation")
        (k1, k2, coef), = mine
        partner = k2 if k1 == S else k1
        vp = self.chi(partner)
        if vp == UNKNOWN:
            raise ReplayError("GP: the partner of the forced basis is unknown")
        others = []
        for t in rel:
            if t is mine[0]:
                continue
            v1, v2 = self.chi(t[0]), self.chi(t[1])
            if v1 == UNKNOWN or v2 == UNKNOWN:
                raise ReplayError("GP: another term is unknown")
            others.append(t[2] * v1 * v2)
        if others[0] != others[1]:
            raise ReplayError("GP: the two known terms disagree, nothing is forced")
        if coef * chi_sign * vp != -others[0]:
            raise ReplayError("GP: the relation forces the other sign")

    def _parse_r4(self, why: str):
        m = self._R4.match(why)
        P, Q = _tuple(m.group(1)), _tuple(m.group(2))
        sv, sg = _tuple(m.group(3)), _tuple(m.group(4))
        s = {i: int(v) for i, v in zip(Q, sv)}
        sigma = {i: int(v) for i, v in zip(Q, sg)}
        clause = m.group(5) or ""
        dPP, D, A = self.schur_block(P, Q)
        if any(v == UNKNOWN for v in D.values()) or any(v == UNKNOWN for v in A.values()):
            raise ReplayError("R4: the Schur block is not fully known")
        if not self.is_z_twist(Q, D, s, sigma):
            raise ReplayError("R4: (s, σ) is not a Z-twist of the block")
        twisted = [s[i] * A[i] for i in Q]
        return P, Q, s, sigma, clause, dPP, twisted

    def _verify_r4(self, why: str, key: MinorKey, sign: int, conflict: bool):
        P, Q, s, sigma, clause, dPP, twisted = self._parse_r4(why)
        if clause.startswith("M-matrix"):
            if not all(t > 0 for t in twisted):
                raise ReplayError("R4: s·A_Q is not positive")
            R = tuple(x for x in key[0] if x not in P)
            Cp = tuple(x for x in key[1] if x not in P)
            inv = {v: k for k, v in sigma.items()}
            J = tuple(sorted(inv[c] for c in Cp))
            if "principal minor > 0" in clause and "almost" not in clause:
                if R != J:
                    raise ReplayError("R4: not a principal minor of S D P_σ")
                coef, k2 = self.twisted_minor(P, s, sigma, R, R)
                if k2 != key or sign != coef * dPP:
                    raise ReplayError("R4: wrong sign for the principal minor")
            elif "almost-principal" in clause:
                I = tuple(sorted(set(R) & set(J)))
                i, j = set(R) - set(I), set(J) - set(I)
                if len(i) != 1 or len(j) != 1:
                    raise ReplayError("R4: not an almost-principal minor")
                (i,), (j,) = i, j
                coef, k2 = self.twisted_minor(P, s, sigma, tuple(sorted(I + (i,))), tuple(sorted(I + (j,))))
                coef *= _merge_sign(I, (i,)) * _merge_sign(I, (j,))
                if k2 != key or sign != -coef * dPP:
                    raise ReplayError("R4: wrong sign for the almost-principal minor")
            else:
                raise ReplayError(f"R4: unknown clause {clause!r}")
        elif clause.startswith("not an M-matrix"):
            if not all(t < 0 for t in twisted):
                raise ReplayError("R4: s·A_Q is not negative")
            m = self._ORDER.search(clause)
            if not m:
                raise ReplayError("R4: no ordering in the clause")
            order = _tuple(m.group(1))
            leading = []
            for mm in range(1, len(Q) + 1):
                I = tuple(sorted(order[:mm]))
                coef, k2 = self.twisted_minor(P, s, sigma, I, I)
                v = self.sign(k2)
                leading.append((k2, coef, UNKNOWN if v == UNKNOWN else coef * v * dPP))
            if conflict:
                if not all(v == 1 for _, _, v in leading):
                    raise ReplayError("R4: not every leading minor of the ordering is known positive")
                return
            mine = [t for t in leading if t[0] == key]
            if len(mine) != 1:
                raise ReplayError("R4: the forced minor is not a leading minor of the ordering")
            if not all(v == 1 for t in leading if t is not mine[0] for v in [t[2]]):
                raise ReplayError("R4: the other leading minors are not known positive")
            if sign != -mine[0][1] * dPP:
                raise ReplayError("R4: the leading minor must be forced ≤ 0")
        else:
            raise ReplayError(f"R4: unknown clause {clause!r}")

    def _verify_r6(self, why: str, key: MinorKey, sign: int):
        m = self._R6.match(why)
        P, Q = _tuple(m.group(1)), _tuple(m.group(2))
        clause = m.group(3) or ""
        dPP, D, A = self.schur_block(P, Q)
        if clause.startswith("minor determined"):
            I = tuple(x for x in key[0] if x not in P)
            J = tuple(x for x in key[1] if x not in P)
            if len(I) != len(J) or not I or any(x not in Q for x in I + J):
                raise ReplayError("R6: not a minor of the Schur block")
            if any(D[(i, j)] == UNKNOWN for i in I for j in J):
                raise ReplayError("R6: an entry of the block is unknown")
            monomials = set()
            for perm in permutations(range(len(I))):
                t = perm_sign(perm)
                for i, pj in zip(I, perm):
                    t *= D[(i, J[pj])]
                monomials.add(t)
            if len(monomials) != 1:
                raise ReplayError("R6: the permutation monomials do not agree")
            if sign != _merge_sign(P, I) * _merge_sign(P, J) * dPP * monomials.pop():
                raise ReplayError("R6: wrong sign for the block minor")
        elif clause.startswith("row "):
            i = int(re.match(r"row (\d+) of D x = A", clause).group(1))
            self.check_sum_step(self._r6_row_terms(P, Q, i, dPP), A[i], None, key, sign, "R6 row")
        else:
            raise ReplayError(f"R6: unknown clause {clause!r}")

    def _r6_row_terms(self, P, Q, i, dPP):
        return [(_merge_sign(P, (i,)) * _merge_sign(P, (j,)) * dPP,
                 (tuple(sorted(P + (i,))), tuple(sorted(P + (j,))))) for j in Q]

    # ------------------------------------------------------------------ the conflict

    _SET_CONFLICT = re.compile(r"^(.*): det B_(\(.*?\)),(\(.*?\)) must be ([<>]) 0 but is ([<>]) 0$")
    _SUM_CONFLICT = re.compile(r"^(.*): every term is ([<>]) 0 but the sum must be ([<>]) 0$")

    @_malformed_is_replay_error
    def verify_conflict(self, reason: str, premises: Sequence) -> str:
        """Check that ``reason`` is a contradiction with the established signs; returns the rule name."""
        if reason == "3-term Grassmann–Plücker relation":
            bases = [p for p in premises if isinstance(p[0], int)]
            C = set(bases[0]).intersection(*map(set, bases[1:]))
            rest = sorted(set().union(*map(set, bases)) - C)
            if len(C) != self.n - 2 or len(rest) != 4:
                raise ReplayError("GP conflict: not a 3-term relation")
            a, b, c, d = rest
            Cl = sorted(C)
            vals = []
            for (u, v), (w, z), coef in (((a, b), (c, d), 1), ((a, c), (b, d), -1), ((a, d), (b, c), 1)):
                k1, k2 = tuple(sorted([u, v] + Cl)), tuple(sorted([w, z] + Cl))
                x1, x2 = self.chi(k1), self.chi(k2)
                if x1 == UNKNOWN or x2 == UNKNOWN:
                    raise ReplayError("GP conflict: a term is unknown")
                vals.append(coef * perm_sign([u, v] + Cl) * perm_sign([w, z] + Cl) * x1 * x2)
            if len(set(vals)) != 1:
                raise ReplayError("GP conflict: the three terms do not all agree")
            return "GP"
        m = self._SET_CONFLICT.match(reason)
        if m:
            why, P, Q = m.group(1), _tuple(m.group(2)), _tuple(m.group(3))
            must = 1 if m.group(4) == ">" else -1
            key = (P, Q)
            if self.sign(key) != -must:
                raise ReplayError("conflict: the table does not hold the opposite sign")
            # the derivation must force `must`; verify it as a step on a scratch copy
            S = self.basis_of[key]
            saved = dict(self.minor)
            del self.minor[key]
            try:
                self.verify_step(S, self.c_of[S] * must, why, [])
            finally:
                self.minor = saved
            return why.split()[0]
        m = self._SUM_CONFLICT.match(reason)
        if m:
            why = m.group(1)
            if self._R1C.match(why):
                k = self.n - int(self._R1C.match(why).group(1))
                self.check_sum_conflict(self.r1_terms(k), +1, None, "R1")
                return "R1"
            if self._R2.match(why):
                mm = self._R2.match(why)
                j, J = int(mm.group(1)), _tuple(mm.group(2))
                self.check_sum_conflict(self.r2_terms(J, j), UNKNOWN, (J, J), "R2")
                return "R2"
            if self._R5.match(why):
                mm = self._R5.match(why)
                P5, Q5 = _tuple(mm.group(1)), _tuple(mm.group(2))
                for terms in self.r5_expansions(P5, Q5):
                    try:
                        self.check_sum_conflict(terms, UNKNOWN, (P5, Q5), "R5")
                        return "R5"
                    except ReplayError:
                        pass
                raise ReplayError("R5 conflict: no expansion is one-sided")
            if self._R6.match(why):
                mm = self._R6.match(why)
                P, Q = _tuple(mm.group(1)), _tuple(mm.group(2))
                i = int(re.search(r"row (\d+) of D x = A", why).group(1))
                dPP, D, A = self.schur_block(P, Q)
                self.check_sum_conflict(self._r6_row_terms(P, Q, i, dPP), A[i], None, "R6")
                return "R6"
            raise ReplayError(f"unsupported sum conflict: {why}")
        if self._R4.match(reason) and "not an M-matrix, but ordering" in reason:
            self._verify_r4(reason, ((), ()), 0, conflict=True)
            return "R4"
        if reason.startswith("R7 tope: circuit"):
            keys = [self.key_of[tuple(p)] if isinstance(p[0], int) else (tuple(p[0]), tuple(p[1])) for p in premises]
            terms = self.circuit_terms(keys)
            signs = self._term_signs(terms)
            if UNKNOWN in signs or len(set(signs)) != 1:
                raise ReplayError("R7 conflict: the circuit is not one-sided")
            return "R7"
        m = self._R8.match(reason)
        if m:
            s = {i: int(v) for i, v in zip(self.species, _tuple(m.group(1)))}
            sigma = {i: int(v) for i, v in zip(self.species, _tuple(m.group(2)))}
            D = {(i, j): int(self.B[i, j]) for i in self.species for j in self.species}
            if not self.is_z_twist(self.species, D, s, sigma):
                raise ReplayError("R8: (s, σ) is not a Z-twist of B")
            moved = [i for i in self.species if sigma[i] != i]
            if not moved:
                raise ReplayError("R8: σ is the identity")
            wrong = [i for i in self.species if s[i] * self.a[i] < 0]
            if len(wrong) == 0:
                return "R8"
            if len(wrong) == 1:
                seen, cycles = set(), 0
                for i in moved:
                    if i in seen:
                        continue
                    cycles += 1
                    j = i
                    while j not in seen:
                        seen.add(j)
                        j = sigma[j]
                if cycles % 2 == 0:                       # det B > 0 is the standing hypothesis (R1)
                    return "R8"
                raise ReplayError("R8: one wrong growth sign needs an even number of cycles")
            raise ReplayError("R8: more than one wrong growth sign")
        raise ReplayError(f"unsupported conflict: {reason}")

    # ------------------------------------------------------------------ the whole certificate

    def replay(self, steps: Sequence[Tuple], conflict: str, conflict_premises: Sequence) -> ReplayReport:
        """Replay steps (basis, chi_sign, why, premise_bases) in a premise-respecting order, then the conflict."""
        pending = list(steps)
        rules: Dict[str, int] = {}
        done = 0
        try:
            while pending:
                progress = False
                rest = []
                for S, chi_sign, why, prem in pending:
                    if self.chi(S) != UNKNOWN:
                        if self.chi(S) != chi_sign:
                            raise ReplayError(f"step {why}: basis {S} already has the opposite sign")
                        progress = True
                        continue
                    if any(self.chi(p) == UNKNOWN for p in prem):
                        rest.append((S, chi_sign, why, prem))
                        continue
                    self.verify_step(S, chi_sign, why, prem)
                    name = "GP" if why == GP_RULE else why.split()[0]
                    rules[name] = rules.get(name, 0) + 1
                    done += 1
                    progress = True
                pending = rest
                if pending and not progress:
                    raise ReplayError(f"{len(pending)} steps have premises that are never established")
            name = self.verify_conflict(conflict, conflict_premises)
            rules[name] = rules.get(name, 0) + 1
        except ReplayError as e:
            return ReplayReport(False, str(e), done, rules, conflict)
        except (ValueError, TypeError) as e:          # a step that is not (basis, sign, why, premises)
            return ReplayReport(False, f"malformed step ({type(e).__name__}: {e})", done, rules, conflict)
        return ReplayReport(True, f"{done} steps and the conflict replayed", done, rules, conflict)


def replay_certificate(a_signs, B_signs, steps, conflict: str, conflict_premises) -> ReplayReport:
    """Replay ``steps`` (as returned by ``certificate_steps``, in basis form) and the final ``conflict`` with its
    ``conflict_premises`` for the pattern ``(a_signs, B_signs)``, re-deriving every step from the identity or
    inequality it names."""
    return Replayer(a_signs, B_signs).replay(steps, conflict, conflict_premises)


def replay_kill(eco, result=None, use=None) -> ReplayReport:
    """Run the root propagation of ``eco`` (unless ``result`` is given) and replay its certificate."""
    from .search import complete
    from .certificate import certificate_steps
    if result is None:
        result = complete(eco, max_nodes=0, max_completions=1, keep_derivations=False) if use is None else \
            complete(eco, use=use, max_nodes=0, max_completions=1, keep_derivations=False)
    if result.root_conflict is None:
        return ReplayReport(False, "no kill at the root to replay", 0, {}, "")
    steps, conflict_premises, _ = certificate_steps(eco, result, backward_slice=True)
    return replay_certificate(eco.a, eco.B, steps, result.root_conflict, conflict_premises)
