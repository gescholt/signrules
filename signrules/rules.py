"""
Sign rules on the ecological Grassmannian.

Setting
-------
A Lotka–Volterra ecology  x' = diag(x)(a − B x)  with a sign pattern
(sign a, sign B), b_ii > 0.  The pattern is *possible* if some choice of
magnitudes has a feasible equilibrium x* = B⁻¹a > 0 whose Jacobian
−diag(x*) B is Hurwitz, and *impossible* otherwise.

Everything here works at the level of signs of the minors of B — equivalently
(see plucker_index.py) of the Plücker coordinates of [diag(a) | B] in Gr(n, 2n).
A *sign table* ``chi`` maps each basis of Gr(n, 2n) to +1, −1 or UNKNOWN.

Each rule inspects a partial sign table and either

* raises ``Conflict``   – the pattern is impossible (with a reason), or
* forces some signs     – they are consequences of the standing hypotheses, or
* does nothing.

Standing hypotheses (used freely by the rules)
----------------------------------------------
* feasibility:   B x* = a with x* > 0;
* stability:     −diag(x*) B is Hurwitz, hence det B > 0 and all coefficients
                 of the characteristic polynomial of diag(x*) B are positive;
* genericity:    no minor of B vanishes (harmless: feasible + Hurwitz is an open
                 condition, so a possible pattern has a generic realisation).

Notation used below
-------------------
* P, Q ⊆ [n] index sets (always sorted tuples);  B_{P,Q} the submatrix.
* w_J = a_J − B_{J,J^c} x*_{J^c}  is the right-hand side of the rows J of
  B x* = a after moving the other species across.  Its sign is known when
  a_k and all −b_kr (r ∉ J) agree — that is the only way x* enters.
* Schur block: for disjoint P, Q,  D = (B/B_PP)_{QQ} = B_QQ − B_QP B_PP⁻¹ B_PQ
  satisfies  D x*_Q = A_Q  with  A_Q = w_Q − B_QP B_PP⁻¹ w_P.  By the quotient
  formula,  det D_{I,J} = det B_{(P,I),(P,J)} / det B_PP,  so entries and minors
  of D are ratios of minors of B and live in the same sign table.

The rules
---------
R1  stability      det B > 0; for each k some principal k-minor is positive.
R2  cramer         B_JJ x*_J = w_J  ⇒  x*_j det B_JJ = Σ_k cof_kj(B_JJ) w_k  has
                   the sign of det B_JJ  (J = [n]: the usual feasibility
                   inequalities;  "R2n" restricts to J = [n]).
R3  balanced       if S B S is a Z-matrix for a signature S, then diag(x*) S B S
                   is an M-matrix, so all principal minors of B are positive and
                   every almost-principal minor det B_{(I,i),(I,j)} has the sign
                   of b_ij  (Schur complements of M-matrices are M-matrices).
R4  z_block        if a row-sign twist s and a column permutation σ make the
                   Schur block D a Z-matrix with positive diagonal, then
                   s·A_Q > 0  ⇒  s D P_σ is an M-matrix (principal minors > 0,
                   almost-principal minors < 0), while
                   s·A_Q < 0  ⇒  s D P_σ is not an M-matrix, so no ordering of Q
                   has all leading principal minors positive.
R5  laplace        expansion of every minor of B along every row and column.
R6  schur          (a) sign of a minor of a Schur block when all permutation
                   monomials in its entries agree;  (b) each row of D x*_Q = A_Q
                   with x*_Q > 0: some term d_ij x*_j has the sign of A_i.
R7  tope           the complete sign-level content of feasibility.  The sign
                   vectors (sign y, sign By), y ∈ R^n, are the covectors of the
                   oriented matroid of N = [I | Bᵀ], whose chirotope is
                   χ_N(T ∪ (n+Q)) = ε(T) · sign det B_{Q,P}, P = [n] ∖ T.  So
                   x* > 0 with sign(Bx*) = sign a says that (+ⁿ, sign a) is a
                   tope, i.e. orthogonal to every circuit; in the uniform case
                   every (n+1)-subset C = {e_0 < … < e_n} is a circuit with
                   C(e_k) = (−1)^k χ_N(C ∖ e_k), and orthogonality means the
                   n+1 products T(e_k) C(e_k) are not all equal.  This is the
                   oriented-matroid Farkas lemma; R2 is a special case.
R8  real_eig       if a twist s and a permutation σ ≠ id make S B P_σ a Z-matrix with
                   s·a > 0 (R4 on the whole matrix, M-matrix case), or with
                   s_i a_i < 0 for exactly one species i and σ having an even
                   number of cycles, then −diag(x*) B has a positive real
                   eigenvalue: the pattern is impossible.  (D_t = S B X P_σ + t S P_σ
                   stays a Z-matrix for t ≥ 0; its least real eigenvalue is > 0 at
                   t = 0 and ≤ 0 for large t, so D_{t₁} = S (B X + t₁ I) P_σ is
                   singular for some t₁ > 0 and −t₁ is an eigenvalue of B X.  With
                   one wrong sign D_0 is still an M-matrix when k is even.)

R7 and R8 are not in the papers; R1–R6 correspond to Çelik–Isoldi–Portakal–Zucal, "Stable Coexistence in
Ecologies and Games" (arXiv:2609.09286) as follows: R3 ↔ Theorems 3.1–3.2, R4 ↔ Theorems 3.4–3.6 and the
Schur steps of Lemmas 3.8–3.10, R2/R5/R6 ↔ the Cramer, Laplace and "forced
sign" steps inside those proofs.

R4 and R6 are subsumed (docs/SIGNRULES.md, "R4 and R6 are subsumed"): R6 is a special case of
the 3-term relations (part a, which can only fire on 2×2 minors of D) and of R7 (part b, one
circuit); R4's conclusions are reached by R2, R3, R5, R7 and R8 — proved for the whole matrix
and for M-matrix blocks with P = ∅; for the rest R4 is stronger than one round of the other rules
(at n = 5 it fixes a 3×3 leading minor of Z they leave open) but redundant up to one branching: the
opposite sign conflicts at their fixpoint, checked exhaustively from R4's premises at n ≤ 5 and up
to symmetry at n = 6 (scripts/r4_subsumption.py), and the search never fires it (censuses).
``search.complete`` therefore idles them by default (their propagators stay registered for the
branching order); they remain here for reduced rule sets, the ladders, the certificates that cite
them and the proofs of that paper.
"""

from functools import lru_cache
from itertools import combinations, permutations, product
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .plucker_index import PluckerIndex, UNKNOWN, perm_sign, index, signed_perms

Index = Tuple[int, ...]
MinorKey = Tuple[Index, Index]          # (rows, cols), both sorted
IdTerm = Tuple[int, Optional[int], int]  # (coef, basis id or None, c): term sign = coef · c · chi[id]


class Conflict(Exception):
    """
    Raised when a rule derives a contradiction: the pattern is impossible.
    ``premises`` lists the minors (as (P, Q) keys) whose signs were used.
    """

    def __init__(self, message: str, premises=()):
        super().__init__(message)
        # premises are (P, Q) minor keys from the rules, or bases (tuples of ints) from the
        # 3-term propagation; keep bases as they are, normalise keys to tuples of tuples.
        self.premises = [k if (isinstance(k, tuple) and k and isinstance(k[0], int)) else tuple(map(tuple, k))
                         for k in premises]


def merge_sign(P: Iterable[int], I: Iterable[int]) -> int:
    """Sign of the permutation that sorts the concatenation (P, I) of two sorted, disjoint lists."""
    return _merge_sign(tuple(P), tuple(I))


@lru_cache(maxsize=None)
def _merge_sign(P: Tuple[int, ...], I: Tuple[int, ...]) -> int:
    inversions = sum(1 for p in P for i in I if i < p)
    return -1 if inversions % 2 else 1


class Ecology:
    """
    A sign pattern (a, B) together with the machinery to read and write signs
    of minors of B in a Plücker sign table.

    Parameters
    ----------
    a_signs : sequence of ±1 (signs of the growth rates)
    B_signs : n×n array of ±1, with +1 on the diagonal
    """

    def __init__(self, a_signs: Sequence[int], B_signs, name: str = ""):
        self.a = [int(s) for s in a_signs]
        self.B = np.array(B_signs, dtype=int)
        self.n = len(self.a)
        self.name = name
        if self.B.shape != (self.n, self.n):
            raise ValueError("B must be n×n")
        if any(self.B[i, i] != 1 for i in range(self.n)):
            raise ValueError("b_ii > 0 is assumed")
        self.species: Index = tuple(range(self.n))
        self.idx = index(self.n)                 # shared per n: bases, ids, decompositions

        # (P, Q) ↦ (basis id S, c) with  sign det B_{P,Q} = c · chi[S].
        self._minor_key: Dict[MinorKey, Tuple[int, int]] = {}
        for S, (_, P, Q, eps) in enumerate(self.idx.decomps):
            self._minor_key[(P, Q)] = (S, eps * self.idx.a_factor(self.a, P))
        self._trail: Optional[list] = None   # set_minor appends each newly set basis id here (engine.py)
        self._rhs_cache: Dict[Index, Dict[int, int]] = {}
        self._block_cache: Dict[MinorKey, tuple] = {}
        self._schur_cache: Dict[MinorKey, tuple] = {}
        self._laplace_expansions = self._build_laplace_expansions()
        self._circuits = self._build_circuits()

    # ------------------------------------------------------------------ tables

    def initial_table(self) -> List[int]:
        """
        Signs of the minors whose permutation expansion has monomials of one sign only, as a
        flat table indexed by basis id (+1, -1, UNKNOWN = 0; ``self.idx.to_dict`` names the bases).
        """
        return self.idx.initial_table(self.a, self.B)

    def minor(self, chi, P: Iterable[int], Q: Iterable[int]):
        """Sign of det B_{P,Q} in the table (rows and columns sorted), or UNKNOWN."""
        S, c = self._minor_key[(tuple(P), tuple(Q))]
        return c * chi[S]

    def set_minor(self, chi, P: Iterable[int], Q: Iterable[int], sign: int,
                  log: Optional[list], why: str, premises=()) -> bool:
        """
        Impose sign det B_{P,Q} = sign.  Returns True if the sign was new,
        False if it was already known, and raises Conflict if it contradicts
        the table.  Every forced sign is appended to ``log`` as
        (why, P, Q, sign, premises), where premises are the (P', Q') minors
        whose signs the rule used — a kill therefore comes with its certificate
        (see certificate.py).
        """
        P, Q = tuple(P), tuple(Q)
        S, c = self._minor_key[(P, Q)]
        current = chi[S]
        if current == UNKNOWN:
            chi[S] = c * sign
            if self._trail is not None:
                self._trail.append(S)
            if log is not None:
                log.append((why, P, Q, sign, [tuple(map(tuple, k)) for k in premises]))
            return True
        if current != c * sign:
            raise Conflict(f"{why}: det B_{P},{Q} must be {_pm(sign)} but is {_pm(-sign)}",
                           list(premises) + [(P, Q)])
        return False

    def rhs_signs(self, chi, J: Iterable[int]) -> Dict[int, int]:
        """
        Signs of w_k = a_k − Σ_{r∉J} b_kr x*_r for k ∈ J (UNKNOWN when the terms disagree).
        Depends on (a, B) only; cached per J (do not modify the returned dict).
        """
        J = tuple(J)
        w = self._rhs_cache.get(J)
        if w is None:
            outside = [r for r in self.species if r not in J]
            w = {k: self.a[k] if all(-self.B[k, r] == self.a[k] for r in outside) else UNKNOWN
                 for k in J}
            self._rhs_cache[J] = w
        return w

    def compile_terms(self, terms) -> Tuple[List[IdTerm], List[Optional[MinorKey]]]:
        """The terms (coef, key) of a signed sum as (coef, id, c) triples, with the keys alongside."""
        mk = self._minor_key
        iterms, keys = [], []
        for coef, key in terms:
            if key is None:
                iterms.append((coef, None, 1))
            else:
                S, c = mk[key]
                iterms.append((coef, S, c))
            keys.append(key)
        return iterms, keys

    def enforce_sum_ids(self, chi, iterms: List[IdTerm], keys, target, target_key: Optional[MinorKey],
                        log, why: str, extra_premises=(), scale: int = 1) -> bool:
        """
        ``enforce_sum`` on compiled terms (``compile_terms``): the sign of term k is
        scale · coef · c · chi[id], an UNKNOWN coefficient (0) makes the term opaque, and
        ``keys[k]`` is its minor for the premises.  Same consequences, conflicts and logs.
        """
        signs = []
        n_open = 0
        for coef, S, c in iterms:
            if S is None:
                s = coef * scale
            else:
                s = coef * c * chi[S] * scale
            if s == 0:
                n_open += 1
            signs.append(s)

        def used():
            out = list(extra_premises)
            out += [k for k, s in zip(keys, signs) if k is not None and s != 0]
            if target != UNKNOWN and target_key is not None:
                out.append(target_key)
            return out

        if target == UNKNOWN:
            if target_key is not None and n_open == 0 and len(set(signs)) == 1:
                return self.set_minor(chi, *target_key, signs[0], log, why, used())
            return False
        if n_open == 0:
            if all(s == -target for s in signs):
                raise Conflict(f"{why}: every term is {_pm(-target)} but the sum must be {_pm(target)}", used())
            return False
        if n_open == 1:
            k = signs.index(0)
            if all(s == -target for i, s in enumerate(signs) if i != k):
                coef, S, c = iterms[k]
                if coef != UNKNOWN and S is not None:
                    return self.set_minor(chi, *keys[k], coef * scale * target, log, why, used())
        return False

    # ------------------------------------------------ signed-sum propagation

    def enforce_sum(self, chi, terms: List[Tuple[object, Optional[MinorKey]]], target,
                    target_key: Optional[MinorKey], log, why: str, extra_premises=()) -> bool:
        """
        Propagate the identity   Σ_terms  =  target   at the level of signs.

        Each term is (coef, key): its sign is coef · sign det B_key, or coef alone
        when key is None; coef may itself be UNKNOWN (then the term is opaque).
        ``target`` is the sign of the right-hand side, or UNKNOWN; when it is
        UNKNOWN and ``target_key`` is given, the right-hand side may be forced.

        Consequences used:
        * target known, every term known and of sign −target      → Conflict;
        * target known, exactly one term open, the rest −target    → that term
          has sign target (forces its minor if its coefficient is known);
        * target unknown, every term known and of one sign         → target gets it.
        """
        signs = []
        for coef, key in terms:
            if coef == UNKNOWN:
                signs.append(UNKNOWN)
            elif key is None:
                signs.append(coef)
            else:
                m = self.minor(chi, *key)
                signs.append(UNKNOWN if m == UNKNOWN else coef * m)
        open_terms = [i for i, s in enumerate(signs) if s == UNKNOWN]
        used = list(extra_premises)
        used += [key for (coef, key), s in zip(terms, signs) if key is not None and s != UNKNOWN]
        if target != UNKNOWN and target_key is not None:
            used.append(target_key)

        if target == UNKNOWN:
            if target_key is not None and not open_terms and len(set(signs)) == 1:
                return self.set_minor(chi, *target_key, signs[0], log, why, used)
            return False
        if not open_terms:
            if all(s == -target for s in signs):
                raise Conflict(f"{why}: every term is {_pm(-target)} but the sum must be {_pm(target)}", used)
            return False
        if len(open_terms) == 1 and all(s == -target for i, s in enumerate(signs) if i != open_terms[0]):
            coef, key = terms[open_terms[0]]
            if coef != UNKNOWN and key is not None:
                return self.set_minor(chi, *key, coef * target, log, why, used)
        return False

    # ------------------------------------------------------------- R1

    def rule_stability(self, chi, log) -> bool:
        """det B > 0, and for each k some principal k-minor is positive (coefficient positivity)."""
        forced = self.set_minor(chi, self.species, self.species, +1, log, "R1 det B > 0")
        for terms, target, target_key, why in self.stability_sums():
            forced |= self.enforce_sum(chi, terms, target, target_key, log, why)
        return forced

    def stability_sums(self):
        """R1 as signed sums (terms, target, target_key, why): Σ principal k-minors > 0."""
        for k in range(2, self.n):
            terms = [(+1, (I, I)) for I in combinations(self.species, k)]
            yield terms, +1, None, f"R1 coefficient of λ^{self.n - k} (sum of principal {k}-minors) > 0"

    # ------------------------------------------------------------- R2

    def rule_cramer(self, chi, log, full_system_only: bool = False) -> bool:
        """
        Cramer's rule on the subsystem B_JJ x*_J = w_J:
            x*_j · det B_JJ = Σ_{k∈J} (−1)^{pos(j)+pos(k)} w_k · det B_{J∖k, J∖j} ,
        and the left-hand side has the sign of det B_JJ because x*_j > 0.
        """
        forced = False
        for terms, _, target_key, why in self.cramer_sums(full_system_only):
            forced |= self.enforce_sum(chi, terms, self.minor(chi, *target_key), target_key, log, why)
        return forced

    def cramer_sums(self, full_system_only: bool = False, sizes=None):
        """R2 as signed sums (terms, None, (J, J), why); the target is the sign of det B_JJ."""
        if sizes is None:
            sizes = [self.n] if full_system_only else range(1, self.n + 1)
        for size in sizes:
            for J in combinations(self.species, size):
                w = self.rhs_signs(None, J)          # depends on (a, B) only
                for pj, j in enumerate(J):
                    cols = tuple(c for c in J if c != j)
                    terms = []
                    for pk, k in enumerate(J):
                        rows = tuple(r for r in J if r != k)
                        cof = -1 if (pj + pk) % 2 else 1
                        coef = UNKNOWN if w[k] == UNKNOWN else cof * w[k]
                        terms.append((coef, (rows, cols)))
                    yield terms, None, (J, J), f"R2 Cramer: x*_{j} det B_JJ, J={J}"

    # ------------------------------------------------------------- R3

    def balancing_signature(self) -> Optional[List[int]]:
        """s ∈ {±1}^n with s_i s_j b_ij < 0 for all i ≠ j, if the pattern is symmetric and balanced."""
        if not np.array_equal(self.B, self.B.T):
            return None
        s = [1] + [-int(self.B[0, j]) for j in range(1, self.n)]
        for i, j in permutations(self.species, 2):
            if s[i] * s[j] * self.B[i, j] != -1:
                return None
        return s

    def rule_balanced(self, chi, log) -> bool:
        """
        Balanced pattern: S B S is a Z-matrix, diag(x*) S B S is Hurwitz, hence an
        M-matrix.  All principal minors of B are positive; the almost-principal
        minor det B_{(I,i),(I,j)} equals (Π_I s)² s_i s_j · det(SBS)_{(I,i),(I,j)}
        ≤ 0 · s_i s_j, i.e. has the sign of b_ij (up to the row/column sorting parity).
        """
        if self.balancing_signature() is None:
            return False
        forced = False
        for k in range(self.n + 1):
            for I in combinations(self.species, k):
                forced |= self.set_minor(chi, I, I, +1, log, "R3 balanced: principal minor > 0")
                others = [i for i in self.species if i not in I]
                for i, j in permutations(others, 2):
                    sign = int(self.B[i, j]) * merge_sign(I, (i,)) * merge_sign(I, (j,))
                    forced |= self.set_minor(chi, sorted(I + (i,)), sorted(I + (j,)), sign, log,
                                             "R3 balanced: almost-principal minor")
        return forced

    # ---------------------------------------------------- Schur blocks (R4, R6)

    def schur_block(self, chi, P: Index, Q: Index):
        """
        Signs of the Schur block D = (B/B_PP)_{QQ} and of the reduced right-hand
        side A_Q in D x*_Q = A_Q, as dicts D[(i,j)], A[i] with UNKNOWN entries
        allowed; returns None when det B_PP is unknown.

        Entries:  d_ij = det B_{(P,i),(P,j)} / det B_PP  (rows/cols in that order).
        RHS:      det B_PP · A_i = det [[B_PP, w_P], [B_iP, w_i]], expanded along
                  the w-column, with w = a − B_{·,R} x*_R, R = [n] ∖ (P ∪ Q).
        """
        dPP_key, entries, rhs = self._block_data(P, Q)
        if dPP_key is None:
            dPP = 1
            used = []
        else:
            S, c = dPP_key
            dPP = c * chi[S]
            if dPP == UNKNOWN:
                return None
            used = [(P, P)]
        D = {}
        umask = 0                                   # bit e set: entry e (row-major in Q × Q) is unknown
        for i, j, S, c, key, bit in entries:
            v = c * chi[S]
            if v:
                D[(i, j)] = v * dPP
                used.append(key)
            else:
                D[(i, j)] = UNKNOWN
                umask |= bit
        self._block_unknown_mask = umask
        A = {}
        for i, (mi, terms) in zip(Q, rhs):
            if terms is None:                       # some w_r is unknown: A_i is never determined
                A[i] = UNKNOWN
                continue
            term_signs = set()
            for par, S, c, key in terms:
                term_signs.add(par * c * chi[S] if S is not None else par)
            if len(term_signs) == 1 and UNKNOWN not in term_signs:
                A[i] = mi * term_signs.pop() * dPP
                used.extend(key for _, _, _, key in terms if key is not None)
            else:
                A[i] = UNKNOWN
        self._block_premises = used
        return dPP, D, A

    def _block_data(self, P: Index, Q: Index):
        """
        The pattern-level data behind ``schur_block(·, P, Q)``, computed once per block:
        the (id, c) of det B_PP (None when P is empty); for every entry (i, j) of D its
        (i, j, id, c, key) with the merge signs folded into c; and for every i in Q the
        merge sign and the terms (parity · w_r, id, c, key) of the right-hand side, or None
        when some w_r is unknown.
        """
        data = self._block_cache.get((P, Q))
        if data is not None:
            return data
        mk = self._minor_key
        dPP_key = mk[(P, P)] if P else None
        entries = []
        for qi, i in enumerate(Q):
            mi = merge_sign(P, (i,))
            for qj, j in enumerate(Q):
                key = (tuple(sorted(P + (i,))), tuple(sorted(P + (j,))))
                S, c = mk[key]
                entries.append((i, j, S, mi * merge_sign(P, (j,)) * c, key, 1 << (qi * len(Q) + qj)))
        w = self.rhs_signs(None, P + Q)
        rhs = []
        for i in Q:
            rows = tuple(sorted(P + (i,)))
            terms = []
            for pos, r in enumerate(rows):
                if w[r] == UNKNOWN:
                    terms = None
                    break
                par = (-1 if (pos + len(P)) % 2 else 1) * w[r]
                if P:
                    key = (tuple(x for x in rows if x != r), P)
                    S, c = mk[key]
                    terms.append((par, S, c, key))
                else:
                    terms.append((par, None, 1, None))
            rhs.append((merge_sign(P, (i,)), terms))
        data = (dPP_key, entries, rhs)
        self._block_cache[(P, Q)] = data
        return data

    def block_minor(self, chi, P: Index, s: Dict[int, int], sigma: Dict[int, int],
                    I: Index, J: Index) -> Tuple[int, MinorKey]:
        """
        The minor det (S D P_σ)_{I,J} of the twisted Schur block, as
        (coef, key) with  sign = coef · sign det B_key · sign det B_PP.
        """
        cols = [sigma[c] for c in J]
        coef = perm_sign(cols) * merge_sign(P, I) * merge_sign(P, sorted(cols))
        for i in I:
            coef *= s[i]
        return coef, (tuple(sorted(P + I)), tuple(sorted(P + tuple(cols))))

    @staticmethod
    def z_twists(Q: Index, D: Dict[Tuple[int, int], int]):
        """
        All (s, σ) such that s_i d_{i,σ(i)} > 0 and s_i d_ij < 0 for j ≠ σ(i):
        each row of D must have exactly one entry of one sign, and these
        exceptional positions must form a permutation of Q.
        """
        per_row = []
        for i in Q:
            options = []
            for s in (1, -1):
                positive = [j for j in Q if s * D[(i, j)] > 0]
                if len(positive) == 1:
                    options.append((s, positive[0]))
            if not options:
                return []
            per_row.append(options)
        twists = []
        for choice in product(*per_row):
            sigma = {i: c[1] for i, c in zip(Q, choice)}
            if len(set(sigma.values())) == len(Q):
                twists.append(({i: c[0] for i, c in zip(Q, choice)}, sigma))
        return twists

    def block_pairs(self, min_P: int):
        """All disjoint (P, Q) with |P| ≥ min_P and |Q| ≥ 2."""
        for kP in range(min_P, self.n - 1):
            for P in combinations(self.species, kP):
                rest = [i for i in self.species if i not in P]
                for kQ in range(2, len(rest) + 1):
                    for Q in combinations(rest, kQ):
                        yield P, Q

    def block_bases(self, P: Index, Q: Index) -> set:
        """
        Every basis that R4 or R6 on the Schur block (P, Q) reads or writes: det B_PP, the
        minors det B_{P∪I, P∪J} (I, J ⊆ Q) and the minors in the reduced right-hand side.
        """
        keys = [(P, P)] if P else []
        for k in range(1, len(Q) + 1):
            for I in combinations(Q, k):
                for J in combinations(Q, k):
                    keys.append((tuple(sorted(P + I)), tuple(sorted(P + J))))
        if P:
            for i in Q:
                rows = tuple(sorted(P + (i,)))
                keys += [(tuple(x for x in rows if x != r), P) for r in rows]
        return {self._minor_key[k][0] for k in keys}

    # ------------------------------------------------------------- R4

    def rule_z_block(self, chi, log) -> bool:
        """Z-semipositivity on twisted Schur blocks (see module docstring, R4)."""
        forced = False
        for P, Q in self.block_pairs(min_P=0):
            forced |= self.z_block(chi, log, P, Q)
        return forced

    def z_block(self, chi, log, P: Index, Q: Index) -> bool:
        """R4 on the single Schur block (P, Q)."""
        data = self.schur_block(chi, P, Q)
        if data is None:
            return False
        dPP, D, A = data
        if self._block_unknown_mask or any(v == UNKNOWN for v in A.values()):
            return False
        premises = list(self._block_premises)
        forced = False
        for s, sigma in self.z_twists(Q, D):
            twisted_rhs = [s[i] * A[i] for i in Q]
            why = (f"R4 Z-block P={P} Q={Q} s={tuple(s[i] for i in Q)} "
                   f"σ={tuple(sigma[i] for i in Q)}")
            if all(t > 0 for t in twisted_rhs):
                forced |= self._m_matrix_signs(chi, P, Q, dPP, s, sigma, log, why, premises)
            elif all(t < 0 for t in twisted_rhs):
                forced |= self._not_m_matrix_signs(chi, P, Q, dPP, s, sigma, log, why, premises)
        return forced

    def _m_matrix_signs(self, chi, P, Q, dPP, s, sigma, log, why, premises) -> bool:
        """S D P_σ is an M-matrix: principal minors > 0, almost-principal minors < 0."""
        forced = False
        for k in range(1, len(Q) + 1):
            for I in combinations(Q, k):
                coef, key = self.block_minor(chi, P, s, sigma, I, I)
                forced |= self.set_minor(chi, *key, coef * dPP, log, why + ": M-matrix, principal minor > 0", premises)
        for k in range(len(Q) - 1):
            for I in combinations(Q, k):
                others = [q for q in Q if q not in I]
                for i, j in permutations(others, 2):
                    rows, cols = tuple(sorted(I + (i,))), tuple(sorted(I + (j,)))
                    coef, key = self.block_minor(chi, P, s, sigma, rows, cols)
                    coef *= merge_sign(I, (i,)) * merge_sign(I, (j,))   # rows/cols in the order (I,i),(I,j)
                    forced |= self.set_minor(chi, *key, -coef * dPP, log,
                                             why + ": M-matrix, almost-principal minor < 0", premises)
        return forced

    def _not_m_matrix_signs(self, chi, P, Q, dPP, s, sigma, log, why, premises) -> bool:
        """S D P_σ is a Z-matrix but not an M-matrix: no ordering has all leading minors positive."""
        forced = False
        for order in permutations(Q):
            leading = []
            for m in range(1, len(Q) + 1):
                I = tuple(sorted(order[:m]))
                coef, key = self.block_minor(chi, P, s, sigma, I, I)
                value = self.minor(chi, *key)
                leading.append((key, coef, UNKNOWN if value == UNKNOWN else coef * value * dPP))
            values = [v for _, _, v in leading]
            known_keys = [k for k, _, v in leading if v != UNKNOWN]
            if all(v == 1 for v in values):
                raise Conflict(why + f": not an M-matrix, but ordering {order} has all leading minors > 0",
                               premises + known_keys)
            unknown = [t for t in leading if t[2] == UNKNOWN]
            if len(unknown) == 1 and all(v == 1 for v in values if v != UNKNOWN):
                key, coef, _ = unknown[0]
                forced |= self.set_minor(chi, *key, -coef * dPP, log,
                                         why + f": not an M-matrix, leading minor of ordering {order} ≤ 0",
                                         premises + known_keys)
        return forced

    # ------------------------------------------------------------- R6

    def rule_schur(self, chi, log) -> bool:
        """Monomial-determined minors of Schur blocks, and the rows of D x*_Q = A_Q."""
        forced = False
        for P, Q in self.block_pairs(min_P=1):
            forced |= self.schur(chi, log, P, Q)
        return forced

    def schur(self, chi, log, P: Index, Q: Index) -> bool:
        """R6 on the single Schur block (P, Q)."""
        forced = False
        data = self.schur_block(chi, P, Q)
        if data is not None:
            dPP, D, A = data
            why, minors, rows = self._schur_data(P, Q)
            block_premises = list(self._block_premises)
            umask = self._block_unknown_mask
            # (a) det D_{I,J} when every permutation monomial in its entries has the same sign
            for k, I, J, mask, coef, rowsP, colsP, premises in minors:
                if umask & mask:                    # some entry of D_{I,J} is unknown
                    continue
                monomials = set()
                for t, perm in signed_perms(k):
                    for i, pj in zip(I, perm):
                        t *= D[(i, J[pj])]
                    monomials.add(t)
                    if len(monomials) == 2:
                        break
                if len(monomials) == 1:
                    forced |= self.set_minor(chi, rowsP, colsP, coef * dPP * monomials.pop(), log,
                                             why + ": minor determined by its entries", premises)
            # (b) row i of D x*_Q = A_Q
            for i, iterms, keys, why_row in rows:
                forced |= self.enforce_sum_ids(chi, iterms, keys, A[i], None, log, why_row,
                                               extra_premises=block_premises, scale=dPP)
        return forced

    def _schur_data(self, P: Index, Q: Index):
        """The pattern-level data of R6 on the block (P, Q): rule name, the minors of D, the rows."""
        data = self._schur_cache.get((P, Q))
        if data is not None:
            return data
        why = f"R6 Schur block P={P} Q={Q}"
        pos = {q: k for k, q in enumerate(Q)}
        minors = []
        for k in range(2, len(Q) + 1):
            for I in combinations(Q, k):
                for J in combinations(Q, k):
                    entry_keys = [(tuple(sorted(P + (i,))), tuple(sorted(P + (j,)))) for i in I for j in J]
                    mask = sum(1 << (pos[i] * len(Q) + pos[j]) for i in I for j in J)
                    minors.append((k, I, J, mask, merge_sign(P, I) * merge_sign(P, J),
                                   tuple(sorted(P + I)), tuple(sorted(P + J)), [(P, P)] + entry_keys))
        rows = []
        for i in Q:
            terms = [(merge_sign(P, (i,)) * merge_sign(P, (j,)),
                      (tuple(sorted(P + (i,))), tuple(sorted(P + (j,))))) for j in Q]
            iterms, keys = self.compile_terms(terms)
            rows.append((i, iterms, keys, why + f": row {i} of D x = A"))
        data = (why, minors, rows)
        self._schur_cache[(P, Q)] = data
        return data

    # ------------------------------------------------------------- R5

    def _build_laplace_expansions(self):
        """All (P, Q, [(coef, (rows, cols)), ...]) row and column Laplace expansions of det B_{P,Q}."""
        B = [[int(v) for v in row] for row in self.B]
        return [(P, Q, [(par * B[i][j], key) for par, i, j, key in terms])
                for P, Q, terms in _laplace_keys(self.n)]

    def rule_laplace(self, chi, log) -> bool:
        """det B_{P,Q} = Σ_j ± b_ij det B_{P∖i, Q∖j} along every row i (and every column)."""
        forced = False
        for terms, _, target_key, why in self.laplace_sums():
            forced |= self.enforce_sum(chi, terms, self.minor(chi, *target_key), target_key, log, why)
        return forced

    def laplace_sums(self):
        """R5 as signed sums (terms, None, (P, Q), why); the target is the sign of det B_PQ."""
        for P, Q, terms in self._laplace_expansions:
            yield terms, None, (P, Q), f"R5 Laplace expansion of det B_{P},{Q}"

    # ------------------------------------------------------------- R7

    def _build_circuits(self):
        """
        For every (n+1)-subset C of the 2n elements of N = [I | Bᵀ]: the list of
        (coef, (Q, P)) with  T(e_k) C(e_k) = coef · sign det B_{Q,P}  for its n+1
        elements, where T = (+ⁿ, sign a) is the feasibility tope.
        """
        T = [1] * self.n + list(self.a)
        circuits = []
        for keyed in _circuit_keys(self.n):
            terms = [(T[e] * coef, key) for coef, e, key in keyed]
            iterms, _ = self.compile_terms(terms)
            circuits.append((terms, iterms))
        return circuits

    def rule_tope(self, chi, log) -> bool:
        """(+ⁿ, sign a) must be orthogonal to every circuit of the oriented matroid of [I | Bᵀ]."""
        forced = False
        for circuit in self._circuits:
            forced |= self.tope_circuit(chi, log, circuit)
        return forced

    def tope_circuit(self, chi, log, circuit) -> bool:
        """R7 on one circuit: (terms, iterms) from _build_circuits, terms = [(coef, (Q, P))]."""
        terms, iterms = circuit
        signs = [coef * c * chi[S] for coef, S, c in iterms]
        n_open = signs.count(UNKNOWN)
        known = {v for v in signs if v != UNKNOWN}
        if len(known) != 1 or n_open > 1:
            return False
        (common,) = known
        used = [k for (_, k), v in zip(terms, signs) if v != UNKNOWN]
        if not n_open:
            raise Conflict(f"R7 tope: circuit {[k for _, k in terms]} is one-sided against (+, sign a)", used)
        coef, key = terms[signs.index(UNKNOWN)]
        return self.set_minor(chi, *key, -common * coef, log, "R7 tope: circuit orthogonality", used)

    # ------------------------------------------------------------- R8

    def rule_real_eigenvalue(self, chi, log) -> bool:
        """
        R8: S B P_σ is a Z-matrix for some σ ≠ id (so s = −1 exactly on the species σ moves,
        +1 on its fixed points), and either s·a > 0, or s_i a_i < 0 for exactly one species i and
        σ has an even number of cycles ⇒ −diag(x*) B has a positive real eigenvalue.

        Proof.  With X = diag(x*), D = S B X P_σ is a Z-matrix and D (P_σ⁻¹ 1) = S a.  If S a > 0,
        D is a nonsingular M-matrix.  Otherwise let G be the species with s_i a_i > 0: the G-block
        of D has positive row sums (the entries outside the block are ≤ 0), so it is a nonsingular
        M-matrix with positive principal minors; and det D = (−1)^k det B det X > 0 for k even
        under R1 (det B > 0).  Ordering G first, every leading principal minor of the Z-matrix D
        is positive, so D is again a nonsingular M-matrix.  For t ≥ 0, D_t = D + t S P_σ adds −t
        at the off-diagonal positions (σ(j), j) of the moved species and +t on the diagonal at
        fixed points, so it stays a Z-matrix; its least real eigenvalue τ(t) is continuous with
        τ(0) > 0.  For y the indicator of the moved species, (D_t y)_i = (D y)_i − t on moved rows
        and ≤ 0 on fixed rows, so D_t y ≤ 0 for large t and τ(t) ≤ 0.  Hence
        D_{t₁} = S (B X + t₁ I) P_σ is singular for some t₁ > 0: −t₁ is an eigenvalue of B X, so
        t₁ > 0 is an eigenvalue of the Jacobian −X B (similar to −B X).
        """
        Q = self.species
        data = self.schur_block(chi, (), Q)
        if data is None:
            return False
        _, D, A = data
        if any(v == UNKNOWN for v in D.values()) or any(v == UNKNOWN for v in A.values()):
            return False
        for s, sigma in self.z_twists(Q, D):
            moved = [i for i in Q if sigma[i] != i]
            if not moved:
                continue
            wrong = [i for i in Q if s[i] * A[i] < 0]
            if len(wrong) == 0:
                why = "S B P_σ is an M-matrix with σ ≠ id"
            elif len(wrong) == 1 and self._cycle_count(sigma, moved) % 2 == 0:
                why = (f"S B P_σ is a Z-matrix whose block off species {wrong[0]} is an M-matrix and "
                       f"det(S B P_σ) = det B > 0 (σ has an even number of cycles)")
            else:
                continue
            raise Conflict(f"R8 real eigenvalue s={tuple(s[i] for i in Q)} σ={tuple(sigma[i] for i in Q)}: "
                           f"{why}, so −diag(x*) B has a positive real eigenvalue",
                           list(self._block_premises))
        return False

    @staticmethod
    def _cycle_count(sigma: Dict[int, int], moved: Sequence[int]) -> int:
        """Number of cycles of the permutation σ on the species it moves."""
        seen, k = set(), 0
        for i in moved:
            if i in seen:
                continue
            k += 1
            j = i
            while j not in seen:
                seen.add(j)
                j = sigma[j]
        return k

    # ---------------------------------------------------------- driver

    RULES = {
        "R1": "rule_stability",
        "R2": "rule_cramer",
        "R2n": "rule_cramer_full",
        "R3": "rule_balanced",
        "R4": "rule_z_block",
        "R5": "rule_laplace",
        "R6": "rule_schur",
        "R7": "rule_tope",
        "R8": "rule_real_eigenvalue",
    }
    ALL_RULES = ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")

    def rule_cramer_full(self, chi, log) -> bool:
        """R2 restricted to J = [n]: the feasibility inequalities of the algorithm of the Annals paper (R1 + R2n)."""
        return self.rule_cramer(chi, log, full_system_only=True)

    def apply_rules(self, chi, log, use: Sequence[str] = ALL_RULES):
        """Run the selected rules on ``chi`` until no sign changes. Raises Conflict."""
        rules = [getattr(self, self.RULES[name]) for name in use]
        changed = True
        while changed:
            changed = False
            for rule in rules:
                changed |= bool(rule(chi, log))
        return chi


def _pm(sign: int) -> str:
    return "> 0" if sign > 0 else "< 0"


# ------------------------------------------------- pattern-independent structure, cached per n

@lru_cache(maxsize=None)
def _laplace_keys(n: int):
    """(P, Q, [(parity, i, j, (rows, cols))]) for every row and column expansion of det B_{P,Q}."""
    species = tuple(range(n))
    expansions = []
    for k in range(2, n + 1):
        for P in combinations(species, k):
            for Q in combinations(species, k):
                for along_rows in (True, False):
                    for pr, r in enumerate(P if along_rows else Q):
                        terms = []
                        for pc, c in enumerate(Q if along_rows else P):
                            i, j = (r, c) if along_rows else (c, r)
                            rows = tuple(x for x in P if x != i)
                            cols = tuple(x for x in Q if x != j)
                            terms.append(((-1 if (pr + pc) % 2 else 1), i, j, (rows, cols)))
                        expansions.append((P, Q, terms))
    return expansions


@lru_cache(maxsize=None)
def _circuit_keys(n: int):
    """For every (n+1)-subset C of the 2n elements: [(coef without the tope factor T(e), e, (Q, P))]."""
    circuits = []
    for C in combinations(range(2 * n), n + 1):
        terms = []
        for k, e in enumerate(C):
            S = tuple(x for x in C if x != e)
            T1 = tuple(x for x in S if x < n)
            Q = tuple(x - n for x in S if x >= n)
            P = tuple(x for x in range(n) if x not in T1)
            eps = -1 if (sum(T1) - len(T1) * (len(T1) - 1) // 2) % 2 else 1
            terms.append(((-1 if k % 2 else 1) * eps, e, (Q, P)))
        circuits.append(terms)
    return circuits
