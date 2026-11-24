"""
Three-term Grassmann–Plücker sign propagation for uniform chirotopes.

For every (n-2)-set C and a<b<c<d outside C the Plücker coordinates satisfy

    p_{abC} p_{cdC} - p_{acC} p_{bdC} + p_{adC} p_{bcC} = 0 .

With all coordinates nonzero, the three terms cannot all have the same sign.
So if two terms are known and agree, the third is forced to the opposite
sign; if its factors are known and disagree with that, the table is
inconsistent.  For uniform chirotopes the 3-term relations imply the full
chirotope axioms (Björner et al., Thm 3.6.2), so this is the complete
combinatorial content of "the signs come from a point of Gr(n,2n)" that
sign propagation can extract.

``build_relations(n)`` lists the relations with bases (for display and for code that
keeps its own tables); ``build_relation_ids(n)`` lists them with basis ids, which is what
``propagate`` and ``violates`` take together with a flat table (see plucker_index.py).
"""

from functools import lru_cache
from itertools import combinations
from typing import List, Optional, Tuple

from .plucker_index import perm_sign, Basis, index

# A term is (basis1, basis2, coefficient) meaning coefficient * p_basis1 * p_basis2
Term = Tuple[Basis, Basis, int]
Relation = Tuple[Term, Term, Term]
IdTerm = Tuple[int, int, int]
IdRelation = Tuple[IdTerm, IdTerm, IdTerm]


@lru_cache(maxsize=None)
def build_relations(n: int) -> List[Relation]:
    """
    All 3-term relations of Gr(n, 2n), as triples of (basis, basis, ±1) terms with sorted bases.
    Cached per n: the list is shared, do not modify it.
    """
    N = 2 * n
    rels: List[Relation] = []

    def key(seq):
        return tuple(sorted(seq)), perm_sign(seq)

    for C in combinations(range(N), n - 2):
        rest = [e for e in range(N) if e not in C]
        Cl = list(C)
        for a, b, c, d in combinations(rest, 4):
            t = []
            for (u, v), (w, z), coef in (((a, b), (c, d), 1), ((a, c), (b, d), -1), ((a, d), (b, c), 1)):
                k1, s1 = key([u, v] + Cl)
                k2, s2 = key([w, z] + Cl)
                t.append((k1, k2, coef * s1 * s2))
            rels.append((t[0], t[1], t[2]))
    return rels


@lru_cache(maxsize=None)
def build_relation_ids(n: int) -> List[IdRelation]:
    """The relations of ``build_relations(n)``, in the same order, with basis ids."""
    bid = index(n).id
    return [tuple((bid[k1], bid[k2], c) for k1, k2, c in rel) for rel in build_relations(n)]


def propagate(chi: List[int], relations: List[IdRelation], log: Optional[list] = None,
              bases: Optional[List[Basis]] = None):
    """
    Unit propagation over the 3-term relations, in place, to a fixed point.
    Returns (consistent, n_inferred); consistent=False means all three terms
    of some relation are known and have the same sign — the offending relation
    (with ids) is then available as ``propagate.last_conflict``.  If ``log`` is given,
    every forced sign is appended as (basis, sign, premises) with premises the
    four other bases of the relation, as bases when ``bases`` is given and as ids otherwise.
    """
    inferred = 0
    changed = True
    propagate.last_conflict = None
    name = bases.__getitem__ if bases is not None else (lambda i: i)
    while changed:
        changed = False
        for rel in relations:
            (a1, a2, ca), (b1, b2, cb), (c1, c2, cc) = rel
            ta = ca * chi[a1] * chi[a2]
            tb = cb * chi[b1] * chi[b2]
            tc = cc * chi[c1] * chi[c2]
            if ta and tb and tc:
                if ta == tb == tc:
                    propagate.last_conflict = rel
                    return False, inferred
                continue
            if ta == 0:
                if tb == 0 or tc == 0 or tb != tc:
                    continue
                k1, k2, coef, forced = a1, a2, ca, -tb
            elif tb == 0:
                if tc == 0 or ta != tc:
                    continue
                k1, k2, coef, forced = b1, b2, cb, -ta
            else:
                if ta != tb:
                    continue
                k1, k2, coef, forced = c1, c2, cc, -ta
            v1, v2 = chi[k1], chi[k2]
            if v1 == 0 and v2 == 0:
                continue
            if v1 == 0:
                chi[k1] = forced * coef * v2
                new = k1
            else:
                chi[k2] = forced * coef * v1
                new = k2
            if log is not None:
                others = [name(b) for t in rel for b in t[:2] if b != new]
                log.append((name(new), chi[new], others))
            inferred += 1
            changed = True
    return True, inferred


def violates(chi: List[int], relations: List[IdRelation]) -> bool:
    """True if a complete table violates some 3-term relation."""
    for (a1, a2, ca), (b1, b2, cb), (c1, c2, cc) in relations:
        ta = ca * chi[a1] * chi[a2]
        if ta == cb * chi[b1] * chi[b2] == cc * chi[c1] * chi[c2]:
            return True
    return False
