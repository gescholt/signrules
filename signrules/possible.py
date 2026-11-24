"""
Certificates of possibility, checked in exact rational arithmetic.

A sign pattern ``(a, B)`` is *possible* when some matrix with the signs of ``B`` and some growth rates with
the signs of ``a`` have a positive, locally asymptotically stable equilibrium.  A certificate is a rational
matrix ``C`` with the signs of ``B`` whose row sums have the signs of ``a``: then ``a = C·1`` makes
``x* = (1, …, 1)`` an equilibrium of ``ẋ = diag(x)(a − Cx)``, its Jacobian is ``−C``, and it is stable when
every eigenvalue of ``C`` has positive real part.  ``check_matrix`` decides all of this exactly, with the
Routh–Hurwitz criterion on ``det(λI + C)``; nothing is trusted to floating point.

The shipped certificates (``signrules.data.possible(n)``, n = 5, 6) are tropical witnesses: integer exponents
``U`` and an integer ``T``, standing for ``c_ij = b_ij · T^(−u_ij)``.  They certify the extension candidates
that the rules leave open, 59 at n = 5 and 102 at n = 6, so that with the 24 and 41 root kills every candidate
is decided.  How such witnesses are found is the subject of the second paper; checking one needs only this
module.

    from signrules.data import possible
    from signrules.possible import check_witness
    assert all(check_witness(r) for r in possible(5))
"""
from fractions import Fraction
from typing import List, Sequence

Matrix = List[List[Fraction]]


def _sign(x) -> int:
    return (x > 0) - (x < 0)


def _det(M: Matrix) -> Fraction:
    """Determinant by fraction-exact Gaussian elimination."""
    A = [row[:] for row in M]
    n, d = len(A), Fraction(1)
    for k in range(n):
        p = next((i for i in range(k, n) if A[i][k] != 0), None)
        if p is None:
            return Fraction(0)
        if p != k:
            A[k], A[p] = A[p], A[k]
            d = -d
        d *= A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            if f:
                A[i] = [x - f * y for x, y in zip(A[i], A[k])]
    return d


def charpoly(C: Matrix) -> List[Fraction]:
    """Coefficients ``[1, c_1, …, c_n]`` of ``det(λI + C)`` (Faddeev–LeVerrier, exact)."""
    n = len(C)
    coeffs = [Fraction(1)]
    M = [[Fraction(0)] * n for _ in range(n)]                       # M_0 = 0
    for k in range(1, n + 1):
        # det(λI − A) with A = −C:  M_k = A M_{k-1} + c_{k-1} I,  c_k = −tr(A M_k) / k
        M = [[-sum(C[i][l] * M[l][j] for l in range(n)) + (coeffs[-1] if i == j else 0) for j in range(n)]
             for i in range(n)]
        coeffs.append(sum(C[i][l] * M[l][i] for i in range(n) for l in range(n)) / k)
    return coeffs


def hurwitz_stable(C: Matrix) -> bool:
    """Every eigenvalue of ``C`` has positive real part: the Hurwitz determinants of ``det(λI + C)`` are > 0."""
    c = charpoly(C)
    n = len(C)
    if any(x <= 0 for x in c):
        return False
    H = [[c[2 * j - i + 1] if 0 <= 2 * j - i + 1 <= n else Fraction(0) for j in range(n)] for i in range(n)]
    return all(_det([row[:k] for row in H[:k]]) > 0 for k in range(1, n + 1))


def check_matrix(a: Sequence[int], B, C: Matrix) -> bool:
    """``C`` certifies that the pattern ``(a, B)`` is possible (see the module docstring)."""
    n = len(a)
    C = [[Fraction(x) for x in row] for row in C]
    if any(_sign(C[i][j]) != B[i][j] for i in range(n) for j in range(n)):
        return False
    if any(_sign(sum(C[i])) != a[i] for i in range(n)):
        return False
    return hurwitz_stable(C)


def witness_matrix(B, U, T: int) -> Matrix:
    """The matrix ``c_ij = b_ij · T^(−u_ij)`` of a tropical witness."""
    T = Fraction(T)
    return [[s * T ** (-u) for s, u in zip(srow, urow)] for srow, urow in zip(B, U)]


def check_witness(rec: dict) -> bool:
    """A record with ``a``, ``B``, ``U``, ``T`` (as in ``signrules.data.possible``) certifies its pattern."""
    return check_matrix(rec["a"], rec["B"], witness_matrix(rec["B"], rec["U"], rec["T"]))


__all__ = ["check_matrix", "check_witness", "witness_matrix", "hurwitz_stable", "charpoly"]
