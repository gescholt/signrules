"""
Infinite families of impossible ecologies, checked at the root for several n.

Benefactor-permutation family (not in the papers).  Let σ be a permutation
of [n] with fixed set F and k cycles of length ≥ 2 on N = [n] ∖ F.  Species
i ∈ F benefits from every other species and has a_i > 0 ("universal
beneficiary"); species i ∈ N is obligate (a_i < 0), benefits only from σ(i)
and is harmed by everyone else.  With s = +1 on F and −1 on N, the matrix
D = S B P_σ is a Z-matrix with positive diagonal and D (P_σ⁻¹ x*) = S a > 0,
so D is a nonsingular M-matrix: det D > 0 and D⁻¹ ≥ 0.

(a) No fixed points (every species obligate, σ a derangement, any k):
    B⁻¹ = P_σ D⁻¹ S gives (B⁻¹)_jj = −(D⁻¹)_{σ⁻¹(j), j} ≤ 0 for every j, so every
    principal (n−1)-minor det B_{[n]∖j} = (B⁻¹)_jj det B has the sign opposite
    to det B, and the coefficient c₁ = Σ_j x*_{[n]∖j} det B_{[n]∖j} of the
    characteristic polynomial has the sign of −det B.  Stability needs
    det B > 0 and c₁ > 0.  Impossible.   (k = 1 is Theorem 3.4 of
    Çelik–Isoldi–Portakal–Zucal.)

(b) With fixed points and k odd:  det D = (Π s) sgn(σ) det B = (−1)^k det B > 0
    contradicts det B > 0.  Impossible.   (k = 1, |F| = 1 is Theorem 3.6.)

(c) In general (any k, any F, N ≠ ∅): D_t = S(B X + t I) P_σ is a Z-matrix for t ≥ 0 that is an
    M-matrix at t = 0 and has D_t 1_N ≤ 0 for large t, so it is singular at some t₁ > 0, and
    t₁ > 0 is an eigenvalue of the Jacobian −X B.  Impossible.  This covers k even with fixed
    points, where (a) and (b) say nothing and R1–R7 leave the patterns open; rule R8 is (c).

(d) One wrong growth sign and k even: D_0 = S B X P_σ is still an M-matrix (its block off the
    wrong species has positive row sums, and det D_0 = det B det X > 0), so (c) goes through.
    Rule R8 includes (d); it decides the last open patterns at n = 5 and n = 6.

    (The pattern constructor ``parity_pattern`` and the shapes ``cycle_types`` live in ``signrules.patterns``.)
"""
import argparse
import time

from signrules import Ecology, complete
from signrules.patterns import cycle_types, parity_pattern


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=6)
    args = ap.parse_args()
    for n in range(4, args.nmax + 1):
        print(f"\nn = {n}")
        for cycles, fixed in cycle_types(n):
            a, B = parity_pattern(n, cycles, fixed)
            t = time.time()
            res = complete(Ecology(a, B), use=("R1", "R4"), max_completions=1)
            r8 = complete(Ecology(a, B), use=("R8",), max_nodes=0).root_conflict is not None
            k = len(cycles)
            verdict = "R1+R4 kill at the root" if res.root_conflict else f"R1+R4 open ({res.root_unknowns} unknown)"
            verdict += ", R8 kills" if r8 else ", R8 silent"
            if not fixed:
                expected = "(a) all obligate: impossible"
            elif k % 2:
                expected = "(b) odd k: impossible"
            else:
                expected = "(c) even k with fixed points: impossible"
            print(f"  cycles {[len(c) for c in cycles]!s:14} fixed {len(fixed)}  k={k}  ->  {verdict:45s} [{expected}] {time.time()-t:.1f}s")


if __name__ == "__main__":
    main()
