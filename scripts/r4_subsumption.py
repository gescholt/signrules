"""
R4 is subsumed, exhaustively per n.

The trigger engine's propagators are monotone: for tables T ⊆ T' (T' consistent), the closure of T is
contained in the closure of T'.  So to show that on every table closed under the other rules R4 forces
nothing, it suffices to start from the *minimal* table on which R4 applies — exactly its premises — and
check that the other rules' closure already contains R4's conclusions.  The premises are pattern-level
data: a block (P, Q), the sign of det B_PP, the entries of the Schur block D = (B/B_PP)_QQ along a
Z-twist (τ, π), the terms of the reduced right-hand side A_Q (which fix the minors det B_{P+i-r, P} and,
when R = [n] ∖ (P ∪ Q) is nonempty, the entries b_{rt} = -a_r for t ∈ R), and the growth signs a.  The
3-term relations and R7 read nothing else of the pattern, so a configuration

    (n, P, Q, π, case ∈ {M, notM}, sgn det B_PP, a ∈ {±1}^n)

determines the minimal table, and the check is finite: up to relabelling P = {0..p-1}, Q = {p..p+q-1}.

    M case:    the closure under GP + R7 contains every principal and almost-principal minor of Z;
    notM case: for every ordering of Q and every position j0, the closure of (premises + the other
               leading minors positive) contains "leading minor j0 negative" or is inconsistent, and
               the premises + all leading minors positive are inconsistent.

R4 is *not* subsumed at the fixpoint: in the notM case the last unknown leading minor may be a small
one that the other rules leave open (at n = 5 already, P = ∅, |Q| = 4).  What holds instead is
subsumption one branching step deeper: a conclusion of R4 missing from the closure is "ok by probe"
when its opposite sign conflicts at the fixpoint of the other rules (for the notM case this is one
table per ordering, all leading minors positive, which R4 itself rejects).  The search never needs R4
then: branching on that minor refutes the wrong side immediately, and the census of whole trees
(``scripts/rule_ablation.py --census``, ``tests/test_subsumption.py``) shows that the DFS, which decides small minors first, never even reaches such a table.

A configuration that fails is printed; with ``--others`` the failing ones are retried under the full rule
set on explicit patterns (R2, R3, R5 and R8 read the pattern).  Writes signrules/data/r4_subsumption/<n>_*.json.

    python scripts/r4_subsumption.py 3 4 --workers 8 --wlog --rules R1,R7 --others
    python scripts/r4_subsumption.py 5 --workers 16 --wlog --rules R1,R7 --others --patterns 256
    python scripts/r4_subsumption.py 6 --workers 16 --symmetry --rules R1,R7 --others --patterns 64
"""
import json
import time
import random
import argparse
from itertools import combinations, permutations, product
from multiprocessing import Pool

import numpy as np
from signrules.data import DATA_DIR

from signrules import Ecology, Conflict
from signrules.engine import Engine
from signrules.plucker_index import UNKNOWN

OUT = DATA_DIR / "r4_subsumption"


def premises(eco, P, Q, sigma, case, dPP):
    """
    The minimal table of R4 on the block (P, Q) with twist column map sigma and the given case:
    a flat table with only det B_PP, the entries of D, the minors of the reduced right-hand side and
    the entries b_rt = -a_r (t in R) set.  Returns (chi, s, A) or None if the premises clash
    (a right-hand side term that cannot be read).
    """
    n = eco.n
    chi = [UNKNOWN] * eco.idx.m
    S0, c0 = eco._minor_key[((), ())]                 # the empty minor is 1 in every table
    chi[S0] = c0
    for i in range(n):                                # b_ii > 0 in every table
        S, c = eco._minor_key[((i,), (i,))]
        chi[S] = c
    dPP_key, entries, rhs = eco._block_data(P, Q)

    def put(S, v):
        if chi[S] == UNKNOWN:
            chi[S] = v
        return chi[S] == v

    if dPP_key is not None:
        S, c = dPP_key
        if not put(S, c * dPP):
            return None
    else:
        dPP = 1
    A = {}
    for i, (mi, terms) in zip(Q, rhs):
        if terms is None:                               # some w_r unreadable: the pattern must have b_rt = -a_r
            return None
        t_i = next((par if key is None else par * c * chi[S]
                    for par, S, c, key in terms if key is None or key == (P, P)), None)
        if t_i is None or not all(S is None or put(S, par * c * t_i) for par, S, c, key in terms):
            return None
        A[i] = mi * t_i * dPP
    s = {i: (A[i] if case == "M" else -A[i]) for i in Q}
    for i, j, S, c, key, bit in entries:
        if not put(S, c * s[i] * (1 if sigma[i] == j else -1) * dPP):
            return None
    R = [t for t in range(n) if t not in P and t not in Q]
    for r in list(P) + list(Q):
        for t in R:
            S, c = eco._minor_key[((r,), (t,))]
            if not put(S, c * (-eco.a[r])):
                return None
    return chi, s, A


def z_minor_keys(eco, P, Q, s, sigma):
    """(key, sign) of every principal and almost-principal minor of Z = S D P_sigma, in minors of B."""
    out = []
    dPP = 1
    if P:
        S, c = eco._minor_key[(P, P)]
        dPP = None                                      # filled by the caller from the table
    for k in range(1, len(Q) + 1):
        for I in combinations(Q, k):
            coef, key = eco.block_minor(None, P, s, sigma, I, I)
            out.append((key, coef, +1))
    for k in range(len(Q) - 1):
        for I in combinations(Q, k):
            others = [q for q in Q if q not in I]
            for i, j in permutations(others, 2):
                rows, cols = tuple(sorted(I + (i,))), tuple(sorted(I + (j,)))
                coef, key = eco.block_minor(None, P, s, sigma, rows, cols)
                coef *= Ecology_merge(I, i) * Ecology_merge(I, j)
                out.append((key, coef, -1))
    return out


def Ecology_merge(I, i):
    from signrules.rules import merge_sign
    return merge_sign(I, (i,))


def closure(eco, use, chi, engine=None):
    """Close ``chi`` in place under GP and ``use``; None if inconsistent.  ``engine`` is reused when given."""
    try:
        (engine or Engine(eco, use=use)).fixpoint(chi, [])
    except Conflict:
        return None
    return chi


def check_config(args):
    n, P, Q, perm, case, dPP, a, use = args
    sigma = dict(zip(Q, perm))
    B = np.ones((n, n), dtype=int)                      # entries other than the premises play no role,
    for r in P + Q:                                     # but w_r must be readable in the pattern: b_rt = -a_r
        for t in range(n):
            if t not in P and t not in Q:
                B[r, t] = -a[r]
    eco = Ecology(list(a), B)
    pre = premises(eco, P, Q, sigma, case, dPP)
    if pre is None:
        return dict(status="no premises")
    chi0, s, A = pre
    gp = closure(eco, (), list(chi0))
    if gp is None:
        return dict(status="premises inconsistent")
    log = []
    try:
        fired = eco.z_block(list(gp), log, P, Q)       # R4's conclusions on the minimal table
    except Conflict:
        return dict(status="R4 conflicts on the premises")
    r = verdict(eco, chi0, P, Q, sigma, case, dPP, s, use)
    if r == "inconsistent":
        return dict(status="inconsistent under the others")
    if r in ("ok", "ok by probe"):
        return dict(status=r)
    return dict(status="R4 forces more", missing=r)


def check_patterns(args):
    """
    A failing configuration retried on explicit patterns consistent with its premises (the free entries
    of B run over ±1: all of them when there are at most ``max_patterns``, else that many at random),
    from the pattern's initial table plus the premises, under the given rules (all but R4: R2, R3, R5
    and R8 read the pattern).  Returns counts and the bad patterns.
    """
    n, f, use, max_patterns = args
    P, Q, perm, case, dPP, a = (tuple(f["P"]), tuple(f["Q"]), tuple(f["sigma"]), f["case"], f["dPP"], tuple(f["a"]))
    sigma = dict(zip(Q, perm))
    B0 = np.ones((n, n), dtype=int)
    for r in P + Q:
        for t in range(n):
            if t not in P and t not in Q:
                B0[r, t] = -a[r]
    eco0 = Ecology(list(a), B0)
    pre = premises(eco0, P, Q, sigma, case, dPP)
    if pre is None:
        return dict(patterns=0, bad=[])
    chi0, s, A = pre
    fixed = {}
    for i in range(n):
        for j in range(n):
            S, c = eco0._minor_key[((i,), (j,))]
            if chi0[S] != UNKNOWN:
                fixed[(i, j)] = c * chi0[S]
    free = [(i, j) for i in range(n) for j in range(n) if i != j and (i, j) not in fixed
            and not (j not in P and j not in Q and i in P + Q)]
    counts = dict(patterns=0, ok=0, probe=0, inconsistent=0, bad=[], free=len(free), sampled=False)
    if 2 ** len(free) <= max_patterns:
        pattern_signs = product((1, -1), repeat=len(free))
    else:                                               # a random sample of the explicit patterns
        rng = random.Random(repr((n, P, Q, perm, case, dPP, a)))
        pattern_signs = (tuple(rng.choice((1, -1)) for _ in free) for _ in range(max_patterns))
        counts["sampled"] = True
    for signs in pattern_signs:
        B = B0.copy()
        for (i, j), v in fixed.items():
            B[i, j] = v
        for (i, j), v in zip(free, signs):
            B[i, j] = v
        eco = Ecology(list(a), B)
        chi = eco.initial_table()
        clash = False
        for S, v in enumerate(chi0):
            if v == UNKNOWN:
                continue
            if chi[S] == UNKNOWN:
                chi[S] = v
            elif chi[S] != v:
                clash = True
        if clash:
            continue
        counts["patterns"] += 1
        r = verdict(eco, chi, P, Q, sigma, case, dPP, s, use)
        if r == "ok":
            counts["ok"] += 1
        elif r == "ok by probe":
            counts["probe"] += 1
        elif r == "inconsistent":
            counts["inconsistent"] += 1
        else:
            counts["bad"].append((B.tolist(), r))
    counts["bad"] = counts["bad"][:3]
    return counts


def refuted(eco, use, engine, table, S, v):
    """Does setting chi[S] = v on the closed ``table`` conflict at the fixpoint of ``use``?"""
    chi = list(table)
    chi[S] = v
    return closure(eco, use, chi, engine) is None


def verdict(eco, chi0, P, Q, sigma, case, dPP, s, use):
    """
    The M or notM check on a given starting table.  Returns 'ok' when R4's conclusions are in the
    closure under ``use``, 'ok by probe' when the missing ones are refuted one step deeper (their
    opposite sign conflicts at the fixpoint), 'inconsistent' when the premises already conflict, and
    otherwise the reason (a conclusion the other rules neither derive nor refute by a probe).
    """
    engine = Engine(eco, use=use)
    probed = False
    if case == "M":
        table = closure(eco, use, list(chi0), engine)
        if table is None:
            return "inconsistent"
        dPP_sign = 1
        if P:
            S, c = eco._minor_key[(P, P)]
            dPP_sign = c * table[S]
        for key, coef, zsign in z_minor_keys(eco, P, Q, s, sigma):
            S, c = eco._minor_key[key]
            want = c * coef * dPP_sign * zsign          # chi value R4 forces
            if table[S] == want:
                continue
            if table[S] != UNKNOWN or not refuted(eco, use, engine, table, S, -want):
                return f"missing {key}"
            probed = True
        return "ok by probe" if probed else "ok"
    dPP_sign = dPP if P else 1
    for order in permutations(Q):
        leading = []
        for m in range(1, len(Q) + 1):
            I = tuple(sorted(order[:m]))
            leading.append(eco.block_minor(None, P, s, sigma, I, I))
        tables = []
        for j0 in range(len(leading) + 1):              # j0 == len(leading): all leading minors positive
            chi = list(chi0)
            ok = True
            for j, (coef, key) in enumerate(leading):
                if j == j0:
                    continue
                S, c = eco._minor_key[key]
                v = c * coef * dPP_sign
                if chi[S] == UNKNOWN:
                    chi[S] = v
                elif chi[S] != v:
                    ok = False
            if not ok:
                tables.append(None)
                continue
            tables.append(closure(eco, use, chi, engine))
        if tables[-1] is not None:
            return f"all leading minors positive is consistent for {order}"
        for j0, table in enumerate(tables[:-1]):        # all positive conflicts: every j0 is at least probed
            if table is None:
                continue
            coef, key = leading[j0]
            S, c = eco._minor_key[key]
            if c * table[S] != -coef * dPP_sign:
                probed = True
    return "ok by probe" if probed else "ok"


def canonical_pairs(n, P, Q):
    """
    Representatives of the pairs (sigma, a), sigma a permutation of Q, a in {±1}^n, under the
    relabellings that fix P, Q and R = [n] - P - Q setwise: sigma up to conjugation (one permutation per
    cycle type), then a up to the centralizer of sigma in Sym(Q) times Sym(P) times Sym(R).
    """
    R = tuple(i for i in range(n) if i not in P and i not in Q)
    reps = []
    seen_types = set()
    for perm in permutations(Q):
        sigma = dict(zip(Q, perm))
        # cycle type
        rest, cycles = set(Q), []
        while rest:
            i = min(rest); c = 0
            while i in rest:
                rest.remove(i); i = sigma[i]; c += 1
            cycles.append(c)
        ctype = tuple(sorted(cycles))
        if ctype in seen_types:
            continue
        seen_types.add(ctype)
        centralizer = [g for g in permutations(Q) if all(dict(zip(Q, g))[sigma[i]] == sigma[dict(zip(Q, g))[i]] for i in Q)]
        group = [(dict(zip(P, gp)), dict(zip(Q, gq)), dict(zip(R, gr)))
                 for gp in permutations(P) for gq in centralizer for gr in permutations(R)]
        seen_a = set()
        for a in product((1, -1), repeat=n):
            key = min(tuple(a[{**gp, **gq, **gr}[i]] for i in range(n)) for gp, gq, gr in group)
            if key in seen_a:
                continue
            seen_a.add(key)
            reps.append((perm, a))
    return reps


def configs(n, wlog, symmetry=False):
    species = tuple(range(n))
    for p in range(0, n - 1):
        for q in range(2, n - p + 1):
            if wlog or symmetry:
                blocks = [(tuple(range(p)), tuple(range(p, p + q)))]
            else:
                blocks = [(P, Q) for P in combinations(species, p)
                          for Q in combinations([i for i in species if i not in P], q)]
            for P, Q in blocks:
                pairs = canonical_pairs(n, P, Q) if symmetry else [(perm, a) for perm in permutations(Q)
                                                                   for a in product((1, -1), repeat=n)]
                for perm, a in pairs:
                    for case in ("M", "notM"):
                        for dPP in ((1, -1) if P else (1,)):
                            yield n, P, Q, perm, case, dPP, a


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ns", nargs="+", type=int)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--wlog", action="store_true", help="P, Q as initial segments (relabelling)")
    ap.add_argument("--symmetry", action="store_true",
                    help="also sigma up to conjugation and a up to its centralizer (relabelling within P, Q, R)")
    ap.add_argument("--rules", default="R7", help="rules besides the 3-term relations, e.g. R7 or R2,R7")
    ap.add_argument("--others", action="store_true",
                    help="retry failing configurations on explicit patterns under every rule but R4")
    ap.add_argument("--patterns", type=int, default=1024,
                    help="with --others: at most this many explicit patterns per configuration (random beyond)")
    ap.add_argument("--sample", type=int, default=0,
                    help="with --others: retry only this many failing configurations, chosen at random (seed 0)")
    args = ap.parse_args()
    use = tuple(r for r in args.rules.split(",") if r)
    OUT.mkdir(exist_ok=True)
    for n in args.ns:
        t0 = time.time()
        jobs = [c + (use,) for c in configs(n, args.wlog, args.symmetry)]
        with Pool(args.workers) as pool:
            res = pool.map(check_config, jobs, chunksize=64)
        counts = {}
        failures = []
        for job, r in zip(jobs, res):
            counts[r["status"]] = counts.get(r["status"], 0) + 1
            if r["status"] == "R4 forces more":
                n_, P, Q, perm, case, dPP, a, _ = job
                failures.append(dict(P=P, Q=Q, sigma=perm, case=case, dPP=dPP, a=a, missing=str(r["missing"])))
        print(f"n={n}: {len(jobs)} configurations in {time.time() - t0:.0f} s: {counts}")
        for f in failures[:10]:
            print("   ", f)
        retried = []
        if args.others and failures:
            t1 = time.time()
            chosen = failures
            if args.sample and args.sample < len(failures):
                chosen = sorted(random.Random(0).sample(failures, args.sample), key=failures.index)
            jobs2 = [(n, f, tuple(r for r in Ecology.ALL_RULES if r != "R4"), args.patterns) for f in chosen]
            with Pool(args.workers) as pool:
                res2 = pool.map(check_patterns, jobs2, chunksize=1)
            for f, r in zip(chosen, res2):
                f["patterns"] = r
                if r["bad"]:
                    retried.append(f)
            totals = {k: sum(f["patterns"][k] for f in chosen) for k in ("patterns", "ok", "probe", "inconsistent")}
            print(f"   retried on explicit patterns under all other rules ({time.time() - t1:.0f} s): "
                  f"{len(chosen) - len(retried)} of {len(chosen)} configurations ok on every pattern; "
                  f"patterns {totals}")
            for f in retried[:10]:
                print("    still failing:", f)
        stem = f"{n}_{'_'.join(('GP',) + use)}" + ("_sym" if args.symmetry else "_wlog" if args.wlog else "")
        with open(OUT / f"{stem}.json", "w") as f:
            json.dump(dict(n=n, rules=["GP"] + list(use), wlog=args.wlog, symmetry=args.symmetry, counts=counts,
                           failures=failures, retried=len(chosen) if args.others and failures else 0,
                           still_failing=len(retried) if args.others else None), f)


if __name__ == "__main__":
    main()
