"""Sign-level disconnection certificates, exact for any number of completions.

Two distinct minors f, g of B that satisfy sgn f = eps * sgn g in EVERY completion of a pattern (neither
sign constant) vanish nowhere on the feasible-stable set S_sigma, which is then disconnected as soon as
both sides are realized (the repository feasible_ecologies, docs/STABILITY_SETS_DIRECTIONS.md §7, whose
``examples/sign_correlations.py`` reads the relation off an enumerated list of completions, exact only
when the list is complete; at n >= 5 it never is).  This script proves the relation without enumerating:

1. screen: a sample of completions (``--sample``, the first ones the search finds), made to vary in every
   free coordinate: a coordinate constant on the sample is either forced or unseen, and one search for a
   completion with the opposite sign (``assume``) decides which, adding the completion it finds to the
   sample.  The free coordinates are then grouped by their sign vector up to a global sign; a group with
   two or more members is a candidate class.  A pair that differs in one sampled completion is not
   correlated, so with every free coordinate varying the screen misses no class;
2. coarse cells: for the class representative S and each member T, ``complete(decide=[S, T])`` branches
   on those two coordinates only and lists the sign combinations propagation allows (at most four
   leaves); if the two forbidden combinations (chi[T] = -eps chi[S]) are absent the relation holds in
   every completion, by soundness of the rules;
3. search: a forbidden combination that survives propagation is written into the initial table
   (``complete(assume={S: s, T: t}, max_completions=1)``): a completion found refutes the pair (and
   says the sample was incomplete), an exhausted search proves it, a node limit leaves it undecided.

Every basis of the Plücker index is a minor of B up to the fixed signs of the growth rates
(``eco.idx.decomps[S] = (T_a, P, Q, eps)``), so correlation classes of coordinates are classes of minors.

    python scripts/sign_correlations_exact.py 3                       # the 512 labelled patterns
    python scripts/sign_correlations_exact.py 4 --workers 20          # the 3044 classes up to relabelling
    python scripts/sign_correlations_exact.py 4 --file patterns.jsonl    # rows with a, B (e.g. the sampling files of feasible_ecologies)
Output: one JSON line per pattern that is open under the rules, in --out (default
signrules/data/n<n>/sign_correlations_exact.jsonl), and a summary on stdout.
"""
import argparse
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from signrules.data import DATA_DIR

from signrules import Ecology, complete

FORBIDDEN = {1: ((1, -1), (-1, 1)), -1: ((1, 1), (-1, -1))}     # eps -> the (chi[S], chi[T]) with chi[T] != eps chi[S]


def prove_pair(eco, S, T, eps, max_nodes):
    """Is chi[T] = eps * chi[S] in every completion?  Returns (verdict, how, nodes, witness)."""
    r = complete(eco, decide=[S, T], keep_derivations=False)
    cells = {(t[S], t[T]) for t in r.completions}
    present = [c for c in FORBIDDEN[eps] if c in cells]
    nodes = r.nodes
    if not present:
        return "proved", "propagation", nodes, None
    for s, t in present:
        r = complete(eco, assume={S: s, T: t}, max_completions=1, max_nodes=max_nodes, keep_derivations=False)
        nodes += r.nodes
        if r.completions:
            return "refuted", "search", nodes, (s, t)
        if r.stopped:
            return "undecided", "node limit", nodes, (s, t)
    return "proved", "search", nodes, None


def minor_of(eco, S):
    _, P, Q, _ = eco.idx.decomps[S]
    return [list(P), list(Q)]


def analyse(a, B, sample, max_nodes):
    eco = Ecology(a, B)
    t0 = time.time()
    res = complete(eco, max_completions=sample, keep_derivations=False)
    rec = dict(a=list(a), B=np.asarray(B).tolist())
    if res.root_conflict:
        rec.update(status="killed", seconds=time.time() - t0)
        return rec
    comps = list(res.completions)
    n_sampled = len(comps)
    if not comps:
        rec.update(status="impossible", seconds=time.time() - t0)
        return rec
    # freeness, exactly: a coordinate constant on the sample is forced or merely unseen; one search with the
    # opposite sign decides, and a completion found joins the sample (so the pair screen below sees every
    # free coordinate vary)
    forced, undecided_free, searches = 0, [], 0
    for S in range(eco.idx.m):
        v = comps[0][S]
        if any(t[S] != v for t in comps):
            continue
        r = complete(eco, assume={S: -v}, max_completions=1, max_nodes=max_nodes, keep_derivations=False)
        searches += 1
        if r.completions:
            comps.append(r.completions[0])
        elif r.stopped:
            undecided_free.append(S)
        else:
            forced += 1
    cols = list(zip(*comps))
    groups = {}
    for S in range(eco.idx.m):
        v = cols[S]
        if 1 in v and -1 in v:
            key, sign = (v, 1) if v[0] == 1 else (tuple(-x for x in v), -1)
            groups.setdefault(key, []).append((S, sign))
    classes = []
    for key, members in groups.items():
        if len(members) < 2:
            continue
        rep, s_rep = members[0]
        pairs = []
        for T, s_T in members[1:]:
            eps = s_rep * s_T
            verdict, how, nodes, witness = prove_pair(eco, rep, T, eps, max_nodes)
            pairs.append(dict(ids=[rep, T], minors=[minor_of(eco, rep), minor_of(eco, T)], eps=eps,
                              verdict=verdict, how=how, nodes=nodes, witness=witness))
        verdicts = {p["verdict"] for p in pairs}
        status = "proved" if verdicts == {"proved"} else ("refuted" if verdicts == {"refuted"} else "mixed")
        classes.append(dict(size=len(members), plus=sum(x == 1 for x in key), minus=sum(x == -1 for x in key),
                            status=status, pairs=pairs))
    proved = [c for c in classes if any(p["verdict"] == "proved" for p in c["pairs"])]
    rec.update(status="open", sampled=n_sampled, exhaustive=not res.stopped,
               free=sum(len(g) for g in groups.values()), forced=forced, undecided_free=undecided_free,
               freeness_searches=searches, candidate_classes=len(classes), proved_classes=len(proved),
               refuted_pairs=sum(p["verdict"] == "refuted" for c in classes for p in c["pairs"]),
               undecided_pairs=sum(p["verdict"] == "undecided" for c in classes for p in c["pairs"]),
               disconnected=bool(proved), seconds=time.time() - t0,
               classes=[c for c in classes if c["status"] != "refuted"])        # refuted classes: counts only
    return rec


def _job(args):
    a, B, sample, max_nodes = args
    return analyse(a, B, sample, max_nodes)


def patterns(n, file, splits):
    if file:
        rows = [json.loads(line) for line in open(file)]
        if splits:
            rows = [r for r in rows if r.get("samples", 1) >= 50 and r.get("components", 1) > 1]
        return [(r["a"], r["B"]) for r in rows]
    if n == 3:
        out = []
        for a in itertools.product((1, -1), repeat=n):
            for off in itertools.product((1, -1), repeat=n * (n - 1)):
                B = np.eye(n, dtype=int)
                B[~np.eye(n, dtype=bool)] = off
                out.append((list(a), B))
        return out
    from signrules.patterns import canonical_classes
    return [(list(a), np.asarray(B)) for a, B in canonical_classes(n)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", type=int)
    ap.add_argument("--file", help="JSONL rows with a, B")
    ap.add_argument("--splits", action="store_true", help="with --file: only the rows whose sample graph split")
    ap.add_argument("--sample", type=int, default=2000, help="completions in the screening sample (0 = all)")
    ap.add_argument("--max-nodes", type=int, default=200000, help="node limit of each deciding search")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out")
    args = ap.parse_args()
    pats = patterns(args.n, args.file, args.splits)
    out = Path(args.out) if args.out else DATA_DIR / f"n{args.n}/sign_correlations_exact.jsonl"
    sample = args.sample or None
    jobs = [(a, B, sample, args.max_nodes) for a, B in pats]
    summary = dict(patterns=len(pats), killed=0, impossible=0, open=0, exhaustive=0, candidate=0, disconnected=0,
                   refuted_pairs=0, undecided_pairs=0)
    t0 = time.time()
    with open(out, "w") as fh, Pool(args.workers) as pool:
        for k, rec in enumerate(pool.imap_unordered(_job, jobs, chunksize=4), 1):
            summary[rec["status"]] += 1
            if rec["status"] == "open":
                summary["exhaustive"] += rec["exhaustive"]
                summary["candidate"] += rec["candidate_classes"] > 0
                summary["disconnected"] += rec["disconnected"]
                summary["refuted_pairs"] += rec["refuted_pairs"]
                summary["undecided_pairs"] += rec["undecided_pairs"]
                fh.write(json.dumps(rec) + "\n")
                if rec["classes"]:
                    print(f"a={rec['a']} B={rec['B']}: {rec['candidate_classes']} candidate class(es), "
                          f"{rec['proved_classes']} proved, {rec['refuted_pairs']} refuted pair(s), "
                          f"{rec['undecided_pairs']} undecided; sample {rec['sampled']}"
                          f"{' (all)' if rec['exhaustive'] else ''}", flush=True)
            if k % 200 == 0:
                print(f"... {k}/{len(pats)} in {time.time() - t0:.0f} s", file=sys.stderr, flush=True)
    summary["seconds"] = round(time.time() - t0, 1)
    print("SUMMARY", json.dumps(summary))


if __name__ == "__main__":
    main()
