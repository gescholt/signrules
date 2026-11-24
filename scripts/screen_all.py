"""
Screen every non-trivial five-species ecology (up to relabeling) with the sign rules.

Stage 1: canonical representatives of the 2^20 off-diagonal sign patterns of B under S_5.
Stage 2: for each, the 32 growth-rate sign vectors up to the stabilizer of B; drop trivial
         ecologies (a species with a_i < 0 whose interactions are all deleterious).
Stage 3: run chirotope propagation + the rules (``--rules``, default all of R1-R8) at the root only (no branching).

Every pattern reported as killed is a proved impossible ecology (the rules are sound);
"open" only means root propagation did not close it.

Output: one JSON line per ecology with a, B, code (20-bit encoding of the off-diagonal
signs), status ("killed" / "open") and info (the kill reason, or the number of signs still
unknown after root propagation).  The default ``--out`` is signrules/data/n5/screen_results.jsonl (not
shipped: 248,436 lines).  The run of 2026-09 with ``--rules R1,R2,R3,R4,R5,R6`` took about 100 CPU-minutes
and killed 21 patterns (data/n5/killed_patterns.json); with R8 three more are killed.  Resumable.

    python scripts/screen_all.py --workers 2 --out results.jsonl --seconds 3000
    python scripts/screen_all.py --workers 2 --out results.jsonl --resume      # continue
    python scripts/screen_all.py --rules R1,R2,R3,R4,R5,R6 --workers 8          # the historical run
"""

import json
import time
import argparse
from itertools import permutations, product
from pathlib import Path
from multiprocessing import Pool


import numpy as np

RECORDS = Path(__file__).resolve().parents[1] / "data"          # the records of the repository (not in the package)

N = 5
RULES = None                                                         # set from --rules; None = Ecology.ALL_RULES
OFF = [(i, j) for i in range(N) for j in range(N) if i != j]          # 20 off-diagonal slots
PERMS = list(permutations(range(N)))                                 # 120


def canonical_B_classes():
    """Representatives (as 20-bit ints, bit=1 means +) of B patterns modulo relabeling."""
    m = len(OFF)
    codes = np.arange(1 << m, dtype=np.int64)
    bits = ((codes[:, None] >> np.arange(m)) & 1).astype(np.int8)     # (2^20, 20)
    slot = {ij: k for k, ij in enumerate(OFF)}
    best = codes.copy()
    for p in PERMS:
        # relabel species: entry (i,j) moves to (p[i], p[j])
        idx = [slot[(p[i], p[j])] for (i, j) in OFF]
        permuted = bits[:, idx]
        code = (permuted.astype(np.int64) * (1 << np.arange(m))).sum(axis=1)
        np.minimum(best, code, out=best)
    reps = np.unique(best)
    return reps, bits


def decode_B(code):
    B = np.eye(N, dtype=int)
    for k, (i, j) in enumerate(OFF):
        B[i, j] = 1 if (code >> k) & 1 else -1
    return B


def stabilizer(B):
    return [p for p in PERMS if all(B[p[i], p[j]] == B[i, j] for i, j in OFF)]


def trivial(a, B):
    return any(a[i] < 0 and all(B[i, j] > 0 for j in range(N) if j != i) for i in range(N))


def a_classes(B):
    stab = stabilizer(B)
    seen, out = set(), []
    for a in product([1, -1], repeat=N):
        if a in seen:
            continue
        orbit = {tuple(a[p.index(i)] for i in range(N)) for p in stab}   # a_{p(i)} = a_i
        seen |= orbit
        if not trivial(a, B):
            out.append(min(orbit))
    return out


_ctx = {}


def _init():
    from signrules import build_relation_ids
    _ctx["rels"] = build_relation_ids(N)


def _root(job):
    from signrules import Ecology, Conflict, propagate
    from signrules.plucker_index import UNKNOWN
    a, code = job
    eco = Ecology(list(a), decode_B(code))            # the index of Gr(N, 2N) is shared per N
    chi = eco.initial_table()
    try:
        while True:
            ok, _ = propagate(chi, _ctx["rels"])
            if not ok:
                return a, int(code), "killed", "3-term Grassmann-Plücker"
            before = chi.count(UNKNOWN)
            eco.apply_rules(chi, [], RULES or eco.ALL_RULES)
            after = chi.count(UNKNOWN)
            if after == before:
                return a, int(code), "open", after
    except Conflict as c:
        return a, int(code), "killed", str(c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default=str(RECORDS / "n5" / "screen_results.jsonl"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--resume", action="store_true", help="skip ecologies already present in --out")
    ap.add_argument("--seconds", type=float, default=None, help="stop cleanly after this many seconds")
    ap.add_argument("--rules", default=None, help="comma-separated rule names (default: all of R1-R8)")
    args = ap.parse_args()
    global RULES
    if args.rules:
        RULES = tuple(r for r in args.rules.split(",") if r)

    t0 = time.time()
    reps, _ = canonical_B_classes()
    print(f"{len(reps)} B-patterns up to relabeling [{time.time()-t0:.0f}s]", flush=True)
    jobs = []
    for code in reps:
        B = decode_B(int(code))
        for a in a_classes(B):
            jobs.append((a, int(code)))
    print(f"{len(jobs)} non-trivial ecologies up to relabeling [{time.time()-t0:.0f}s]", flush=True)
    if args.limit:
        jobs = jobs[: args.limit]

    done = set()
    if args.resume and Path(args.out).exists():
        with open(args.out) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                done.add((tuple(r["a"]), r["code"]))
        jobs = [j for j in jobs if j not in done]
        print(f"resuming: {len(done)} done, {len(jobs)} remaining", flush=True)

    killed = 0
    mode = "a" if args.resume else "w"
    with open(args.out, mode) as f, Pool(args.workers, initializer=_init) as pool:
        it = pool.imap_unordered(_root, jobs, chunksize=32)
        for k, (a, code, status, info) in enumerate(it, 1):
            f.write(json.dumps({"a": list(a), "code": code, "B": decode_B(code).tolist(),
                                "status": status, "info": info}) + "\n")
            killed += status == "killed"
            if k % 5000 == 0:
                print(f"  {k}/{len(jobs)} done, {killed} killed [{time.time()-t0:.0f}s]", flush=True)
            if args.seconds and time.time() - t0 > args.seconds:
                f.flush()
                print(f"time budget reached after {k} ecologies ({killed} killed); rerun with --resume", flush=True)
                pool.terminate()
                return
    print(f"done: {len(jobs)} ecologies, {killed} proved impossible at the root, "
          f"{len(jobs)-killed} open [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
