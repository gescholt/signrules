"""
Candidate screen for impossible n-species patterns by one-species extension of the impossible
(n-1)-species classes, and the sign rules on every candidate.

Deletion lemma.  If a non-trivial pattern (a, B) is impossible, then deleting any one species k
leaves an impossible pattern (trivial ones count as impossible): if the deletion were possible,
add k back with x_k = eps, weak effects of k on the others and a strong self-regulation
b_kk = M / eps, M chosen so that (Bx)_k = a_k (for a_k < 0 this uses a benefactor of k, which
exists since the pattern is non-trivial); the Jacobian is block triangular up to O(eps), so the
extension is possible.

So every non-trivial impossible n-species pattern is one of
  (i)  an extension (a_n, row n, column n: 2^(2n-1) choices) of a non-trivial impossible
       (n-1)-species class, or
  (ii) a pattern all of whose deletions are trivial: then every species has a < 0 and exactly
       one benefactor, the benefactor map is a derangement (the obligate family of R8).
``candidates`` generates both, keeps the classes up to relabelling whose n deletions are all
impossible or trivial; ``run`` runs the sign rules (root, then branching) on each.

The (n-1)-species impossible classes (``base``):
  census  every class of canonical_classes(n-1) on which complete() finds no completion; complete
          when n-1 <= 4, since the rules and the tropical witnesses of the repository feasible_ecologies
          decide every pattern there
  known   the 24 minimal impossible five-species classes, the killed_root rows of
          signrules/data/n5/extension_screen.jsonl (n = 6 only); complete, since the n = 5 screen
          decided every candidate
  FILE    a JSON-lines file with fields a, B (e.g. the impossible classes of a previous run)

    signrules screen 5 --workers 12                  # from the n=4 census
    signrules screen 6 --base known --workers 12
    signrules screen 6 --candidates-only             # count candidates only

The output goes to signrules/data/n<n>/extension_screen.jsonl (``--out`` for another file): one record
per candidate with ``code``, ``a``, ``B``, ``status`` (killed_root / killed_search / open / undecided),
``info``, ``seconds``. What the rules leave open is decided by certificates of possibility outside this
package (the tropical witnesses of the repository feasible_ecologies).
"""
import argparse
import json
import time
from multiprocessing import Pool
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .data import DATA_DIR, read_jsonl, root_kills
from .patterns import canonical, canonical_classes, decode, encode, key, perms, trivial

Base = List[Tuple[list, np.ndarray]]


def _census_job(ab):
    from . import Ecology, complete
    a, B = ab
    return not complete(Ecology(a, B), max_completions=1, keep_derivations=False).completions


def base_classes(n: int, base: str, workers: int = 4) -> Base:
    """Non-trivial impossible (n-1)-species classes, as (a, B) pairs."""
    m = n - 1
    if base == "census":
        if m > 4:
            raise SystemExit("--base census needs n - 1 <= 4 (the census is decided there)")
        classes = [(list(a), np.array(B)) for a, B in canonical_classes(m) if not trivial(a, B)]
        with Pool(workers) as pool:
            imp = pool.map(_census_job, classes, chunksize=8)
        return [c for c, i in zip(classes, imp) if i]
    if base == "known":
        if m != 5:
            raise SystemExit("--base known is the five-species list (n = 6)")
        return [(p["a"], np.array(p["B"])) for p in root_kills(5)]
    pats = read_jsonl(base)
    return [(p["a"], np.array(p["B"])) for p in pats if len(p["a"]) == m and not trivial(p["a"], p["B"])]


def candidates(n: int, base: Base):
    """(raw count, classes up to relabelling, candidate codes) for a list of base classes."""
    keys = {key(a, B) for a, B in base}
    m = n - 1
    raw = []
    for a0, B0 in base:                                              # (i) extensions
        for bits in range(1 << (2 * m + 1)):
            a = list(a0) + [1 if bits & 1 else -1]
            B = np.eye(n, dtype=int)
            B[:m, :m] = B0
            for j in range(m):
                B[m, j] = 1 if (bits >> (1 + j)) & 1 else -1
                B[j, m] = 1 if (bits >> (1 + m + j)) & 1 else -1
            raw.append(encode(a, B))
    for p in perms(n):                                               # (ii) derangement family
        if all(p[i] != i for i in range(n)):
            B = np.ones((n, n), dtype=int)
            for i in range(n):
                B[i, p[i]] = -1
            raw.append(encode([-1] * n, B))
    reps = np.unique(canonical(raw, n))
    out = []
    for code in reps:
        a, B = decode(int(code), n)
        if trivial(a, B):
            continue
        ok = True
        for k in range(n):
            keep = [i for i in range(n) if i != k]
            a1, B1 = [a[i] for i in keep], B[np.ix_(keep, keep)]
            if not (trivial(a1, B1) or key(a1, B1) in keys):
                ok = False
                break
        if ok:
            out.append(int(code))
    return len(raw), len(reps), out


def screen_job(args):
    """The sign rules on one candidate: root propagation, then a search of at most max_nodes nodes."""
    code, n, max_nodes = args
    from . import Ecology, complete
    a, B = decode(code, n)
    t0 = time.time()
    res = complete(Ecology(a, B), max_completions=1, max_nodes=0, keep_derivations=False)
    rec = dict(code=code, a=a, B=B.tolist())
    if res.root_conflict is not None:
        rec.update(status="killed_root", info=res.root_conflict)
    else:
        res = complete(Ecology(a, B), max_completions=1, max_nodes=max_nodes, keep_derivations=False)
        if res.completions:
            rec.update(status="open", info="sign-consistent completion found")
        elif res.stopped:
            rec.update(status="undecided", info=f"node limit {max_nodes}")
        else:
            rec.update(status="killed_search", info=f"{res.nodes} nodes")
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def run(n: int, base: Optional[str] = None, workers: int = 4, max_nodes: int = 20000, out=None,
        candidates_only: bool = False, verbose: bool = True) -> List[dict]:
    """The whole screen; returns the records (and writes them to ``out`` unless candidates_only)."""
    t0 = time.time()
    base_list = base_classes(n, base or ("census" if n <= 5 else "known"), workers)
    if verbose:
        print(f"{len(base_list)} non-trivial impossible {n - 1}-species classes [{time.time() - t0:.0f}s]", flush=True)
    n_raw, n_reps, cands = candidates(n, base_list)
    if verbose:
        print(f"{n_raw} raw extensions, {n_reps} up to relabelling, {len(cands)} non-trivial candidates "
              f"with all deletions impossible or trivial [{time.time() - t0:.0f}s]", flush=True)
    if candidates_only:
        return [dict(code=c, n=n) for c in cands]
    out = DATA_DIR / f"n{n}" / "extension_screen.jsonl" if out is None else out
    out.parent.mkdir(parents=True, exist_ok=True)
    counts, recs = {}, []
    with Pool(workers) as pool, open(out, "w") as f:
        for i, rec in enumerate(pool.imap_unordered(screen_job, [(c, n, max_nodes) for c in cands]), 1):
            f.write(json.dumps(rec) + "\n")
            f.flush()
            recs.append(rec)
            counts[rec["status"]] = counts.get(rec["status"], 0) + 1
            if verbose and (i % 50 == 0 or i == len(cands)):
                print(f"{i}/{len(cands)} {counts} [{time.time() - t0:.0f}s]", flush=True)
    return recs


def add_arguments(ap: argparse.ArgumentParser) -> None:
    """The options of the screen (shared by ``signrules screen`` and ``signrules-screen``)."""
    ap.add_argument("n", type=int, help="number of species of the candidates (5 or 6)")
    ap.add_argument("--base", default=None, help="census | known | FILE (default: census for n <= 5, known for n = 6)")
    ap.add_argument("--workers", type=int, default=4, help="parallel processes (default 4)")
    ap.add_argument("--max-nodes", type=int, default=20000, help="node limit of the search on each candidate (default 20000)")
    ap.add_argument("--candidates-only", action="store_true", help="count the candidates, do not run the rules")
    ap.add_argument("--out", default=None, help="output JSONL (default: signrules/data/n<n>/extension_screen.jsonl, "
                                                 "which overwrites the shipped list)")


def run_args(args) -> None:
    """Run the screen for parsed ``args`` (see ``add_arguments``)."""
    from pathlib import Path
    run(args.n, args.base, args.workers, args.max_nodes, Path(args.out) if args.out else None, args.candidates_only)


def main(argv: Optional[Sequence[str]] = None):
    """Command line of the screen (``signrules-screen``, the same as ``signrules screen``)."""
    ap = argparse.ArgumentParser(prog="signrules-screen", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    run_args(ap.parse_args(argv))


if __name__ == "__main__":
    main()
