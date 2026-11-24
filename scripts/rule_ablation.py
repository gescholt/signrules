"""
Ablation of the Schur-block rules R4 and R6: rerun the searches of the
extension screens under variants of the rule set and compare the trees.

Patterns: the eleven of scripts/impossible_eleven.py and every row of
signrules/data/n<n>/extension_screen.jsonl for n = 5, 6 (83 + 143).  Each search is what the screen
runs (``signrules.screen.screen_job``): complete(max_nodes=0), then, if open at the root,
complete(max_completions=1, max_nodes=20000).  Variants:

    all       every rule (the baseline)
    root46    R4 and R6 at the root only (complete(root_only=("R4", "R6"))): the propagators stay
              registered, so the branching order is unchanged and only the pruning can differ
    root4, root6
    drop46    R4 and R6 left out of ``use``: also changes the branching order (speed ceiling)
    drop4, drop6

    python scripts/rule_ablation.py 4 5 6 --variant all root46 drop46 --workers 8
    python scripts/rule_ablation.py --compare data/rule_ablation/all.jsonl data/rule_ablation/root46.jsonl
    python scripts/rule_ablation.py 3 4 --census --max-nodes 3000 --workers 20
    python scripts/rule_ablation.py 5 --census --random 20000 --max-nodes 2000 --workers 20

Writes data/rule_ablation/<variant>.jsonl, one record per pattern: status, nodes, conflicts,
the hash of the first completion, seconds and the per-rule counts (runs, signs set, conflicts
raised; at the root and below).  ``--compare`` joins two files on (n, name), lists every pattern
whose verdict or tree differs, sums the seconds and the R4/R6 counts, and exits 1 on a changed
verdict.

``--census`` runs every pattern — all 512 at n = 3, the 3044 relabelling classes at
n = 4, or ``--random COUNT`` patterns — under the full rule set with R4 and R6 really running
(``idle=()``), whole trees up to ``--max-nodes`` and ``--max-completions``, and counts what R4 and R6
contribute: they sit at the highest priority levels, so each of their runs sees a table closed under
GP, R1, R2, R3, R5, R7, R8, and a nonzero count is a pattern where the other rules do not subsume
them (docs/SIGNRULES.md, "R4 and R6 are subsumed"; none found).  Records census<n>[_random].jsonl,
with the R4/R6 counts only; exit code 1 if any pattern fired.  Shipped records: all.jsonl, root46.jsonl,
drop46.jsonl, the single-rule variants root4.jsonl, root6.jsonl, drop4.jsonl, drop6.jsonl, census3.jsonl,
census4.jsonl; the random census files for n = 5, 6 are not committed (10 MB),
their summary is in docs/SIGNRULES.md §3 ("R4 and R6 are subsumed").
"""

import sys
import json
import time
import hashlib
import argparse
from multiprocessing import Pool

import numpy as np
from pathlib import Path

from signrules.data import DATA_DIR

RECORDS = Path(__file__).resolve().parents[1] / "data"          # the records of the repository (not in the package)

OUT = RECORDS / "rule_ablation"
ALL = ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
VARIANTS = {                                 # idle=(): every variant really runs R4 and R6 where it names them
    "all": dict(use=ALL, root_only=(), idle=()),
    "root46": dict(use=ALL, root_only=("R4", "R6"), idle=()),
    "root4": dict(use=ALL, root_only=("R4",), idle=()),
    "root6": dict(use=ALL, root_only=("R6",), idle=()),
    "drop46": dict(use=tuple(r for r in ALL if r not in ("R4", "R6")), root_only=(), idle=()),
    "drop4": dict(use=tuple(r for r in ALL if r != "R4"), root_only=(), idle=()),
    "drop6": dict(use=tuple(r for r in ALL if r != "R6"), root_only=(), idle=()),
}
TREE_KEYS = ("status", "root_unknowns", "nodes", "conflicts", "completion")


def patterns(n):
    if n == 4:
        from signrules.patterns import ELEVEN
        return [(name, a, np.array(B)) for name, (a, B) in ELEVEN.items()]
    recs = [json.loads(l) for l in open(DATA_DIR / f"n{n}" / "extension_screen.jsonl") if l.strip()]
    return [(str(r["code"]), r["a"], np.array(r["B"])) for r in recs]


def merge_stats(*stats):
    out = {}
    for st in stats:
        for rule, counts in (st or {}).items():
            tot = out.setdefault(rule, {w: dict(runs=0, signs=0, conflicts=0) for w in ("root", "below")})
            for where, c in counts.items():
                for k, v in c.items():
                    tot[where][k] += v
    return out


def _job(args):
    n, variant, name, a, B, max_nodes = args
    from signrules import Ecology, complete
    opts = dict(VARIANTS[variant], stats=True, keep_derivations=False, max_completions=1)
    eco = Ecology(a, B, name)
    t0 = time.time()
    res = complete(eco, max_nodes=0, **opts)
    rec = dict(n=n, name=name, variant=variant, root_unknowns=res.root_unknowns,
               nodes=0, conflicts=0, stopped=False, completion=None)
    stats = [res.rule_stats]
    if res.root_conflict is not None:
        rec.update(status="killed_root", info=res.root_conflict.split(":")[0])
    else:
        res = complete(eco, max_nodes=max_nodes, **opts)
        stats.append(res.rule_stats)
        rec.update(nodes=res.nodes, conflicts=res.conflicts, stopped=res.stopped)
        if res.completions:
            rec.update(status="open", info=f"completion after {res.nodes} nodes",
                       completion=hashlib.sha1(bytes(x % 256 for x in res.completions[0])).hexdigest()[:12])
        elif res.stopped:
            rec.update(status="undecided", info=f"node limit {max_nodes}, {res.conflicts} pruned")
        else:
            rec.update(status="killed_search", info=f"{res.nodes} nodes, {res.conflicts} pruned")
    rec["seconds"] = round(time.time() - t0, 3)
    rec["rule_stats"] = merge_stats(*stats)
    return rec


def run(ns, variants, workers, max_nodes):
    OUT.mkdir(exist_ok=True)
    for variant in variants:
        jobs = [(n, variant, name, a, B, max_nodes) for n in ns for name, a, B in patterns(n)]
        t0 = time.time()
        with Pool(workers) as pool:
            recs = pool.map(_job, jobs, chunksize=1)
        path = OUT / f"{variant}.jsonl"
        with open(path, "w") as f:
            for r in recs:
                f.write(json.dumps(r) + "\n")
        print(f"{variant:8s} {len(recs)} patterns, {sum(r['seconds'] for r in recs):.1f} s of search "
              f"({time.time() - t0:.1f} s wall), -> {path}")
        summary(recs)


def summary(recs):
    tot = merge_stats(*(r["rule_stats"] for r in recs))
    for rule in ("R4", "R6"):
        if rule in tot:
            c = tot[rule]
            print(f"    {rule}: root {c['root']['runs']} runs / {c['root']['signs']} signs / {c['root']['conflicts']} conflicts;"
                  f" below {c['below']['runs']} runs / {c['below']['signs']} signs / {c['below']['conflicts']} conflicts")
    by = {}
    for r in recs:
        by[r["status"]] = by.get(r["status"], 0) + 1
    print("    " + ", ".join(f"{k} {v}" for k, v in sorted(by.items())))


def load(path):
    return {(r["n"], r["name"]): r for r in (json.loads(l) for l in open(path) if l.strip())}


def compare(base_path, other_path):
    base, other = load(base_path), load(other_path)
    keys = sorted(set(base) & set(other))
    missing = set(base) ^ set(other)
    if missing:
        print(f"only in one file: {sorted(missing)}")
    verdicts, trees = [], []
    for k in keys:
        b, o = base[k], other[k]
        if b["status"] != o["status"]:
            verdicts.append(k)
        if any(b[f] != o[f] for f in TREE_KEYS):
            trees.append(k)
    bs, os_ = sum(base[k]["seconds"] for k in keys), sum(other[k]["seconds"] for k in keys)
    bv, ov = next(iter(base.values()))["variant"], next(iter(other.values()))["variant"]
    print(f"{len(keys)} patterns: {bv} {bs:.1f} s, {ov} {os_:.1f} s ({os_ / bs:.2f}x)")
    print(f"verdict changes: {len(verdicts)}; tree changes (status, root_unknowns, nodes, conflicts, first completion): {len(trees)}")
    for k in trees:
        b, o = base[k], other[k]
        print(f"  n={k[0]} {k[1]}: " + "; ".join(f"{f} {b[f]} -> {o[f]}" for f in TREE_KEYS if b[f] != o[f]))
    for name, recs in ((bv, base), (ov, other)):
        print(f"  {name}:")
        summary([recs[k] for k in keys])
    return 1 if verdicts else 0


# ------------------------------------------------------------------ census

def census_patterns(n, random_count, seed):
    """All 512 patterns at n = 3, the 3044 relabelling classes at n = 4, or a random sample."""
    if random_count:
        from signrules.patterns import random_pattern
        seen, out = set(), []
        for i in range(random_count):
            a, B = random_pattern(n, seed + i)
            key = (tuple(a), tuple(map(tuple, B)))
            if key not in seen:
                seen.add(key)
                out.append((f"seed{seed + i}", a, np.array(B)))
        return out
    if n == 3:
        import itertools
        out = []
        for a in itertools.product((1, -1), repeat=3):
            for off in itertools.product((1, -1), repeat=6):
                B = np.eye(3, dtype=int)
                B[~np.eye(3, dtype=bool)] = off
                out.append(("".join("+-"[x < 0] for x in a) + "/" + "".join("+-"[x < 0] for x in off), list(a), B))
        return out
    if n == 4:
        from signrules.patterns import canonical_classes
        return [(f"class{i}", list(a), np.asarray(B)) for i, (a, B) in enumerate(canonical_classes(4))]
    raise SystemExit(f"no exhaustive census at n = {n}: use --random COUNT")


def _census_job(args):
    n, name, a, B, max_nodes, max_completions = args
    from signrules import Ecology, complete
    t0 = time.time()
    res = complete(Ecology(a, B, name), use=ALL, idle=(), stats=True, keep_derivations=False,
                   max_nodes=max_nodes, max_completions=max_completions)
    st = res.rule_stats
    fired = {r: st[r]["root"]["signs"] + st[r]["below"]["signs"] + st[r]["root"]["conflicts"] + st[r]["below"]["conflicts"]
             for r in ("R4", "R6") if r in st}
    status = ("killed_root" if res.root_conflict is not None else "open" if res.completions
              else "undecided" if res.stopped else "killed_search")
    return dict(n=n, name=name, a=list(a), B=np.asarray(B).tolist(), status=status, nodes=res.nodes,
                conflicts=res.conflicts, completions=len(res.completions), stopped=res.stopped,
                seconds=round(time.time() - t0, 3), fired=fired,
                rule_stats={r: st[r] for r in ("R4", "R6") if r in st})     # the other rules' counts are not kept


def census(n, workers, max_nodes, max_completions, random_count, seed):
    """
    Run every pattern under the full rule set and count the signs and conflicts that R4 and R6
    contribute: they sit at the highest priority levels, so each of their runs sees a table closed
    under GP, R1, R2, R5, R7, R8, and a nonzero count is a pattern where the other rules do not
    subsume them.  Writes data/rule_ablation/census<n>[_random].jsonl and lists the offenders.
    """
    OUT.mkdir(exist_ok=True)
    pats = census_patterns(n, random_count, seed)
    jobs = [(n, name, a, B, max_nodes, max_completions) for name, a, B in pats]
    t0 = time.time()
    with Pool(workers) as pool:
        recs = pool.map(_census_job, jobs, chunksize=4)
    path = OUT / (f"census{n}_random.jsonl" if random_count else f"census{n}.jsonl")
    with open(path, "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    print(f"n={n}: {len(recs)} patterns, {sum(r['nodes'] for r in recs)} nodes, "
          f"{sum(r['completions'] for r in recs)} completions, {sum(r['stopped'] for r in recs)} hit the node limit, "
          f"{sum(r['seconds'] for r in recs):.1f} s of search ({time.time() - t0:.1f} s wall) -> {path}")
    summary(recs)
    bad = [r for r in recs if any(r["fired"].values())]
    for r in bad:
        print(f"  R4/R6 fired: {r['name']} a={r['a']} B={r['B']} {r['fired']} status {r['status']}")
    print(f"    patterns where R4 or R6 set a sign or raised a conflict: {len(bad)}")
    return len(bad)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ns", nargs="*", type=int, help="4, 5 and/or 6 (3, 4 or any n with --random for --census)")
    ap.add_argument("--variant", nargs="+", default=["all", "root46"], choices=sorted(VARIANTS))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-nodes", type=int, default=20000)
    ap.add_argument("--compare", nargs=2, metavar=("BASE", "OTHER"))
    ap.add_argument("--census", action="store_true",
                    help="exhaustive census of R4/R6 activity under the full rule set")
    ap.add_argument("--max-completions", type=int, default=10 ** 9, help="census: completions per pattern")
    ap.add_argument("--random", type=int, default=0, metavar="COUNT", help="census: random patterns instead")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if args.compare:
        sys.exit(compare(*args.compare))
    if not args.ns:
        ap.error("give n (4 5 6) or --compare")
    if args.census:
        bad = sum(census(n, args.workers, args.max_nodes, args.max_completions, args.random, args.seed)
                  for n in args.ns)
        sys.exit(1 if bad else 0)
    run(args.ns, args.variant, args.workers, args.max_nodes)


if __name__ == "__main__":
    main()
