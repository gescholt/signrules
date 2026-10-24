"""
Which sign rules do the root kills actually use?

For every pattern class of n species (up to relabelling) killed at the root, the certificate
is the backward slice of the root log from the contradiction (``certificate.certificate_steps``).
The script counts, per rule, in how many certificates it occurs and how often it is the rule
that reports the contradiction, and then drops one rule at a time to see how many root kills
(and how many impossibility verdicts, root or search) survive without it.

    python scripts/rule_usage.py 3
    python scripts/rule_usage.py 4                 # 3044 classes, about a minute on 4 cores
    python scripts/rule_usage.py 4 --engine sweep  # the same with the sweep engine's certificates
"""

import sys
import argparse
from collections import Counter
from multiprocessing import Pool
from pathlib import Path


from signrules.patterns import canonical_classes  # noqa: E402
from signrules import Ecology, complete  # noqa: E402
from signrules.certificate import certificate_steps, GP_RULE  # noqa: E402

RULES = Ecology.ALL_RULES


def _rule(why: str) -> str:
    return "GP" if why == GP_RULE or why.startswith("3-term") else why.split()[0]


def _job(args):
    (a, B), engine = args
    eco = Ecology(a, B)
    res = complete(eco, max_completions=1, keep_derivations=False, engine=engine)
    out = dict(root=res.root_conflict is not None, impossible=not res.completions)
    if out["root"]:
        steps, _, _ = certificate_steps(eco, res)
        out["used"] = sorted({_rule(why) for _, _, why, _ in steps} | {_rule(res.root_conflict)})
        out["kill"] = _rule(res.root_conflict)
        out["size"] = len(steps)
    if out["impossible"]:
        out["without"] = {}
        for r in RULES:
            use = tuple(x for x in RULES if x != r)
            root = complete(eco, use=use, max_nodes=0, engine=engine).root_conflict is not None
            full = root or not complete(eco, use=use, max_completions=1, keep_derivations=False,
                                        engine=engine).completions
            out["without"][r] = (root, full)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int, nargs="?", default=3)
    ap.add_argument("--engine", default="trigger", choices=["trigger", "sweep"])
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    classes = canonical_classes(args.n)
    with Pool(args.jobs) as pool:
        rows = pool.map(_job, [(c, args.engine) for c in classes], chunksize=8)
    roots = [r for r in rows if r["root"]]
    imp = [r for r in rows if r["impossible"]]
    print(f"n={args.n}: {len(rows)} classes, {len(imp)} impossible by the rules, {len(roots)} of them at the root"
          f"  (engine: {args.engine})")
    if roots:
        sizes = sorted(r["size"] for r in roots)
        print(f"certificate size (forced signs in the slice): median {sizes[len(sizes) // 2]}, max {sizes[-1]}\n")
    used = Counter(x for r in roots for x in r["used"])
    kill = Counter(r["kill"] for r in roots)
    print(f"{'rule':5s} {'in certificates':>16s} {'reports the kill':>17s} {'root kills lost':>16s} {'verdicts lost':>14s}"
          "   (lost: when this rule alone is dropped)")
    for r in ("GP",) + RULES:
        lost_root = sum(1 for x in roots if not x["without"][r][0]) if r != "GP" else "-"
        lost_full = sum(1 for x in imp if not x["without"][r][1]) if r != "GP" else "-"
        print(f"{r:5s} {used[r]:16d} {kill[r]:17d} {lost_root:>16} {lost_full:>14}")


if __name__ == "__main__":
    main()
