"""
Draw the binary search trees of the sign-completion search.

`complete(eco, record_tree=True)` records every node: the basis branched on, the sign chosen,
how many further signs the propagation forced, how many stay open, and whether the node is a
completion or a conflict (with the rule).  This script plots the tree with matplotlib (PNG) and,
if plotly is installed, as an interactive HTML page with hover details, without Graphviz.

    python scripts/search_trees_n4.py                       # the built-in cases into signrules/data/n4/trees/
    python scripts/search_trees_n4.py --only S-d --rules R1,R2n --out out/
    python scripts/search_trees_n4.py --n3                  # an open three-species pattern
    python scripts/search_trees_n4.py --only A-c --rules R1,R2 --no-labels   # 647 nodes

Built-in cases: the eleven candidates of Meng et al. under a chosen rule set (all rules kill them at
the root, so the interesting trees come from the weaker sets R1,R2n and R1,R2), and an open n = 3 pattern.
"""

import sys
import argparse
from pathlib import Path


from signrules import Ecology, complete, search_tree_text, save_search_tree_png
from signrules.patterns import ELEVEN

try:
    from signrules import save_search_tree_html
    import plotly  # noqa: F401
    HAVE_PLOTLY = True
except ImportError:
    HAVE_PLOTLY = False

N3 = {"n3 open (+,-,+)": ([1, -1, 1], [[1, -1, 1], [1, 1, -1], [-1, 1, 1]])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", default="R1,R2n", help="comma-separated rule names (default R1,R2n)")
    ap.add_argument("--only", default=None, help="run only the patterns whose name contains this string")
    ap.add_argument("--n3", action="store_true", help="the open n=3 pattern instead of the eleven")
    ap.add_argument("--max", type=int, default=None, help="stop after this many completions")
    ap.add_argument("--out", default="out/n4_trees")
    ap.add_argument("--no-labels", action="store_true")
    ap.add_argument("--minors", action="store_true", help="label nodes by minors of B instead of Plücker coordinates")
    ap.add_argument("--text-depth", type=int, default=3, help="depth of the text tree printed (0 = none)")
    args = ap.parse_args()
    use = tuple(r for r in args.rules.split(",") if r)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cases = N3 if args.n3 else ELEVEN
    for name, (a, B) in cases.items():
        if args.only and args.only not in name:
            continue
        eco = Ecology(a, B, name)
        res = complete(eco, use=use, max_completions=args.max, record_tree=True)
        print(f"{name:45s} {res.summary()}  tree: {len(res.tree)} nodes")
        if res.root_conflict:
            continue
        if args.text_depth:
            print(search_tree_text(eco, res, max_depth=args.text_depth, plucker=not args.minors))
        stem = out / (name.split()[0] + "_" + "-".join(use))
        kw = dict(labels=not args.no_labels if args.no_labels else "auto", plucker=not args.minors, use=use)
        save_search_tree_png(eco, res, f"{stem}.png", **kw)
        print(f"      wrote {stem}.png")
        if HAVE_PLOTLY:
            save_search_tree_html(eco, res, f"{stem}.html", **kw)
            print(f"      wrote {stem}.html")


if __name__ == "__main__":
    main()
