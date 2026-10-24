"""
Search trees and certificates for the five-species ecologies proved impossible.

The 21 kills of the n = 5 screen (`killed_patterns.json`) all die at the root under the full rule
set, so their search tree is a single node and the proof is the certificate DAG.  Under the classical
rules R1+R2n (feasibility + weak stability, the algorithm of the Grassmannian paper) they survive
root propagation with 150–175 open signs and the branching search finds completions without end.
This script writes, for each selected pattern,

  <id>_certificate.png/.pdf/.dot   the root kill under all rules (derivation DAG)
  <id>_R1-R2n.png/.html            the classical search tree, stopped after --max completions
                                   (matplotlib PNG, plotly HTML)
  <id>_<rules>.png/.html      optionally the tree under another rule set (--rules)

    python scripts/search_trees_n5.py                         # the ten kills not covered by the papers' theorems
    python scripts/search_trees_n5.py --ids 21 7 --max 16     # chosen patterns
    python scripts/search_trees_n5.py --rules R1,R2,R5 --max-nodes 300
"""

import sys
import json
import argparse
from pathlib import Path


from signrules import (Ecology, complete, write_certificate, search_tree_text,
from signrules.data import DATA_DIR
                           save_search_tree_png, save_search_tree_html)

NEW_TEN = [2, 3, 4, 5, 6, 15, 16, 17, 18, 19]     # not covered by Theorems 3.1–3.6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", type=int, nargs="*", default=NEW_TEN, help="pattern ids from killed_patterns.json")
    ap.add_argument("--rules", default="R1,R2n", help="rule set for the branching tree (default: the classical R1,R2n)")
    ap.add_argument("--max", type=int, default=32, help="stop the tree after this many completions")
    ap.add_argument("--max-nodes", type=int, default=None, help="or after this many nodes")
    ap.add_argument("--out", default="out/n5_trees")
    ap.add_argument("--no-certificates", action="store_true")
    args = ap.parse_args()
    use = tuple(r for r in args.rules.split(",") if r)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    patterns = {k["id"]: k for k in json.load(open(DATA_DIR / "n5" / "killed_patterns.json"))}

    for pid in args.ids:
        k = patterns[pid]
        label = f"#{pid} {k['structure']}"
        eco = Ecology(k["a"], k["B"], label)
        if not args.no_certificates:
            full = complete(eco)
            assert full.root_conflict, f"#{pid} not killed at the root"
            write_certificate(eco, full, str(out / f"{pid:02d}_certificate"))
        res = complete(eco, use=use, record_tree=True, max_completions=args.max, max_nodes=args.max_nodes)
        tag = "-".join(use)
        print(f"{label:70s} {tag}: {res.summary()}; {res.root_unknowns} open at the root")
        if res.root_conflict:
            continue
        print(search_tree_text(eco, res, max_depth=2))
        stem = out / f"{pid:02d}_{tag}"
        save_search_tree_png(eco, res, f"{stem}.png", use=use)
        try:
            save_search_tree_html(eco, res, f"{stem}.html", use=use)
        except ImportError:
            pass


if __name__ == "__main__":
    main()
