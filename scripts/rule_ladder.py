"""
What each rule buys: for a ladder of rule sets, how many of the minimal impossible patterns of the
extension screens are killed at the root, killed after branching, or left open.

For n = 4 the eleven patterns of scripts/impossible_eleven.py, for n = 5 and 6 the patterns with
status killed_root in signrules/data/n<n>/extension_screen.jsonl (24 and 41).  For each rule set:

    root       complete(max_nodes=0) ends in a conflict at the root;
    search     otherwise the depth-first search with --max-nodes nodes finds no completion and
               terminates: killed, but only by branching;
    open       a sign completion consistent with the rules exists (so the rule set can never kill it);
    stopped    the node budget ran out before a verdict.

    python scripts/rule_ladder.py 4 5 --workers 8                  # JSONL + the LaTeX table
    python scripts/rule_ladder.py 4 5 6 --workers 8 --max-nodes 2000

Writes signrules/data/rule_ladder.jsonl (one record per (n, rule set, pattern)) and
docs/rule_ladder_table.tex (rows = rule sets, columns = root / search / open per n).
"""

import sys
import json
import time
import argparse
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
from signrules.data import DATA_DIR

# R4 and R6 are subsumed by the other rules (docs/SIGNRULES.md "R4 and R6 are subsumed", bead Python-7d3)
# and no longer in the ladder; the earlier rows (R1,R4), (R1,R6) are in git history.
LADDER = [("R1", "R2n"), ("R1", "R2"), ("R1", "R7"),
          ("R1", "R2", "R3", "R5", "R7"), ("R1", "R2", "R3", "R5", "R7", "R8")]


def patterns(n):
    if n == 4:
        from signrules.patterns import ELEVEN
        return [(name, a, np.array(B)) for name, (a, B) in ELEVEN.items()]
    recs = [json.loads(l) for l in open(DATA_DIR / f"n{n}" / "extension_screen.jsonl") if l.strip()]
    return [(str(r["code"]), r["a"], np.array(r["B"])) for r in recs if r["status"] == "killed_root"]


def _job(args):
    n, use, name, a, B, max_nodes = args
    from signrules import Ecology, complete
    eco = Ecology(a, B, name)
    t0 = time.time()
    res = complete(eco, use=use, max_nodes=0, max_completions=1, keep_derivations=False)
    rec = dict(n=n, rules=",".join(use), name=name)
    if res.root_conflict is not None:
        rec.update(verdict="root", info=res.root_conflict.split(":")[0], open_at_root=0)
    else:
        rec["open_at_root"] = res.root_unknowns
        res = complete(eco, use=use, max_nodes=max_nodes, max_completions=1, keep_derivations=False)
        if res.completions:
            rec.update(verdict="open", info=f"completion after {res.nodes} nodes")
        elif res.stopped:
            rec.update(verdict="stopped", info=f"{max_nodes} nodes, {res.conflicts} pruned")
        else:
            rec.update(verdict="search", info=f"{res.nodes} nodes, {res.conflicts} pruned")
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def rules_tex(use):
    if len(use) >= 5:
        return "all but R8" if "R8" not in use else "all rules"
    return ", ".join(use)


def table(recs, ns, max_nodes):
    cols = "@{\\hspace{1.5em}}".join("rrr" for _ in ns)
    head1 = " & ".join(f"\\multicolumn{{3}}{{c}}{{$n={n}$ ($N={sum(1 for r in recs if r['n'] == n and r['rules'] == recs[0]['rules'])}$)}}" for n in ns)
    head2 = " & ".join("root & search & open" for _ in ns)
    rules = "".join(f"\\cmidrule(lr){{{3 * k + 2}-{3 * k + 4}}}" for k in range(len(ns)))
    lines = ["\\begin{tabular}{@{}l" + cols + "@{}}", "\\toprule", "rules & " + head1 + "\\\\", rules, " & " + head2 + "\\\\", "\\midrule"]
    for use in LADDER:
        key = ",".join(use)
        cells = []
        for n in ns:
            sub = [r for r in recs if r["n"] == n and r["rules"] == key]
            root = sum(r["verdict"] == "root" for r in sub)
            search = sum(r["verdict"] == "search" for r in sub)
            open_ = sum(r["verdict"] == "open" for r in sub)
            stopped = sum(r["verdict"] == "stopped" for r in sub)
            cells += [str(root), str(search), str(open_) + (f"$+{stopped}^*$" if stopped else "")]
        lines.append(rules_tex(use) + " & " + " & ".join(cells) + "\\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", type=int, nargs="+")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-nodes", type=int, default=2000)
    ap.add_argument("--out", default=str(DATA_DIR / "rule_ladder.jsonl"))
    ap.add_argument("--tex", default=str(ROOT / "docs" / "rule_ladder_table.tex"))
    ap.add_argument("--tex-only", action="store_true", help="rebuild the table from --out without rerunning the search")
    args = ap.parse_args()
    if args.tex_only:
        recs = [json.loads(l) for l in open(args.out) if l.strip()]
        recs = [r for r in recs if r["n"] in args.ns]
    else:
        jobs = [(n, use, name, a, B, args.max_nodes) for n in args.ns for use in LADDER for name, a, B in patterns(n)]
        t0 = time.time()
        recs = []
        with Pool(args.workers) as pool:
            for k, rec in enumerate(pool.imap_unordered(_job, jobs), 1):
                recs.append(rec)
                if k % 25 == 0 or k == len(jobs):
                    print(f"{k}/{len(jobs)} [{time.time() - t0:.0f}s]", flush=True)
        recs.sort(key=lambda r: (r["n"], LADDER.index(tuple(r["rules"].split(","))), r["name"]))
        with open(args.out, "w") as f:
            for r in recs:
                f.write(json.dumps(r) + "\n")
    for n in args.ns:
        print(f"n={n}")
        for use in LADDER:
            sub = [r for r in recs if r["n"] == n and r["rules"] == ",".join(use)]
            counts = {v: sum(r["verdict"] == v for r in sub) for v in ("root", "search", "open", "stopped")}
            print(f"  {','.join(use):40s} {counts}  max {max(r['seconds'] for r in sub):.1f}s")
    tex = ("% Generated by scripts/rule_ladder.py -- do not edit by hand.\n"
           f"% Minimal impossible patterns of the extension screens (n=4: the eleven; n=5, 6: the root kills of\n"
           f"% extension_screen.jsonl) under each rule set: killed at the root, killed by a search of at most\n"
           f"% {args.max_nodes} nodes, or open (a sign completion exists; * = node budget exhausted).\n")
    with open(args.tex, "w") as f:
        f.write(tex + table(recs, args.ns, args.max_nodes))
    print("wrote", args.tex if args.tex_only else f"{args.out} and {args.tex}")


if __name__ == "__main__":
    main()
