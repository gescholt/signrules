"""
Sign completions of the preset patterns (``signrules.patterns.PATTERNS``, n = 5 and 6), with their derivation graphs.

Root propagation under R1–R6 leaves 248,415 of the 248,436 five-species ecologies open (docs/SIGNRULES.md §6;
R8 closes three more).  This script takes the presets (a certificate is drawn for those killed at the root), runs the branching search until a couple of completions are
found, and writes for each completion its derivation graph: the root propagation, then one band
per branch choice with every Plücker sign it forces and the rule that forced it.  A completion is a
sign-consistent table, not a proof of possibility.

    python scripts/completions.py                          # built-in patterns, 2 completions each
    python scripts/completions.py --max 3 --out out/n5     # write out/n5/<pattern>_completion<k>.png / .pdf / .dot
    python scripts/completions.py --only parity            # one pattern
    python scripts/completions.py --a=-1,-1,-1,-1,1 --B=1,-1,1,1,1/1,1,-1,1,1/... --name mine   # your own (five rows)

Conventions as everywhere: b_ii > 0, B[i, j] < 0 means species j benefits species i.
"""
import argparse
import time
from pathlib import Path

from signrules import Ecology, complete, write_completion, write_certificate
from signrules.patterns import PATTERNS, parse_pattern


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=2, help="completions per pattern")
    ap.add_argument("--out", default="out/n5_completions", help="directory for the .dot / .png files")
    ap.add_argument("--only", default=None, help="run only patterns whose name contains this string")
    ap.add_argument("--rules", default=",".join(Ecology.ALL_RULES))
    ap.add_argument("--a", default=None, help="custom growth-rate signs, e.g. -1,-1,-1,-1,1")
    ap.add_argument("--B", default=None, help="custom B signs, rows separated by /, e.g. 1,-1,1,1,1/...")
    ap.add_argument("--name", default="custom")
    ap.add_argument("--initial-band", action="store_true", help="also draw the ~190 entry-level signs as a band")
    ap.add_argument("--edges", default=None, choices=["all", "band"],
                    help="draw every dependency edge, or only those inside a band (default: band)")
    args = ap.parse_args()
    use = tuple(r for r in args.rules.split(",") if r)
    patterns = dict(PATTERNS)
    if args.a and args.B:
        patterns = {args.name: parse_pattern(args.a, args.B)}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for name, (a, B) in patterns.items():
        if args.only and args.only not in name:
            continue
        eco = Ecology(a, B, name)
        t = time.time()
        res = complete(eco, use=use, max_completions=args.max)
        print(f"{name:26s} a={list(map(int, a))}  {res.summary()}  [{time.time() - t:.1f}s]")
        if res.root_conflict:
            write_certificate(eco, res, str(out / f"{name}_certificate"))
            print(f"      wrote {out / name}_certificate.png (+ .pdf, .dot)")
        for k, d in enumerate(res.derivations):
            branches = "  ".join(f"p{''.join(map(str, S))}={'+' if s > 0 else '-'}" for S, s in d.branches)
            print(f"    completion {k}: {len(d.branches)} branches, {len(d.log)} rule steps, {len(d.gp_log)} GP steps")
            print(f"      branches: {branches}")
            stage_sizes = [len(l) + len(g) for _, l, g in d.stages()]
            print(f"      signs forced per stage (root, then per branch): {stage_sizes}")
            stem = out / f"{name}_completion{k}"
            write_completion(eco, res, str(stem), index=k, initial_band=args.initial_band or None, edges=args.edges)
            print(f"      wrote {stem}.png (+ .pdf, .dot)")


if __name__ == "__main__":
    main()
