"""
The eleven four-species ecologies conjectured impossible in Meng et al. (2025) and proved impossible in
Çelik–Isoldi–Portakal–Zucal, "Stable Coexistence in Ecologies and Games" (2026), run through the sign rules
(the patterns themselves are ``signrules.patterns.ELEVEN``).

    python scripts/impossible_eleven.py                 # all rules: every pattern dies at the root
    python scripts/impossible_eleven.py --rules R1,R2n  # feasibility + weak stability only:
                                                        # the counts quoted in the paper (64, 64, 1996, 1528, ...)
    python scripts/impossible_eleven.py --rules ""      # chirotope propagation alone
    python scripts/impossible_eleven.py --engine sweep  # the kill logs of docs/SIGNRULES.md §4 (R3/R4/R6)
    python scripts/impossible_eleven.py --only S-d --rules R1,R2n --dot out/   # derivation graphs:
                # certificate of every kill, first completion of every survivor (dot -Tpng out/*.dot)
"""
import argparse
from pathlib import Path

from signrules import Ecology, complete, write_certificate, write_completion, UNKNOWN
from signrules.patterns import ELEVEN


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", default=",".join(Ecology.ALL_RULES), help="comma-separated rule names")
    ap.add_argument("--only", default=None, help="run only the patterns whose name contains this string")
    ap.add_argument("--max", type=int, default=None, help="stop after this many completions")
    ap.add_argument("--engine", default="trigger", choices=["trigger", "sweep"],
                    help="propagation engine; 'sweep' gives the kill logs that follow the paper's proofs")
    ap.add_argument("--branching", default="conflict", choices=["conflict", "constrained", "first"])
    ap.add_argument("--dot", default=None, help="directory for Graphviz derivation graphs (certificate or first completion)")
    args = ap.parse_args()
    use = tuple(r for r in args.rules.split(",") if r)
    print(f"rules: {use}\n")
    for name, (a, B) in ELEVEN.items():
        if args.only and args.only not in name:
            continue
        eco = Ecology(a, B, name)
        res = complete(eco, use=use, max_completions=args.max, engine=args.engine, branching=args.branching)
        n_init = eco.initial_table().count(UNKNOWN)
        print(f"{name:45s} unknown signs: {n_init:2d} initially, "
              f"{'-' if res.root_unknowns is None else res.root_unknowns:>2} after root propagation | "
              f"{res.summary()} [{res.seconds:.1f}s]")
        if res.root_conflict is None and res.root_log:
            fired = sorted({w.split()[0] for w, *_ in res.root_log})
            print(f"{'':45s} rules that forced signs at the root: {fired}")
        if args.dot:
            out = Path(args.dot); out.mkdir(parents=True, exist_ok=True)
            stem = out / name.split()[0]
            if res.root_conflict:
                write_certificate(eco, res, f"{stem}_certificate")        # .png + .pdf + .dot
            elif res.derivations:
                write_completion(eco, res, f"{stem}_completion0", index=0)


if __name__ == "__main__":
    main()
