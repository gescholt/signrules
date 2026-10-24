"""
Replay the impossibility certificates of the extension screens independently of the search.

For every pattern killed at the root in signrules/data/n<n>/extension_screen.jsonl (n = 5, 6), or for the
eleven four-species patterns (n = 4), the root propagation is re-run, its certificate (the backward
slice of the log ending in the conflict) is extracted with ``certificate_steps``, and every step is
re-derived by ``signrules/replay.py`` from the pattern and the signs established before it.  The
script prints one line per pattern (steps replayed, rules used, the rule reporting the conflict) and a
summary; a non-zero exit code means some certificate did not replay.

    python scripts/replay_certificates.py 4
    python scripts/replay_certificates.py 5 --workers 8
    python scripts/replay_certificates.py 6 --workers 12 --tamper     # also check that tampered steps fail
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


def patterns(n):
    if n == 4:
        from signrules.patterns import ELEVEN
        return [(name, a, np.array(B)) for name, (a, B) in ELEVEN.items()]
    recs = [json.loads(l) for l in open(DATA_DIR / f"n{n}" / "extension_screen.jsonl") if l.strip()]
    return [(str(r["code"]), r["a"], np.array(r["B"])) for r in recs if r["status"] == "killed_root"]


def _job(args):
    name, a, B, tamper = args
    from signrules import Ecology, complete, replay_certificate
    from signrules.certificate import certificate_steps
    eco = Ecology(a, B, name)
    t0 = time.time()
    res = complete(eco, max_nodes=0, max_completions=1, keep_derivations=False)
    if res.root_conflict is None:
        return dict(name=name, ok=False, message="not killed at the root")
    steps, prem, _ = certificate_steps(eco, res, backward_slice=True)
    rep = replay_certificate(a, B, steps, res.root_conflict, prem)
    out = dict(name=name, ok=rep.ok, message=rep.message, steps=rep.steps, rules=rep.rules,
               conflict=res.root_conflict.split(":")[0].split(" ")[0], seconds=round(time.time() - t0, 2))
    if tamper and steps:
        # flipping the sign of any single step must make the replay fail (the flipped step itself, a
        # later step that depended on it, or the conflict)
        bad = []
        for k in range(len(steps)):
            S, sign, why, p = steps[k]
            tampered = list(steps)
            tampered[k] = (S, -sign, why, p)
            if replay_certificate(a, B, tampered, res.root_conflict, prem).ok:
                bad.append(k)
        out["tamper_ok"] = not bad
        out["tamper_survivors"] = bad
        if bad:
            out["ok"] = False
            out["message"] += f"; tampering survived at steps {bad}"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tamper", action="store_true", help="also flip each step in turn and require failure")
    ap.add_argument("--out", default=None, help="write the per-pattern results as JSON lines")
    args = ap.parse_args()
    pats = patterns(args.n)
    with Pool(args.workers) as pool:
        results = pool.map(_job, [(name, a, B, args.tamper) for name, a, B in pats])
    rules_total = {}
    for r in results:
        flag = "ok  " if r["ok"] else "FAIL"
        rules = " ".join(f"{k}:{v}" for k, v in sorted(r.get("rules", {}).items()))
        print(f"{flag} {r['name']:>12s} {r.get('steps', 0):4d} steps  conflict {r.get('conflict', '?'):4s} [{rules}] {r['message']}")
        for k, v in r.get("rules", {}).items():
            rules_total[k] = rules_total.get(k, 0) + v
    n_ok = sum(r["ok"] for r in results)
    print(f"{n_ok}/{len(results)} certificates replayed; steps by rule: {dict(sorted(rules_total.items()))}")
    if args.out:
        with open(args.out, "w") as f:
            for r in results:
                f.write(json.dumps(r) + "\n")
    sys.exit(0 if n_ok == len(results) else 1)


if __name__ == "__main__":
    main()
