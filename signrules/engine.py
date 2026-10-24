"""
Trigger-driven propagation engine for ``search.complete``.

The sweep engine (``Ecology.apply_rules`` + ``chirotope.propagate``) re-checks every
3-term relation and reruns every rule whenever one sign changes.  Here each rule is
split into small *propagators* — one per 3-term relation, signed sum (R1, R2, R5),
circuit (R7) or Schur block (R4, R6) — and every Plücker basis knows the
propagators that read or write it.  A propagator is re-run only when one of its
bases changes; between branch choices only the propagators touching the branched
basis are woken.  All propagators are monotone (a stronger table yields stronger
conclusions), so the fixed point — or the conflict — does not depend on the order
(Lagerkvist–Schulte 2009), and the completions are exactly those of the sweep.

Propagators sit in priority levels, cheapest first (Schulte–Stuckey 2008): the
queue of a level is drained only when all cheaper levels are empty, so the
expensive Schur-block rules see a table already closed under the cheap ones.

    level 0  3-term Grassmann–Plücker relations
    level 1  R1 (det B > 0, coefficient sums), R2n (Cramer on the full system), R8 (once)
    level 2  R5 Laplace expansions
    level 3  R2 Cramer on proper subsystems, R7 circuits
    level 4  R3 (once, at the root), R6 Schur blocks
    level 5  R4 Z-blocks

Tables are flat lists indexed by basis id (plucker_index.py); every propagator is
compiled to ids, and bases appear only in the logs and conflicts it writes.  The
3-term propagators and their watcher lists depend on n only and are cached.

``Engine.branch_basis`` picks the open basis for the next branch: the one read by
the most propagators, whose choice wakes the most relations and rules at once
(most-constrained first, Finschi–Fukuda).

``root_only`` names rules whose propagators run in the root fixpoint but are never woken
below it (an ablation, scripts/rule_ablation.py).  They stay registered, so the watcher
counts, hence ``order`` and ``branch_basis``, are those of the full rule set: only the
pruning below the root can differ.  ``idle`` names rules that are registered the same way
but never run at all, root included: ``search.complete`` idles R4 and R6 by default, since
on a table closed under the other rules they force nothing (docs/SIGNRULES.md, "R4 and R6
are subsumed"), while their watchers keep the branching order that was tuned with them
(without them the n = 5, 6 screens take 15–21 % more nodes).  ``stats=True`` counts, per
rule and separately at the root and below it, the runs, the signs set and the conflicts
raised (``rule_stats()``); the 3-term relations are not counted.
"""

from collections import deque
from functools import lru_cache
from typing import Dict, List, Optional, Sequence

from .chirotope import build_relation_ids
from .plucker_index import UNKNOWN, index
from .rules import Conflict, Ecology

GP_LEVEL = 0
LEVELS = 6


@lru_cache(maxsize=None)
def _gp_propagators(n: int):
    """(run, bases) for every 3-term relation of Gr(n, 2n); they do not depend on the ecology."""
    bases = index(n).bases
    return [(Engine._gp(rel, bases), tuple({b for t in rel for b in t[:2]})) for rel in build_relation_ids(n)]


@lru_cache(maxsize=None)
def _gp_watchers(n: int):
    """watchers[id] = the 3-term propagators (pids 0..) reading basis id."""
    watchers = [[] for _ in range(index(n).m)]
    for pid, (_, bases) in enumerate(_gp_propagators(n)):
        for S in bases:
            watchers[S].append(pid)
    return watchers


class Engine:
    """Compiled propagators of one ecology and one rule selection."""

    def __init__(self, eco: Ecology, use: Sequence[str] = Ecology.ALL_RULES,
                 root_only: Sequence[str] = (), stats: bool = False, idle: Sequence[str] = ()):
        self.eco = eco
        self.use = tuple(use)
        self.root_only_rules = frozenset(root_only)
        unknown = self.root_only_rules - set(self.use)
        if unknown:
            raise ValueError(f"root_only names rules not in use: {sorted(unknown)}")
        unknown = set(idle) - set(Ecology.RULES)
        if unknown:
            raise ValueError(f"idle names unknown rules: {sorted(unknown)}")
        self.idle_rules = frozenset(idle) & set(self.use)
        gp = _gp_propagators(eco.n)
        self.level: List[int] = [GP_LEVEL] * len(gp)
        self.idempotent = bytearray(b"\x01" * len(gp))   # 1: rerunning right after it fired cannot set anything new
        self.root_only = bytearray(len(gp))               # 1: never woken below the root
        self.idle = bytearray(len(gp))                    # 1: registered (watchers, branching order) but never run
        self.rule_of: List[str] = ["GP"] * len(gp)        # pid -> rule name
        self.run: List = [run for run, _ in gp]             # run(chi, log, gp_log) -> list of newly set bases
        self.watchers: List[List[int]] = [list(w) for w in _gp_watchers(eco.n)]
        self.stats: Optional[Dict[str, Dict[str, List[int]]]] = {} if stats else None
        self._at_root = True
        for name in self.use:
            self._compiling = name
            getattr(self, "_compile_" + name)()
        self._compiling = None
        self._skip_below = bytearray(a | b for a, b in zip(self.root_only, self.idle))
        # branching order: most-watched bases first (ties: table order)
        self.order = sorted(range(eco.idx.m), key=lambda S: -len(self.watchers[S]))

    # ------------------------------------------------------------ compilation

    def _add(self, level: int, run, bases):
        pid = len(self.run)
        self.level.append(level)
        self.idempotent.append(level <= 3)     # GP, sums and circuits; a Schur block may need a second pass
        self.root_only.append(self._compiling in self.root_only_rules)
        self.idle.append(self._compiling in self.idle_rules)
        self.rule_of.append(self._compiling)
        self.run.append(run)
        for S in bases:
            self.watchers[S].append(pid)

    def _rule(self, fn, need=()):
        """
        Wrap a rule call ``fn(chi, log)`` so that it reports the bases it set; it is skipped
        while any basis in ``need`` (inputs without which it cannot act) is open.
        """
        eco = self.eco
        need = tuple(need)

        def run(chi, log, gp_log):
            for S in need:
                if chi[S] == UNKNOWN:
                    return ()
            eco._trail = trail = []
            try:
                fn(chi, log)
            finally:
                eco._trail = None
            return trail
        if self.stats is None:
            return run
        counts = self.stats.setdefault(self._compiling, {"root": [0, 0, 0], "below": [0, 0, 0]})
        root, below = counts["root"], counts["below"]

        def counted(chi, log, gp_log):
            c = root if self._at_root else below
            c[0] += 1
            try:
                trail = run(chi, log, gp_log)
            except Conflict:
                c[2] += 1
                raise
            c[1] += len(trail)
            return trail
        return counted

    def rule_stats(self) -> Optional[Dict[str, Dict[str, Dict[str, int]]]]:
        """{rule: {"root" | "below": {"runs", "signs", "conflicts"}}} when built with stats=True."""
        if self.stats is None:
            return None
        return {name: {where: dict(zip(("runs", "signs", "conflicts"), c)) for where, c in counts.items()}
                for name, counts in self.stats.items()}

    def _keys_bases(self, keys):
        mk = self.eco._minor_key
        return {mk[k][0] for k in keys}

    def _sums(self, level, sums):
        eco = self.eco
        for terms, target, target_key, why in sums:
            keys = [k for c, k in terms if k is not None] + ([target_key] if target_key else [])
            iterms, ikeys = eco.compile_terms(terms)
            if target is None:
                ti, tc = eco._minor_key[target_key]

                def fn(chi, log, iterms=iterms, ikeys=ikeys, target_key=target_key, why=why, ti=ti, tc=tc):
                    eco.enforce_sum_ids(chi, iterms, ikeys, tc * chi[ti], target_key, log, why)
            else:
                def fn(chi, log, iterms=iterms, ikeys=ikeys, target=target, target_key=target_key, why=why):
                    eco.enforce_sum_ids(chi, iterms, ikeys, target, target_key, log, why)
            self._add(level, self._rule(fn), self._keys_bases(keys))

    def _compile_R1(self):
        eco = self.eco
        self._add(1, self._rule(lambda chi, log: eco.set_minor(chi, eco.species, eco.species, +1, log,
                                                               "R1 det B > 0")), ())
        self._sums(1, eco.stability_sums())

    def _compile_R2n(self):
        self._sums(1, self.eco.cramer_sums(full_system_only=True))

    def _compile_R2(self):
        self._sums(1, self.eco.cramer_sums(sizes=[self.eco.n]))
        self._sums(3, self.eco.cramer_sums(sizes=range(1, self.eco.n)))

    def _compile_R5(self):
        self._sums(2, self.eco.laplace_sums())

    def _compile_R7(self):
        eco = self.eco
        for circuit in eco._circuits:
            self._add(3, self._rule(lambda chi, log, circuit=circuit: eco.tope_circuit(chi, log, circuit)),
                      self._keys_bases([k for _, k in circuit[0]]))

    def _compile_R8(self):
        # depends on the entries of B and on sign a only: run once at the root
        self._add(1, self._rule(self.eco.rule_real_eigenvalue), ())

    def _compile_R3(self):
        # the signs R3 forces depend on (a, B) only: run once at the root, never woken again
        self._add(4, self._rule(self.eco.rule_balanced), ())

    def _compile_R6(self):
        eco = self.eco
        for P, Q in eco.block_pairs(min_P=1):
            need = self._keys_bases([(P, P)])       # cheap pre-check: nothing to do while det B_PP is open
            self._add(4, self._rule(lambda chi, log, P=P, Q=Q: eco.schur(chi, log, P, Q), need),
                      eco.block_bases(P, Q))

    def _compile_R4(self):
        eco = self.eco
        for P, Q in eco.block_pairs(min_P=0):
            # R4 needs det B_PP and every entry of the Schur block
            need = self._keys_bases(([(P, P)] if P else []) +
                                    [(tuple(sorted(P + (i,))), tuple(sorted(P + (j,)))) for i in Q for j in Q])
            self._add(5, self._rule(lambda chi, log, P=P, Q=Q: eco.z_block(chi, log, P, Q), need),
                      eco.block_bases(P, Q))

    @staticmethod
    def _gp(rel, bases):
        """One 3-term relation (ids): two known terms of one sign force the third to the other sign."""
        (a1, a2, ca), (b1, b2, cb), (c1, c2, cc) = rel
        ids = (a1, a2, b1, b2, c1, c2)

        def run(chi, log, gp_log):
            ta = ca * chi[a1] * chi[a2]
            tb = cb * chi[b1] * chi[b2]
            tc = cc * chi[c1] * chi[c2]
            if ta and tb and tc:
                if ta == tb == tc:
                    raise Conflict("3-term Grassmann–Plücker relation", [bases[b] for b in ids])
                return ()
            if ta == 0:
                if tb == 0 or tc == 0 or tb != tc:
                    return ()
                k1, k2, coef, forced = a1, a2, ca, -tb
            elif tb == 0:
                if tc == 0 or ta != tc:
                    return ()
                k1, k2, coef, forced = b1, b2, cb, -ta
            else:
                if ta != tb:
                    return ()
                k1, k2, coef, forced = c1, c2, cc, -ta
            v1, v2 = chi[k1], chi[k2]
            if v1 == 0 and v2 == 0:
                return ()
            if v1 == 0:
                new = k1
                chi[k1] = forced * coef * v2
            else:
                new = k2
                chi[k2] = forced * coef * v1
            if gp_log is not None:
                gp_log.append((bases[new], chi[new], [bases[b] for b in ids if b != new]))
            return (new,)
        return run

    # ------------------------------------------------------------ propagation

    def fixpoint(self, chi, log: list, gp_log: Optional[list] = None, changed=None):
        """
        Propagate in place to the fixed point.  ``changed`` lists the basis ids set since the
        table was last at a fixed point (a branch choice); None runs every propagator.
        Raises Conflict.
        """
        queues = [deque() for _ in range(LEVELS)]
        queued = bytearray(len(self.run))
        level, watchers, run, idem = self.level, self.watchers, self.run, self.idempotent
        at_root = changed is None
        self._at_root = at_root
        # idle rules are never queued; root_only rules are not queued below the root
        skip = self.idle if at_root else self._skip_below
        if at_root:
            for pid in range(len(run)):
                if not skip[pid]:
                    queues[level[pid]].append(pid)
                    queued[pid] = 1
        else:
            for S in changed:
                for pid in watchers[S]:
                    if not queued[pid] and not skip[pid]:
                        queued[pid] = 1
                        queues[level[pid]].append(pid)
        lv = 0
        while lv < LEVELS:
            q = queues[lv]
            if not q:
                lv += 1
                continue
            pid = q.popleft()
            queued[pid] = 0
            new = run[pid](chi, log, gp_log)
            for S in new:
                for w in watchers[S]:
                    if not queued[w] and not skip[w] and (w != pid or not idem[pid]):
                        queued[w] = 1
                        queues[level[w]].append(w)
                        if level[w] < lv:
                            lv = level[w]
        return chi

    # ------------------------------------------------------------ branching

    def branch_basis(self, chi, unknown: List = None):
        """
        The open basis (id) that the most propagators read: its choice wakes the most relations
        and rules, so conflicts surface near the root (most-constrained first).  On the
        eleven n = 4 and twenty-one n = 5 kills with R1, R2n, R5, R6 this cuts the pruned
        part of the tree about sevenfold against first-open branching; a dynamic score
        (relations with one other open basis) prunes no better and costs a scan per node.
        """
        for S in self.order:
            if chi[S] == UNKNOWN:
                return S
        return None
