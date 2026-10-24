"""
Depth-first sign completion.

``complete(eco)`` starts from the signs determined by the entries of (a, B),
runs the 3-term Grassmann–Plücker propagation and the ecological sign rules
to a fixed point, and then branches on an undetermined Plücker sign,
repeating the fixed-point computation in every node.  By default the fixed
point is computed by ``engine.Engine`` (only the relations and rules whose
coordinates changed are re-run, cheapest first) and the branch goes to the
open coordinate with the highest conflict activity (ties and the start: the
one read by the most propagators); ``engine="sweep"`` and ``branching="first"``
restore the original full sweeps and table order.  ``symmetry=True`` prunes branches that a
relabelling fixing the pattern maps onto branches already explored (see ``symmetry.py``).

Tables are flat lists indexed by basis id (``eco.idx.bases[i]``, ``eco.idx.id[S]``,
``plucker_index.py``); the search works with ids and writes bases into everything it
records (branches, logs, tree nodes, progress events).

Outcomes
--------
* ``result.root_conflict`` is set: the pattern is impossible, without any
  branching; ``result.root_log`` lists every sign that was forced on the way,
  each with the rule that forced it (the certificate).
* ``result.completions`` is empty after a full search: impossible as well.
* ``result.completions`` is non-empty: each entry is a complete sign table
  consistent with all rules.  This is *not* a proof of possibility — a
  complete table need not be realisable by real parameters.  The matching
  entry of ``result.derivations`` records how the table was reached: the
  branch choices and every sign forced along the path, with its rule
  (rendered by ``certificate.completion_dot``).

``use`` selects the rules by name (see ``Ecology.RULES``); ``("R1", "R2n")``
reproduces the feasibility + weak-stability algorithm of the Grassmannian
paper, and ``()`` runs the chirotope propagation alone.
"""

import random
import time
from typing import Dict, List, Optional, Sequence, Tuple

from .plucker_index import UNKNOWN, Basis
from .chirotope import build_relation_ids, propagate
from .engine import Engine
from .rules import Ecology, Conflict
from .symmetry import stabilizer_actions


class Derivation:
    """
    How one completion was reached.

    * ``branches``: the (basis, sign) choices made on the path, root first;
    * ``log``: (why, P, Q, sign, premises) for every sign forced by a rule on the path
      (root steps first), same format as ``SearchResult.root_log``;
    * ``gp_log``: (basis, sign, premises) for every sign forced by a 3-term relation;
    * ``marks``: for each branch, the lengths of ``log`` and ``gp_log`` when it was made,
      so ``stages()`` can split the derivation into root propagation and one stage per branch.
    """

    def __init__(self, branches: List[Tuple[Basis, int]], log: List, gp_log: List,
                 marks: Optional[List[Tuple[int, int]]] = None):
        self.branches = branches
        self.log = log
        self.gp_log = gp_log
        self.marks = marks if marks is not None else [(len(log), len(gp_log))] * len(branches)

    def stages(self):
        """Yield (branch or None, rule-log slice, gp-log slice): the root stage, then one per branch."""
        bounds = [(0, 0)] + list(self.marks) + [(len(self.log), len(self.gp_log))]
        for k, ((l0, g0), (l1, g1)) in enumerate(zip(bounds, bounds[1:])):
            yield (None if k == 0 else self.branches[k - 1]), self.log[l0:l1], self.gp_log[g0:g1]

    def __repr__(self):
        return f"Derivation({len(self.branches)} branches, {len(self.log)} rule steps, {len(self.gp_log)} GP steps)"


class SearchResult:
    """Outcome of ``complete``."""

    def __init__(self):
        self.completions: List[List[int]] = []     # complete sign tables found (flat, by basis id)
        self.derivations: List[Derivation] = []    # how each completion was reached (same order)
        self.tree: List[dict] = []                 # search-tree nodes when record_tree=True (see complete())
        self.nodes = 0                             # search-tree nodes visited
        self.conflicts = 0                         # branches pruned below the root
        self.root_log: List = []                   # (why, P, Q, sign, premises) forced by rules at the root
        self.root_gp_log: List = []                # (basis, sign, premises) forced by 3-term relations at the root
        self.root_conflict: Optional[str] = None   # reason, when killed at the root
        self.root_conflict_premises: List = []     # minors (P, Q) or bases whose signs the contradiction used
        self.root_unknowns: Optional[int] = None   # undetermined signs after root propagation
        self.probes: List[Tuple[int, Basis, int]] = []  # (depth, basis, sign) fixed by failed-literal probing
        self.symmetry: List = []                   # stabilizer actions used (symmetry=True), identity excluded
        self.symmetry_forced = 0                   # signs forced by symmetry in second branches
        self.symmetry_cut = 0                      # nodes cut because a relabelling maps them into a searched branch
        self.stopped = False                       # search cut short by max_completions / max_nodes
        self.seconds = 0.0
        self.rule_stats: Optional[dict] = None     # {rule: {"root"|"below": {runs, signs, conflicts}}} (stats=True)
        self.assumed: Dict[int, int] = {}          # the extra hypotheses of complete(assume=...), by basis id

    def summary(self) -> str:
        if self.root_conflict:
            return f"impossible at the root ({self.root_conflict})"
        tail = " (search stopped early)" if self.stopped else ""
        return f"{len(self.completions)} completions, {self.nodes} nodes, {self.conflicts} pruned{tail}"


SUBSUMED = {"R6": frozenset({"R5", "R7"}), "R4": frozenset({"R2", "R3", "R5", "R7", "R8"})}


def subsumed_rules(use: Sequence[str]) -> tuple:
    """The rules of ``use`` that the other rules of ``use`` subsume (R6 by R5 and R7, R4 by R2, R3, R5, R7, R8)."""
    present = set(use)
    return tuple(r for r in ("R4", "R6") if r in present and SUBSUMED[r] <= present)


def complete(eco: Ecology, use: Sequence[str] = Ecology.ALL_RULES,
             max_completions: Optional[int] = None, keep_derivations: bool = True,
             record_tree: bool = False, max_nodes: Optional[int] = None,
             progress=None, engine: str = "trigger", branching: str = "conflict",
             probe=None, symmetry: bool = False, decide: Optional[Sequence[int]] = None,
             seed: Optional[int] = None, root_only: Sequence[str] = (), stats: bool = False,
             assume: Optional[Dict[int, int]] = None, idle="auto") -> SearchResult:
    """
    Complete the sign table of ``eco`` by propagation and branching.

    Parameters
    ----------
    eco : Ecology
    use : rule names to apply at every node (default: all)
    idle : names among ``use`` that are never run.  The default ``"auto"`` idles the subsumed
        rules whose subsuming rules are in ``use``: R6 when R5 and R7 are, R4 when R2, R3, R5, R7
        and R8 are (on a table closed under those rules they force nothing: docs/SIGNRULES.md,
        "R4 and R6 are subsumed").  With the trigger engine their propagators stay registered, so
        the branching order, the trees and the certificates are those of the full rule set;
        ``idle=()`` runs them.
    root_only : names among ``use`` to apply at the root only (trigger engine; an ablation, see
        ``engine.Engine``).  Their propagators stay registered, so the branching order is unchanged;
        the root certificate is unchanged; the pruning below the root may be weaker.
    stats : count the runs, signs set and conflicts raised per rule, at the root and below
        (``result.rule_stats``; trigger engine)
    max_completions : stop after this many completions (None = exhaustive)
    keep_derivations : record branch choices and forced signs for every completion
        (``result.derivations``); switch off to save memory in very large searches
    record_tree : record every node of the binary search tree in ``result.tree`` as a dict with
        ``id``, ``parent``, ``depth``, ``basis`` and ``sign`` (the choice that led here),
        ``forced`` (signs forced by propagation at this node), ``unknown`` (signs still open),
        ``status`` ("root", "branch", "completion", "conflict", or "stopped" when the completion
        limit ended the search there), ``branch_on`` (the basis this
        node branches on), ``completion`` (index) and ``reason`` (for conflicts).
        ``certificate.search_tree_text`` / ``search_tree_dot`` draw it.
    max_nodes : stop the search after visiting this many nodes (``result.stopped`` is then True and,
        with ``record_tree``, the unsearched nodes are marked "stopped")
    progress : optional callable, called with one dict per event of the search: ``kind`` ("root",
        "node", "forced" (propagation after a branch choice), "completion", "conflict", "stopped"), ``depth``, ``path`` (the branch choices),
        ``unknown`` (open signs), ``log`` / ``gp_log`` (the signs forced by this event) and, for
        conflicts, ``reason``.  ``livesearch.LiveSearch`` renders these as a terminal dashboard.
    engine : "trigger" (``engine.Engine``: re-run only the propagators whose bases changed, cheap
        rules first) or "sweep" (rerun every relation and rule until nothing changes).  Both reach
        the same fixed points, hence the same verdicts and completions; the rule that reports a
        conflict, and the order of the logs, may differ.
    branching : "conflict" (default), "constrained" (branch on the open basis read by the most
        propagators, ``Engine.branch_basis``), "first" (the first open basis); "conflict" is activity-based
        as VSIDS / dom-wdeg: every conflict bumps the bases among its premises, with a growing
        increment so recent conflicts weigh more, and the open basis of highest activity is chosen;
        ties and the start follow "constrained".  Changes the tree, not the set of completions.
    probe : failed-literal probing (lookahead) after propagation: every open basis is set to each
        sign in turn on a copy of the table and propagated; if one sign conflicts the other is
        forced, if both conflict the node is a conflict; repeated until no basis is forced.
        ``"root"`` (or True) probes at the root only, ``"all"`` at every node, None not at all.
        Sound and order-independent like the rules, so the completions do not change; the probed
        signs are listed in ``result.probes`` and are not part of the rule logs (a kill that
        probing finds is a case split of depth one, not a single derivation).
    symmetry : prune with the stabilizer G of the pattern (relabellings fixing sign a and sign B,
        ``symmetry.py``).  At a node whose table T is fixed by the subgroup H of G, once the branch
        chi[S] = + is searched, every literal g·(chi[S] = +), g in H, has an isomorphic subtree, so
        the branch chi[S] = - also sets chi[gS] to the opposite of that image (logged with rule
        "SYM"); when some g maps chi[S] = + to chi[S] = -, the second branch is cut whole.  More
        generally (SBDS, Gent–Smith 2000) every node below that branch whose table satisfies g(A),
        A the decisions above S, for any g in G, gets chi[gS] set opposite to the image of +.  Verdicts
        are unchanged, and every completion is the image under G of one in ``result.completions``
        (``symmetry.expand_orbits(result.completions, result.symmetry)`` recovers them all), but the
        list may hold a few orbits more than once.  ``max_completions`` then counts these
        representatives.
    decide : basis ids (``eco.idx.id[S]``) to decide; the other signs are only propagated.  A node where
        every id in ``decide`` is known is a leaf and its table (with UNKNOWN elsewhere) is appended to
        ``result.completions``: the *coarse cells*, the distinct projections of the completions onto
        these coordinates that propagation allows (a coarse cell may still have no completion; every
        completion projects into one).  Branching is restricted to open ids in ``decide``.
    seed : when given, the branch to try first at every node (+ or -) is drawn at random with this
        seed, so ``max_completions=1`` returns a random completion (or coarse cell) at the cost of one
        root-to-leaf path.  Verdicts and the full set of completions are unchanged, the tree is not.
        Not with ``symmetry``.
    assume : extra hypotheses {basis id: ±1} written into the initial table before root propagation,
        so the search runs over the completions with these signs (``result.assumed``).  A sign that
        contradicts the initial table is a root conflict.  With ``max_completions=1`` this decides
        whether any completion has the assumed signs: a completion found says yes; an exhausted
        search (``stopped`` False, no completion) says no, exactly (``scripts/sign_correlations_exact.py``).
    """
    result = SearchResult()
    t0 = time.time()
    bases = eco.idx.bases
    decide_set = None if decide is None else set(decide)
    if seed is not None and symmetry:
        raise ValueError("seed (random branch order) cannot be combined with symmetry (its nogoods assume + first)")
    value_rng = None if seed is None else random.Random(seed)
    if engine != "trigger" and (root_only or stats):
        raise ValueError("root_only and stats need engine='trigger'")
    relations = build_relation_ids(eco.n)
    if isinstance(idle, str):
        if idle != "auto":
            raise ValueError("idle is 'auto' or a sequence of rule names")
        idle = subsumed_rules(use)
    trigger = (Engine(eco, use, root_only=root_only, stats=stats, idle=idle)
               if engine == "trigger" or branching in ("constrained", "conflict") else None)
    sweep_use = tuple(r for r in use if r not in set(idle))
    activity: List[float] = [0.0] * eco.idx.m
    any_activity = [False]
    if symmetry:
        result.symmetry = stabilizer_actions(eco)
    bump = [1.0]

    def on_conflict(c):
        if branching != "conflict":
            return
        for k in c.premises:
            S = eco.idx.id.get(k) if isinstance(k[0], int) else eco._minor_key.get(k, (None,))[0]
            if S is not None:
                activity[S] += bump[0]
                any_activity[0] = True
        bump[0] *= 1.05

    def pick(table, unknown):
        if branching == "first":
            return unknown[0]
        if decide_set is None:
            S = trigger.branch_basis(table, unknown)
        else:                                       # most-constrained among the ids to decide
            S = next((T for T in trigger.order if T in decide_set and table[T] == UNKNOWN), unknown[0])
        if branching == "conflict" and any_activity[0]:
            best = max(unknown, key=activity.__getitem__)
            if activity[best] > 0:
                return best
        return S

    def fixpoint(table, log, gp_log=None, changed=None):
        """Propagate to a fixed point; ``changed``: the basis ids set since the last one."""
        if engine == "trigger":
            return trigger.fixpoint(table, log, gp_log, changed)
        while True:
            consistent, _ = propagate(table, relations, gp_log, bases)
            if not consistent:
                rel = propagate.last_conflict
                raise Conflict("3-term Grassmann–Plücker relation",
                               [bases[b] for t in rel for b in t[:2]] if rel else [])
            before = table.count(UNKNOWN)
            eco.apply_rules(table, log, sweep_use)
            if table.count(UNKNOWN) == before:
                return table

    if probe is True:
        probe = "root"

    def lookahead(table, log, gp_log, depth):
        """Failed-literal probing in place; raises Conflict when both signs of a basis conflict."""
        while True:
            forced = False
            for S in (trigger.order if trigger else range(len(table))):
                if table[S] != UNKNOWN:
                    continue
                ok = []
                for sign in (1, -1):
                    trial = table[:]
                    trial[S] = sign
                    try:
                        fixpoint(trial, [], None, [S])
                        ok.append(sign)
                    except Conflict:
                        pass
                if not ok:
                    raise Conflict(f"failed-literal probing: both signs of basis {bases[S]} conflict", [bases[S]])
                if len(ok) == 1:
                    table[S] = ok[0]
                    result.probes.append((depth, bases[S], ok[0]))
                    fixpoint(table, log, gp_log, [S])
                    forced = True
            if not forced:
                return table

    chi = eco.initial_table()
    n_initial_unknown = chi.count(UNKNOWN)
    try:
        if assume:
            result.assumed = {int(S): int(s) for S, s in assume.items()}
            for S, s in result.assumed.items():
                if s not in (1, -1):
                    raise ValueError("assume maps basis ids to ±1")
                if chi[S] == -s:
                    raise Conflict(f"assumption chi[{bases[S]}] = {s:+d} contradicts the initial table", [bases[S]])
                chi[S] = s
        fixpoint(chi, result.root_log, result.root_gp_log)
        if probe:
            lookahead(chi, result.root_log, result.root_gp_log, 0)
    except Conflict as reason:
        if progress:
            progress(dict(kind="conflict", depth=0, path=[], unknown=None, reason=str(reason),
                          log=result.root_log, gp_log=result.root_gp_log))
        result.root_conflict = str(reason)
        result.root_conflict_premises = list(reason.premises)
        result.seconds = time.time() - t0
        result.rule_stats = trigger.rule_stats() if stats else None
        if record_tree:
            result.tree.append(dict(id=0, parent=None, depth=0, basis=None, sign=None, status="conflict",
                                    forced=n_initial_unknown - chi.count(UNKNOWN),
                                    unknown=None, branch_on=None, reason=str(reason).split(":")[0],
                                    reason_full=str(reason), premises=list(reason.premises),
                                    log=list(result.root_log), gp_log=list(result.root_gp_log)))
        return result
    result.root_unknowns = chi.count(UNKNOWN)

    def node(**kw):
        kw["id"] = len(result.tree)
        result.tree.append(kw)
        return kw["id"]

    def dfs(table, branches, branch_ids, log, gp_log, marks, node_id, nogoods=()):
        """``branches``: the (basis, sign) choices on the path; ``branch_ids``: the same with ids."""
        result.nodes += 1
        if (max_completions is not None and len(result.completions) >= max_completions) or \
                (max_nodes is not None and result.nodes > max_nodes):
            result.stopped = True
            if progress:
                progress(dict(kind="stopped", depth=len(branches), path=branches, unknown=None, log=[], gp_log=[]))
            if record_tree:
                result.tree[node_id]["status"] = "stopped"      # limit reached, subtree not searched
            return
        unknown = [T for T, v in enumerate(table) if v == UNKNOWN]
        if decide_set is not None:
            unknown = [T for T in unknown if T in decide_set]
        if progress:
            progress(dict(kind="completion" if not unknown else "node", depth=len(branches), path=branches,
                          unknown=len(unknown), log=[], gp_log=[]))
        if not unknown:
            result.completions.append(table[:])
            if keep_derivations:
                result.derivations.append(Derivation(list(branches), list(log), list(gp_log), list(marks)))
            if record_tree:
                result.tree[node_id]["status"] = "completion"
                result.tree[node_id]["completion"] = len(result.completions) - 1
            return
        S = pick(table, unknown)
        if record_tree:
            result.tree[node_id]["branch_on"] = bases[S]
        first = 1 if value_rng is None else value_rng.choice((1, -1))
        for sign in (first, -first):
            # the branch chi[S] = first is searched once we get to -first: below, g(branches) => not g(chi[S] = first)
            child_nogoods = nogoods + ((tuple(branch_ids), S),) if sign == -first and result.symmetry else nogoods
            child = table[:]
            child[S] = sign
            track = keep_derivations or record_tree or progress is not None
            child_log = list(log) if track else []
            child_gp = list(gp_log) if track else None
            try:
                fixpoint(child, child_log, child_gp, [S])
                if probe == "all":
                    lookahead(child, child_log, child_gp, len(branches) + 1)
                if child_nogoods:
                    symmetric_closure(child, child_log if track else None, child_gp, child_nogoods)
            except Conflict as c:
                result.conflicts += 1
                on_conflict(c)
                if progress:
                    progress(dict(kind="conflict", depth=len(branches) + 1, path=branches + [(bases[S], sign)],
                                  unknown=None, reason=str(c), log=child_log[len(log):], gp_log=child_gp[len(gp_log):]))
                if record_tree:
                    node(parent=node_id, depth=len(branches) + 1, basis=bases[S], sign=sign, status="conflict",
                         forced=None, unknown=None, branch_on=None, reason=str(c).split(":")[0],
                         reason_full=str(c), premises=list(c.premises),
                         log=child_log[len(log):], gp_log=child_gp[len(gp_log):])
                continue
            if progress:
                progress(dict(kind="forced", depth=len(branches) + 1, path=branches + [(bases[S], sign)],
                              unknown=None, log=child_log[len(log):], gp_log=child_gp[len(gp_log):]))
            child_id = None
            if record_tree:
                left = child.count(UNKNOWN)
                child_id = node(parent=node_id, depth=len(branches) + 1, basis=bases[S], sign=sign,
                                status="branch", forced=len(unknown) - 1 - left, unknown=left, branch_on=None,
                                log=child_log[len(log):], gp_log=child_gp[len(gp_log):])
            dfs(child, branches + [(bases[S], sign)], branch_ids + [(S, sign)], child_log,
                child_gp if child_gp is not None else [],
                marks + [(len(log), len(gp_log))], child_id, child_nogoods)

    def symmetric_closure(child, log, gp_log, nogoods):
        """
        SBDS (Gent–Smith 2000).  A nogood (A, S) says that the branch chi[S] = + below the decisions
        A has been searched; for every g in G, a table satisfying g(A) must then violate g(chi[S] = +),
        since g maps its completions to completions already covered.  Forces those signs, propagates,
        and repeats until nothing is forced; raises Conflict when one is violated.
        """
        decomps = eco.idx.decomps
        while True:
            new = []
            for A, S in nogoods:
                for g in result.symmetry:
                    image, gsign = g.image, g.sign
                    if any(child[image[B]] != gsign[B] * s for B, s in A):
                        continue
                    T, v = image[S], -gsign[S]
                    if child[T] == v:
                        continue
                    _, P, Q, _ = decomps[T]
                    why = (f"SYM π={g.perm}: image of the searched branch p_{bases[S]} > 0 "
                           f"is p_{bases[T]} {'>' if v < 0 else '<'} 0")
                    premises = [decomps[image[B]][1:3] for B, _ in A]
                    if child[T] != UNKNOWN:
                        result.symmetry_cut += 1
                        raise Conflict(why + ", already covered", premises + [(P, Q)])
                    child[T] = v
                    new.append(T)
                    result.symmetry_forced += 1
                    if log is not None:
                        log.append((why, P, Q, eco._minor_key[(P, Q)][1] * v, premises))
            if not new:
                return
            fixpoint(child, log if log is not None else [], gp_log, new)

    root_id = None
    if record_tree:
        root_id = node(parent=None, depth=0, basis=None, sign=None, status="root",
                       forced=n_initial_unknown - result.root_unknowns, unknown=result.root_unknowns, branch_on=None,
                       log=list(result.root_log), gp_log=list(result.root_gp_log))
    if progress:
        progress(dict(kind="root", depth=0, path=[], unknown=result.root_unknowns,
                      log=result.root_log, gp_log=result.root_gp_log))
    dfs(chi, [], [], result.root_log, result.root_gp_log, [], root_id)
    result.seconds = time.time() - t0
    result.rule_stats = trigger.rule_stats() if stats else None
    return result
