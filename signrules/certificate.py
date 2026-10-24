"""
Derivation graphs and search trees: which rule forced which Plücker sign.

Every sign forced during a search is logged with the rule that forced it and the signs it
used (``SearchResult.root_log`` / ``root_gp_log`` at the root, ``Derivation.log`` / ``gp_log``
along the path to a completion), and a contradiction records its premises.  This module draws
those records.  All text comes from ``labels.py`` (Plücker coordinates ``$p_S$``, rules as
``R4 · Z-block`` plus their parameters), the layout from Graphviz and the typesetting from
matplotlib (``gvrender.Diagram``), so the same drawing is written as PNG, PDF and a plain-text
``.dot`` file.

* ``certificate_diagram`` — the *certificate* of a kill at the root: the backward slice of the
  record, i.e. only the steps the contradiction depends on,

      signs read off the entries  →  [rule]  →  forced signs  →  …  →  ⊥

* ``completion_diagram`` — the derivation of one *completion*: one band per search stage (the
  root propagation, then each branch choice with the signs it forces);
* ``search_tree_diagram`` — the binary search tree: nodes are the coordinates branched on,
  edges the sign chosen (red +, blue −), leaves completions (green), contradictions (⊥) or
  unexplored subtrees (grey).  ``treeplot.py`` draws the same tree without Graphviz, for large
  trees; ``search_tree_text`` prints it.

Colour: red = positive, blue = negative (the fill of a sign node); the border says where the sign
came from — thin grey: read off the entries; black: forced by a rule; dashed violet: a branch
choice.  Rule boxes are white, the contradiction is black.

Usage::

    result = complete(eco)
    write_certificate(eco, result, "certificate")            # certificate.png / .pdf / .dot
    write_completion(eco, result, "completion", index=0)
    result = complete(eco, use=("R1", "R2n"), record_tree=True)
    write_search_tree(eco, result, "tree")
"""

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .plucker_index import UNKNOWN
from .rules import Ecology
from .search import SearchResult
from . import labels as L
from . import style
from .gvrender import Diagram

GP_RULE = "Grassmann–Plücker 3-term relation"


# --------------------------------------------------------------------------- the record as steps

def _steps(eco: Ecology, log: Sequence, gp_log: Sequence) -> Tuple[List[Tuple], callable]:
    """Both logs as (basis, sign, why, premise_bases); GP steps first, then rule steps."""
    bases = eco.idx.bases
    key_to_basis = {k: bases[i] for k, (i, _) in eco._minor_key.items()}

    def to_basis(k):
        return k if (isinstance(k, tuple) and k and isinstance(k[0], int)) else key_to_basis[tuple(map(tuple, k))]

    steps: List[Tuple] = []
    for S, sign, premises in gp_log:
        steps.append((S, sign, GP_RULE, [to_basis(p) for p in premises]))
    for why, P, Q, sign, premises in log:
        i, c = eco._minor_key[(tuple(P), tuple(Q))]
        steps.append((bases[i], c * sign, why, [to_basis(p) for p in premises]))
    # the two logs interleave in time, but premises always precede their consequences,
    # so any order in which each step follows its premises (a topological order) is fine.
    return steps, to_basis


def _backward_slice(steps: List[Tuple], targets: Iterable) -> List[Tuple]:
    """Only the steps that the signs of ``targets`` depend on."""
    produced = {S: i for i, (S, *_rest) in enumerate(steps)}
    needed, stack = set(), list(targets)
    while stack:
        S = stack.pop()
        if S in needed or S not in produced:
            continue
        needed.add(S)
        stack.extend(steps[produced[S]][3])
    return [st for st in steps if st[0] in needed]


def certificate_steps(eco: Ecology, result: SearchResult, backward_slice: bool = True):
    """Return (steps, conflict_premises, initial) with steps = [(basis, sign, why, premise_bases)]."""
    steps, to_basis = _steps(eco, result.root_log, result.root_gp_log)
    conflict_premises = [to_basis(p) for p in result.root_conflict_premises]
    if backward_slice:
        steps = _backward_slice(steps, conflict_premises)
    return steps, conflict_premises, eco.idx.to_dict(eco.initial_table())


# --------------------------------------------------------------------------- node styles

def _nid(S) -> str:
    return "s" + "_".join(map(str, S))


def _sign_node(d: Diagram, eco: Ecology, S, sign, origin: str, label_style: str, cluster=None,
               note: Optional[str] = None, size: Optional[float] = None) -> str:
    """A sign node: fill = sign, border = origin ('initial', 'forced', 'branch')."""
    lines: List = list(L.basis(eco, S, sign, label_style))
    if label_style == "both" and len(lines) > 1:
        lines[1] = (lines[1], (size or d.fontsize) * 0.75, style.GREY)
    if note:
        lines.append((note, (size or d.fontsize) * 0.7, style.GREY))
    border = dict(initial=(style.LIGHT, 0.9, "-"), forced=(style.INK, 1.3, "-"), branch=(style.CHOICE, 1.4, "--"))[origin]
    return d.node(_nid(S), lines, fill=style.sign_tint(sign) if sign != UNKNOWN else "white",
                  edge=border[0], lw=border[1], ls=border[2], size=size, cluster=cluster)


def _rule_node(d: Diagram, rid: str, lab: L.RuleLabel, cluster=None, with_clause: bool = False) -> str:
    fs = d.fontsize * 0.85
    lines: List = [(lab.headline(), fs, style.rule_color(lab.rule))]
    for extra in lab.lines(with_clause=with_clause)[1:]:
        lines.append((extra, fs * 0.9, style.GREY))
    return d.node(rid, lines, shape="box", fill="white", edge=style.RULE_EDGE, lw=0.8, size=fs, cluster=cluster,
                  pad=(6.0, 3.0))


def _bottom_node(d: Diagram, eco: Ecology, reason: str) -> str:
    lines = L.conflict(eco, reason)
    fs = d.fontsize * 0.85
    styled: List = [(lines[0], d.fontsize * 1.3, "white")] + [(t, fs, "white") for t in lines[1:]]
    rule = L.rule(reason, eco.n).rule
    return d.node("bottom", styled, fill=style.BOTTOM, edge=style.rule_color(rule), lw=2.2, shape="round", pad=(9.0, 5.0))


def _legend(d: Diagram, branch: bool = False) -> None:
    d.legend_item("read off the entries", fill="white", edge=style.LIGHT, lw=0.9)
    d.legend_item("forced by a rule", fill="white", edge=style.INK, lw=1.3)
    if branch:
        d.legend_item("branch choice", fill="white", edge=style.CHOICE, lw=1.4, ls="--")
    d.legend_item("positive", fill=style.POS_TINT, edge=style.LIGHT)
    d.legend_item("negative", fill=style.NEG_TINT, edge=style.LIGHT)


def _default_title(eco: Ecology, what: str) -> List:
    head = f"{what}: {eco.name}" if getattr(eco, "name", "") else what
    return [head, (L.pattern(eco), 9.5, style.GREY)]


# --------------------------------------------------------------------------- certificate DAG

def _flat_diagram(eco: Ecology, steps: List[Tuple], initial: Dict, branches: Sequence[Tuple] = (),
                  conflict: Optional[Tuple[str, List]] = None, title=None, label_style: str = "plucker",
                  legend: bool = True, per_row: int = 9) -> Diagram:
    """
    The DAG with explicit rule boxes:  premises → [rule] → consequences.  ``steps`` as from
    ``_steps``; ``branches`` the (basis, sign) choices; ``conflict`` = (reason, premise_bases).
    One rule box per application (same rule and premises), possibly with several outputs.
    """
    d = Diagram(title=title, nodesep=0.18, ranksep=0.42)
    produced = {S for S, *_ in steps}
    branch_sign = dict(branches)
    used = set()
    for S, sign, why, prem in steps:
        used.update(p for p in prem if p not in produced)
    if conflict:
        used.update(p for p in conflict[1] if p not in produced)
    for S in sorted(used):
        if S in branch_sign:
            _sign_node(d, eco, S, branch_sign[S], "branch", label_style, note="branch")
        else:
            _sign_node(d, eco, S, initial[S], "initial", label_style)
    # many initial signs feeding one rule: wrap them into rows instead of one very wide rank
    first_rule_prem: Dict[Tuple, List] = {}
    for S, sign, why, prem in steps:
        first_rule_prem.setdefault((L.rule(why, eco.n).key(), tuple(sorted(prem))), [p for p in prem if p in used])
    if conflict:
        first_rule_prem[("bottom",)] = [p for p in conflict[1] if p in used]
    for group in first_rule_prem.values():
        if len(group) > per_row:
            rows = [group[i:i + per_row] for i in range(0, len(group), per_row)]
            for row in rows:
                d.same_rank([_nid(S) for S in row])
            for r1, r2 in zip(rows, rows[1:]):
                d.edge(_nid(r1[0]), _nid(r2[0]), ls="invis")
    rule_ids: Dict[Tuple, str] = {}
    for S, sign, why, prem in steps:
        lab = L.rule(why, eco.n)
        key = (lab.key(), tuple(sorted(prem)))
        if key not in rule_ids:
            rid = rule_ids[key] = f"r{len(rule_ids)}"
            _rule_node(d, rid, lab)
            for p in prem:
                d.edge(_nid(p), rid)
        _sign_node(d, eco, S, sign, "forced", label_style)
        d.edge(rule_ids[key], _nid(S))
    if conflict:
        _bottom_node(d, eco, conflict[0])
        for p in conflict[1]:
            d.edge(_nid(p), "bottom")
    if legend:
        _legend(d, branch=bool(branches))
    return d


def certificate_diagram(eco: Ecology, result: SearchResult, backward_slice: bool = True, title=None,
                        label_style: str = "plucker", legend: bool = True) -> Diagram:
    """The certificate DAG of a kill at the root (see the module docstring for the colours)."""
    if result.root_conflict is None:
        raise ValueError("no root conflict: nothing to certify")
    steps, conflict_premises, initial = certificate_steps(eco, result, backward_slice)
    if title is None:
        title = _default_title(eco, "Certificate")
    return _flat_diagram(eco, steps, initial, conflict=(result.root_conflict, conflict_premises),
                         title=title, label_style=label_style, legend=legend)


def certificate_dot(eco: Ecology, result: SearchResult, **kw) -> str:
    """Graphviz source (plain-text labels) of the certificate DAG."""
    kw.pop("plucker", None)
    return certificate_diagram(eco, result, **kw).to_dot()


def _write(d: Diagram, path: str) -> None:
    """``path`` ending in .dot writes the DOT only; otherwise it is a stem (or .png/.pdf) → png + pdf + dot."""
    if path.endswith(".dot"):
        with open(path, "w") as f:
            f.write(d.to_dot())
        return
    stem = path.rsplit(".", 1)[0] if path.endswith((".png", ".pdf")) else path
    d.save(stem)


def write_certificate(eco: Ecology, result: SearchResult, path: str, **kw) -> None:
    kw.pop("plucker", None)
    _write(certificate_diagram(eco, result, **kw), path)


# --------------------------------------------------------------------------- completion derivation

def _structured_diagram(eco: Ecology, d_, index: int, title, label_style: str, initial_band: bool,
                        edges: str, per_row: int = 10, legend: bool = True) -> Diagram:
    """
    One band per search stage: the signs read off the entries (optional), the root propagation,
    then each branch choice with the signs it forces.  A forced sign is one node carrying its rule
    in small print; edges join derived signs (branch choices and forced signs) to what they force.
    Consecutive branches that force nothing share one band.  ``edges="band"`` draws only the
    edges inside a band and names the other premises in small print instead.
    """
    initial = eco.idx.to_dict(eco.initial_table())
    steps_all, _ = _steps(eco, d_.log, d_.gp_log)
    derived = {S for S, *_ in steps_all} | {S for S, _ in d_.branches}
    n_init = sum(1 for v in initial.values() if v != UNKNOWN)
    if title is None:
        title = _default_title(eco, f"Completion {index}") + [
            (f"{n_init} signs from the entries, {len(d_.branches)} branch choices, {len(steps_all)} forced", 9, style.GREY)]
    d = Diagram(title=title, nodesep=0.16, ranksep=0.3)
    fs_small = d.fontsize * 0.68

    rows = []
    if initial_band:
        init_nodes = sorted(S for S, v in initial.items() if v != UNKNOWN)
        d.cluster("initial", label=f"read off the entries ({n_init} signs)", color="#c4cfe6")
        rows = [init_nodes[i:i + per_row] for i in range(0, len(init_nodes), per_row)]
        for row in rows:
            for S in row:
                _sign_node(d, eco, S, initial[S], "initial", label_style, cluster="initial")
            d.same_rank([_nid(S) for S in row])
        for r1, r2 in zip(rows, rows[1:]):
            d.edge(_nid(r1[0]), _nid(r2[0]), ls="invis")

    stages = list(d_.stages())
    bands: List[List[Tuple]] = []
    for k, (branch, log, gp_log) in enumerate(stages):
        steps, _ = _steps(eco, log, gp_log)
        if branch is not None and not steps and bands and all(not st for _, _, st in bands[-1]) \
                and bands[-1][0][1] is not None:
            bands[-1].append((k, branch, steps))
        else:
            bands.append([(k, branch, steps)])
    band_of: Dict = {}
    for b, band in enumerate(bands):
        for _, branch, steps in band:
            if branch:
                band_of[branch[0]] = b
            for S, *_rest in steps:
                band_of[S] = b

    bottom = [_nid(rows[-1][0])] if rows else []
    for b, band in enumerate(bands):
        in_band = {S for _, _, steps in band for S, *_r in steps} | {br[0] for _, br, _ in band if br}
        depth: Dict = {}
        for _, branch, steps in band:
            if branch:
                depth[branch[0]] = 0
            for S, sign, why, prem in steps:
                inner = [depth[p] for p in prem if p in in_band and p in depth]
                depth[S] = 1 + max(inner) if inner else 1
        has_branch = band[0][1] is not None
        top = [S for S in depth if depth[S] == (0 if has_branch else 1)]
        deepest = max(depth.values()) if depth else 0
        n_forced = sum(len(steps) for _, _, steps in band)
        if not has_branch:
            head = f"root propagation ({n_forced} forced)"
        elif len(band) == 1:
            k, (S, sgn), _ = band[0]
            head = f"branch {k}: {L.plucker(S, sgn)} ({n_forced} forced)"
        else:
            head = (f"branches {band[0][0]}–{band[-1][0]}: "
                    + ", ".join(L.plucker(S, sgn) for _, (S, sgn), _ in band) + " (nothing forced)")
        cid = f"band{b}"
        d.cluster(cid, label=head, color=style.LIGHT)
        branch_ids = []
        for k, branch, steps in band:
            if branch is not None:
                S, sgn = branch
                _sign_node(d, eco, S, sgn, "branch", label_style, cluster=cid, note=f"branch {k}")
                branch_ids.append(_nid(S))
            for S, sign, why, prem in steps:
                lab = L.rule(why, eco.n)
                small = lab.headline()
                named = [p for p in prem if p not in derived]
                if edges == "band":
                    named += [p for p in prem if p in derived and band_of.get(p) != b]
                if named:
                    small += " · with " + ", ".join(
                        L.plucker(p) + (f" [band {band_of[p]}]" if p in derived else "") for p in named)
                _sign_node(d, eco, S, sign, "forced", label_style, cluster=cid, note=small)
                for p in prem:
                    if p in derived and (edges != "band" or band_of.get(p) == b):
                        d.edge(_nid(p), _nid(S))
                if branch is not None and depth[S] == 1 and not any(p in in_band for p in prem):
                    d.edge(_nid(branch[0]), _nid(S), ls="invis")
        if len(branch_ids) > 1:
            d.same_rank(branch_ids)
        for x in bottom:
            for t in top:
                d.edge(x, _nid(t), ls="invis")
        bottom = [_nid(S) for S in depth if depth[S] == deepest] or bottom
    if legend:
        _legend(d, branch=True)
    return d


def completion_diagram(eco: Ecology, result: SearchResult, index: int = 0, title=None,
                       targets: Optional[Iterable] = None, label_style: str = "plucker",
                       structured: bool = True, initial_band: Optional[bool] = None,
                       edges: Optional[str] = None, legend: bool = True) -> Diagram:
    """
    The derivation of completion ``index``.

    ``structured`` (default): one band per search stage.  ``initial_band=True`` also draws the
    table of entry-level signs as the first band.  ``edges="all"`` draws every dependency between
    derived signs, ``"band"`` (default) only those inside a band and names the rest in small print.  ``structured=False`` draws the flat DAG with explicit rule boxes, as the
    certificate does; ``targets`` (bases) then restricts it to the steps those signs depend on.
    """
    if not result.derivations:
        raise ValueError("no derivations: the search found no completion, or ran with keep_derivations=False")
    d_ = result.derivations[index]
    if structured and targets is None:
        if initial_band is None:
            initial_band = False
        if edges is None:
            edges = "band"
        return _structured_diagram(eco, d_, index, title, label_style, initial_band, edges, legend=legend)
    steps, _ = _steps(eco, d_.log, d_.gp_log)
    if targets is not None:
        steps = _backward_slice(steps, targets)
    if title is None:
        title = _default_title(eco, f"Completion {index}") + [
            (f"{len(d_.branches)} branch choices, {len(steps)} forced signs", 9, style.GREY)]
    return _flat_diagram(eco, steps, eco.idx.to_dict(eco.initial_table()), branches=d_.branches, title=title,
                         label_style=label_style, legend=legend)


def completion_dot(eco: Ecology, result: SearchResult, index: int = 0, **kw) -> str:
    kw.pop("plucker", None)
    kw.pop("show_unused_initial", None)
    return completion_diagram(eco, result, index=index, **kw).to_dot()


def write_completion(eco: Ecology, result: SearchResult, path: str, index: int = 0, **kw) -> None:
    kw.pop("plucker", None)
    kw.pop("show_unused_initial", None)
    _write(completion_diagram(eco, result, index=index, **kw), path)


# --------------------------------------------------------------------------- search tree

def search_tree_text(eco: Ecology, result: SearchResult, plucker: bool = True, max_depth: Optional[int] = None) -> str:
    """
    The binary search tree of ``complete(..., record_tree=True)`` as an indented text tree.
    Each line: the sign chosen, how many further signs the propagation forced, how many stay
    open, and either the coordinate branched on next, the completion number, or the conflict.
    """
    if not result.tree:
        raise ValueError("no tree: run complete(..., record_tree=True)")
    children: Dict[Optional[int], List[Dict]] = {}
    for nd in result.tree:
        children.setdefault(nd["parent"], []).append(nd)
    name = (lambda S: L.plucker(S, fmt="text")) if plucker else (lambda S: L.minor(eco, S, fmt="text"))

    def describe(nd):
        if nd["status"] == "conflict":
            return f"⊥ {L.rule(nd['reason'], eco.n).headline('text')}"
        if nd["status"] == "stopped":
            return f"{nd['unknown']} open … not searched (limit reached)"
        tail = (f"→ branch on {name(nd['branch_on'])}" if nd.get("branch_on") is not None
                else f"✓ completion {nd.get('completion')}")
        forced = f"{nd['forced']} forced, " if nd["forced"] else ""
        return f"{forced}{nd['unknown']} open {tail}"

    lines = []

    def walk(nd, prefix, last, is_root):
        if is_root:
            lines.append(f"root: {describe(nd)}")
            child_prefix = ""
        else:
            sign = "+" if nd["sign"] > 0 else "−"
            lines.append(f"{prefix}{'└─' if last else '├─'} {sign}  {describe(nd)}")
            child_prefix = prefix + ("   " if last else "│  ")
        if max_depth is not None and nd["depth"] >= max_depth:
            kids = children.get(nd["id"], [])
            if kids:
                lines.append(f"{child_prefix}└─ … ({_subtree_size(children, nd['id'])} nodes below)")
            return
        for i, k in enumerate(children.get(nd["id"], [])):
            walk(k, child_prefix, i == len(children[nd["id"]]) - 1, False)

    walk(result.tree[0], "", True, True)
    return "\n".join(lines)


def _subtree_size(children, node_id) -> int:
    return sum(1 + _subtree_size(children, k["id"]) for k in children.get(node_id, []))


def tree_children(result: SearchResult) -> Dict[Optional[int], List[Dict]]:
    children: Dict[Optional[int], List[Dict]] = {}
    for nd in result.tree:
        children.setdefault(nd["parent"], []).append(nd)
    return children


def repeated_subtrees(result: SearchResult, max_classes: int = 2, min_size: int = 3):
    """
    Isomorphism classes of subtrees that the search visits more than once — same coordinate
    branched on, same signs chosen, same propagation counts and same leaves, recursively.
    Returns [(root_ids, copies, size, root_basis)] for the classes chosen greedily by copies × size,
    skipping a class nested inside one already chosen; at most ``max_classes`` classes.
    """
    children = tree_children(result)
    sig, size = {}, {}

    def walk(nd):
        if nd["status"] == "completion":
            s = ("C",)
        elif nd["status"] == "conflict":
            s = ("X", nd["reason"].split(":")[0])
        elif nd["status"] == "stopped":
            s = ("S",)
        else:
            s = (nd["branch_on"], tuple((k["sign"], k.get("forced"), walk(k)) for k in children.get(nd["id"], [])))
        sig[nd["id"]] = s
        size[nd["id"]] = 1 + sum(size[k["id"]] for k in children.get(nd["id"], []))
        return s

    walk(result.tree[0])
    roots: Dict = {}
    for i, s in sig.items():
        if len(s) == 2 and isinstance(s[1], tuple) and size[i] >= min_size:
            roots.setdefault(s, []).append(i)
    classes = sorted(((ids, len(ids), size[ids[0]], s[0]) for s, ids in roots.items() if len(ids) >= 2),
                     key=lambda t: -(t[1] * t[2]))
    by_id = {nd["id"]: nd for nd in result.tree}

    def inside(i, chosen_roots):
        while i is not None:
            if i in chosen_roots:
                return True
            i = by_id[i]["parent"]
        return False

    chosen, chosen_roots = [], set()
    for ids, copies, sz, root in classes:
        if len(chosen) >= max_classes:
            break
        if any(inside(by_id[i]["parent"], chosen_roots) for i in ids):
            continue
        chosen.append((ids, copies, sz, root))
        chosen_roots.update(ids)
    return chosen


def subtree_ids(result: SearchResult, root_id: int) -> List[int]:
    children = tree_children(result)
    out, stack = [], [root_id]
    while stack:
        i = stack.pop()
        out.append(i)
        stack.extend(k["id"] for k in children.get(i, []))
    return out


def tree_title(eco: Ecology, result: SearchResult, use: Optional[Sequence[str]] = None) -> List:
    """Default title of a search-tree drawing: the pattern, the rules and the outcome."""
    head = f"Search tree: {eco.name}" if getattr(eco, "name", "") else "Search tree"
    parts = [L.pattern(eco)]
    if use:
        parts.append("rules " + ", ".join(use))
    if result.root_unknowns is not None:
        parts.append(f"{result.root_unknowns} coordinates open at the root")
    outcome = (f"{len(result.completions)} completions, {result.conflicts} contradictions"
               + (", search stopped early" if result.stopped else ""))
    return [head, (" · ".join(parts), 9.5, style.GREY), (outcome, 9.5, style.GREY)]


def _open_shade(nd, root_open: int) -> str:
    """Grey level of a branch node by the fraction of coordinates still open (white = as open as the root)."""
    if not root_open or nd.get("unknown") is None:
        return "white"
    f = max(0.0, min(1.0, nd["unknown"] / root_open))
    g = int(255 - (1 - f) * 110)
    return f"#{g:02x}{g:02x}{g:02x}"


def _edge_width(forced, base: float = 1.2, per_sign: float = 0.14, cap: float = 3.4) -> float:
    return min(cap, base + per_sign * (forced or 0))


def search_tree_diagram(eco: Ecology, result: SearchResult, title=None, label_style: str = "plucker",
                        show_forced: bool = False, use: Optional[Sequence[str]] = None,
                        completion_numbers: bool = True, legend: bool = True,
                        width_by_forced: bool = True, repeats: int = 2, shade_open: bool = False) -> Diagram:
    """
    The binary search tree laid out by Graphviz.

    Colour channels: the edge to a child is red when + was chosen and blue when −, and its width
    grows with the number of further signs the propagation forced (``width_by_forced``; with
    ``show_forced`` the number is printed too); a ⊥ leaf is filled with the colour of the rule that
    closed it; completions are green, unexplored subtrees grey.  Subtrees that the search visits
    several times — same coordinates, same choices, same outcome — share a pale hull
    (``repeats`` classes at most, see ``repeated_subtrees``).  ``shade_open`` greys a branch node
    by the number of coordinates still open there (white = as open as the root).
    """
    if not result.tree:
        raise ValueError("no tree: run complete(..., record_tree=True)")
    d = Diagram(title=title if title is not None else tree_title(eco, result, use), nodesep=0.1, ranksep=0.35)
    name = (lambda S: L.plucker(S)) if label_style != "minor" else (lambda S: L.minor(eco, S))
    fs = d.fontsize
    root_open = result.tree[0].get("unknown") or 0
    rules_seen: Dict[str, str] = {}
    for nd in result.tree:
        i = f"n{nd['id']}"
        st = nd["status"]
        if st == "conflict":
            lab = L.rule(nd["reason"], eco.n)
            rules_seen[lab.rule] = lab.headline()
            d.node(i, [(r"$\bot$", fs * 1.2, "white"), (lab.headline(), fs * 0.7, "white")],
                   fill=style.rule_color(lab.rule), edge=style.rule_color(lab.rule), pad=(5.0, 2.0), gap=1.0)
        elif st == "stopped":
            d.node(i, [("…", fs, style.GREY)], fill="white", edge=style.STOPPED, ls="--", shape="circle", pad=(3.0, 3.0))
        elif st == "completion":
            text = str(nd.get("completion")) if completion_numbers else "✓"
            d.node(i, [(text, fs * 0.85, "white")], fill=style.COMPLETION, edge=style.COMPLETION, shape="circle", pad=(3.0, 3.0))
        else:
            d.node(i, [name(nd["branch_on"])], fill=_open_shade(nd, root_open) if shade_open else "white",
                   edge=style.INK if st == "root" else style.GREY, lw=1.6 if st == "root" else 0.9,
                   shape="circle", pad=(4.0, 4.0))
        if nd["parent"] is not None:
            lab = None
            if show_forced and nd.get("forced"):
                lab = f"+{nd['forced']}"
            d.edge(f"n{nd['parent']}", i, color=style.sign_line(nd["sign"]),
                   lw=_edge_width(nd.get("forced")) if width_by_forced else 1.4, label=lab,
                   size=fs * 0.7, label_color=style.GREY)
    classes = repeated_subtrees(result, max_classes=repeats) if repeats else []
    # one cluster per copy of a repeated subtree (dot keeps clusters contiguous); outer classes are
    # created first so that an inner copy can name the copy it sits in as its parent, and a node
    # belongs to the innermost copy containing it
    by_id = {nd["id"]: nd for nd in result.tree}
    copy_root: Dict[int, str] = {}
    tint = {root: style.HULLS[k % len(style.HULLS)] for k, (_, _, _, root) in enumerate(classes)}
    for ids, copies, sz, root in sorted(classes, key=lambda t: -t[2]):
        for r in ids:
            parent, p = None, by_id[r]["parent"]
            while p is not None and parent is None:
                parent, p = copy_root.get(p), by_id[p]["parent"]
            cid = copy_root[r] = f"copy{r}"
            d.cluster(cid, label=f"×{copies}", fill=tint[root], color="none",
                      parent=parent, label_pos="tr", size=fs * 0.8, margin=6.0)
            for j in subtree_ids(result, r):
                nid = f"n{j}"
                old = d.nodes[nid]["cluster"]
                if old is not None:
                    d.clusters[old]["nodes"].remove(nid)
                d.nodes[nid]["cluster"] = cid
                d.clusters[cid]["nodes"].append(nid)
    if legend:
        d.legend_item("chose +", fill=style.POS_LINE, edge=style.POS_LINE)
        d.legend_item("chose −", fill=style.NEG_LINE, edge=style.NEG_LINE)
        if width_by_forced:
            d.legend_item("thicker: more signs forced", fill="white", edge=style.GREY)
        d.legend_item("completion", fill=style.COMPLETION, edge=style.COMPLETION, shape="circle")
        for rule, head in sorted(rules_seen.items()):
            d.legend_item(f"⊥ by {head}", fill=style.rule_color(rule), edge=style.rule_color(rule))
        if any(nd["status"] == "stopped" for nd in result.tree):
            d.legend_item("not searched", fill="white", edge=style.STOPPED, ls="--", shape="circle")
        for k, (ids, copies, sz, root) in enumerate(classes):
            d.legend_item(f"same subtree searched {copies}×  (from {L.plucker(root)}, {sz} nodes)",
                          fill=style.HULLS[k % len(style.HULLS)], edge="none")
        if shade_open:
            d.legend_item("darker: fewer coordinates open", fill="#a0a0a0", edge=style.GREY, shape="circle")
    return d


def search_tree_dot(eco: Ecology, result: SearchResult, **kw) -> str:
    kw.pop("plucker", None)
    return search_tree_diagram(eco, result, **kw).to_dot()


def write_search_tree(eco: Ecology, result: SearchResult, path: str, **kw) -> None:
    kw.pop("plucker", None)
    _write(search_tree_diagram(eco, result, **kw), path)
