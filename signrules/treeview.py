"""
Interactive viewer of a sign-completion search tree: one self-contained HTML file.

    result = complete(eco, use=("R1", "R2n"), record_tree=True)
    save_tree_viewer(eco, result, "tree.html", use=("R1", "R2n"))

    python -m signrules.treeview --a=1,-1,-1 --B=1,1,1/-1,1,1/-1,-1,1 --rules R1,R2n --out tree.html
    python -m signrules.treeview                     # prompts for the pattern (needs the ``live`` extra)

The console script ``signrules-tree`` is the same program; ``signrules-live --html FILE`` writes the
same page after the dashboard.

The page shows the binary search tree (edges red when the branch chose +, blue when it chose −;
completions green, contradictions black squares framed in the colour of the rule that closed them,
unexplored subtrees grey) and, for the selected node,

* the path of branch choices from the root (clickable);
* what the signs known at that node say about **stability** (each coefficient c_k = Σ ξ_I det B_II
  of the characteristic polynomial: holds / possible / violated / needs an open minor, with one
  chip per principal minor) and **feasibility** (Cramer's rule for every x*_j against the sign
  of det B), and the whole table of minors; signs decided at the node are ringed;
* every Plücker sign forced at that node, one row each: rule tag (GP, R1–R8, in the rule colours
  of ``style.RULE_COLORS``), coordinate, minor of B; the premises are in the tooltip;
* for a contradiction, the rule, the clash and its premises.

Interaction: click a node to inspect it; double-click (or the button) to fold / unfold its subtree;
space unfolds the whole subtree below the selected node (folds it if nothing below is folded);
wheel to zoom, drag to pan, ``f`` to fit; arrow keys walk the tree (↑ parent, ← + child,
→ − child); the search box highlights the nodes where a coordinate (e.g. ``0257``) is forced or
branched on.

No external scripts: the data is embedded as JSON and drawn as SVG by plain JavaScript, so the
file works offline and can be shared as it is.
"""

import argparse
import html
import itertools
import json
from typing import Dict, List, Optional, Sequence

from . import labels as L
from . import style
from .certificate import _steps
from .plucker_index import UNKNOWN
from .rules import Ecology
from .search import SearchResult, complete


def _coord(S) -> str:
    return "".join(str(i) for i in S)


def _step_record(eco: Ecology, S, sign, why: str, premises) -> Dict:
    lab = L.rule(why, eco.n)
    return dict(S=_coord(S), sign=int(sign), coord=L.plucker(S, sign, "html"),
                minor=L.minor(eco, S, sign, "html"), rule=lab.rule, head=lab.headline("html"),
                detail=L.convert(lab.detail, "html").replace("\n", " ") if lab.detail else "",
                clause=L.convert(lab.clause, "html") if lab.clause else "",
                premises=[L.plucker(P, None, "html") for P in premises],
                premS=[_coord(P) for P in premises])


def tree_data(eco: Ecology, result: SearchResult, use: Optional[Sequence[str]] = None) -> Dict:
    """Everything the viewer needs, as plain JSON-able data."""
    if not result.tree:
        raise ValueError("run complete(..., record_tree=True) to record the tree")
    bases = eco.idx.bases
    key_to_basis = {k: bases[i] for k, (i, _) in eco._minor_key.items()}

    def to_basis(k):
        return k if (isinstance(k, tuple) and k and isinstance(k[0], int)) else key_to_basis[tuple(map(tuple, k))]

    nodes = []
    for nd in result.tree:
        steps, _ = _steps(eco, nd.get("log", []), nd.get("gp_log", []))
        rec = dict(id=nd["id"], parent=nd["parent"], depth=nd["depth"], status=nd["status"],
                   forced=nd.get("forced"), unknown=nd.get("unknown"), completion=nd.get("completion"),
                   steps=[_step_record(eco, *st) for st in steps])
        if nd.get("basis") is not None:
            rec["choiceS"] = _coord(nd["basis"])
            rec["choice"] = L.plucker(nd["basis"], nd["sign"], "html")
            rec["choiceMinor"] = L.minor(eco, nd["basis"], nd["sign"], "html")
            rec["sign"] = int(nd["sign"])
        if nd.get("branch_on") is not None:
            rec["branchS"] = _coord(nd["branch_on"])
            rec["branch"] = L.plucker(nd["branch_on"], None, "html")
            rec["branchMinor"] = L.minor(eco, nd["branch_on"], None, "html")
        if nd["status"] == "conflict":
            reason = nd.get("reason_full") or nd.get("reason", "")
            lab = L.rule(reason, eco.n)
            rec["conflict"] = dict(rule=lab.rule, lines=L.conflict(eco, reason, "html")[1:],
                                   premises=[L.plucker(to_basis(p), None, "html") for p in nd.get("premises", [])])
        nodes.append(rec)

    title = L.pattern(eco, "html")
    n = eco.n
    coords, init = {}, {}
    table0 = eco.initial_table()
    for (P, Q), (i, c) in eco._minor_key.items():
        S = bases[i]
        code = _coord(S)
        coords[code] = dict(m=L.minor(eco, S, None, "html"), c=int(c), k=len(P), pr=P == Q, P=list(P), Q=list(Q))
        init[code] = int(table0[i])                 # 0 = UNKNOWN
    stab = []
    for k in range(1, n + 1):
        items = [_coord(bases[eco._minor_key[(I, I)][0]]) for I in itertools.combinations(range(n), k)]
        stab.append(dict(k=k, items=items))
    full = tuple(range(n))
    feas = []
    for j in range(n):
        terms = []
        for k in range(n):
            P = tuple(i for i in full if i != k)
            Q = tuple(i for i in full if i != j)
            i, c = eco._minor_key[(P, Q)]
            terms.append(dict(k=k, S=_coord(bases[i]), f=int((-1) ** (k + j) * eco.a[k] * c)))
        feas.append(dict(j=j, terms=terms))
    detI, detc = eco._minor_key[(full, full)]
    detS = bases[detI]
    signs = dict(coords=coords, init=init, stab=stab, feas=feas, det=dict(S=_coord(detS), c=int(detc)))
    return dict(signs=signs,
        name=eco.name or "pattern", n=eco.n, a=[int(x) for x in eco.a],
        B=[[int(x) for x in row] for row in eco.B],
        title=title, rules=list(use) if use else None,
        summary=result.summary(), seconds=round(result.seconds, 3),
        stopped=bool(result.stopped), nodes=nodes,
        colors=dict(pos=style.POS, neg=style.NEG, posLine=style.POS_LINE, negLine=style.NEG_LINE,
                    posTint=style.POS_TINT, negTint=style.NEG_TINT, completion=style.COMPLETION,
                    bottom=style.BOTTOM, stopped=style.STOPPED, rules=style.RULE_COLORS))


def tree_viewer_html(eco: Ecology, result: SearchResult, use: Optional[Sequence[str]] = None) -> str:
    """The viewer page as a string (``result`` from ``complete(..., record_tree=True)``; ``use``: the rules it ran)."""
    data = tree_data(eco, result, use)
    payload = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    return _TEMPLATE.replace("__TITLE__", html.escape(f"Sign tree · {data['name']}")).replace("__DATA__", payload)


def save_tree_viewer(eco: Ecology, result: SearchResult, path: str, use: Optional[Sequence[str]] = None) -> None:
    """Write the self-contained HTML viewer of a recorded search tree to ``path``."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(tree_viewer_html(eco, result, use))


# --------------------------------------------------------------------------- command line

def _parse(a_str: str, B_str: str):
    a = [int(x) for x in a_str.split(",")]
    B = [[int(x) for x in row.split(",")] for row in B_str.split("/")]
    return a, B


def main(argv=None):
    """Command line of the viewer (``signrules-tree``)."""
    ap = argparse.ArgumentParser(description="Write an interactive HTML viewer of a sign-completion search tree.")
    ap.add_argument("--a", help="growth-rate signs, e.g. 1,-1,-1")
    ap.add_argument("--B", help="signs of B, rows separated by /, e.g. 1,1,1/-1,1,1/-1,-1,1")
    ap.add_argument("--rules", default="R1,R2n", help="rules to use (default R1,R2n; 'all' for every rule)")
    ap.add_argument("--max-nodes", type=int, default=20000, help="stop the search after this many nodes")
    ap.add_argument("--name", default="")
    ap.add_argument("--out", default="tree.html")
    args = ap.parse_args(argv)
    if not (args.a and args.B):
        from rich.prompt import Prompt
        args.a = Prompt.ask("growth-rate signs a", default="1,-1,-1")
        args.B = Prompt.ask("signs of B, rows separated by /", default="1,1,1/-1,1,1/-1,-1,1")
        args.rules = Prompt.ask("rules ('all' or a list such as R1,R2n)", default=args.rules)
        args.out = Prompt.ask("output file", default=args.out)
    a, B = _parse(args.a, args.B)
    use = Ecology.ALL_RULES if args.rules == "all" else tuple(r for r in args.rules.split(",") if r)
    eco = Ecology(a, B, args.name)
    res = complete(eco, use=use, record_tree=True, max_nodes=args.max_nodes)
    save_tree_viewer(eco, res, args.out, use=use)
    print(f"{res.summary()}  ->  {args.out}")


_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root {
  --bg: #ffffff; --panel: #f7f7f5; --ink: #222; --muted: #6b6b6b; --line: #dcdcdc; --node: #ffffff;
  --sel: #f2b33d; --hit: #8a5cf6; --bottom: #222222;
  --ok: #2e7d32; --bad: #c62828; --open: #b26a00;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #17181a; --panel: #1f2023; --ink: #e8e8e6; --muted: #9a9a96; --line: #34353a; --node: #26272b;
    --sel: #f2b33d; --hit: #a98bff; --bottom: #e8e8e6;
    --ok: #6fcf73; --bad: #ff7b72; --open: #e3a33b;
  }
}
:root[data-theme="dark"] {
  --bg: #17181a; --panel: #1f2023; --ink: #e8e8e6; --muted: #9a9a96; --line: #34353a; --node: #26272b;
  --bottom: #e8e8e6; --ok: #6fcf73; --bad: #ff7b72; --open: #e3a33b;
}
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; background: var(--bg); color: var(--ink);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
#app { display: grid; grid-template-columns: 1fr 380px; grid-template-rows: auto 1fr; height: 100vh; }
header { grid-column: 1 / 3; display: flex; flex-wrap: wrap; align-items: center; gap: 16px;
  padding: 10px 16px; border-bottom: 1px solid var(--line); }
header h1 { font-size: 16px; margin: 0; font-weight: 600; }
header .meta { color: var(--muted); }
.mat { display: inline-grid; gap: 2px; vertical-align: middle; }
.mat span { width: 16px; height: 16px; border-radius: 3px; display: grid; place-items: center;
  font-size: 11px; color: #fff; font-weight: 600; }
.tools { margin-left: auto; display: flex; gap: 8px; align-items: center; }
input[type=search] { background: var(--panel); color: var(--ink); border: 1px solid var(--line);
  border-radius: 6px; padding: 5px 8px; width: 150px; }
button { background: var(--panel); color: var(--ink); border: 1px solid var(--line); border-radius: 6px;
  padding: 5px 10px; cursor: pointer; }
button:hover { border-color: var(--muted); }
#view { position: relative; overflow: hidden; cursor: grab; }
#view.drag { cursor: grabbing; }
svg { width: 100%; height: 100%; display: block; user-select: none; }
#panel { border-left: 1px solid var(--line); background: var(--panel); overflow-y: auto; padding: 14px 16px; }
#panel h2 { font-size: 15px; margin: 0 0 6px; }
#panel h3 { font-size: 12px; text-transform: uppercase; letter-spacing: .04em; color: var(--muted);
  margin: 16px 0 6px; }
.crumbs { display: flex; flex-wrap: wrap; gap: 4px; }
.crumb { border: 1px solid var(--line); border-radius: 12px; padding: 1px 8px; cursor: pointer; font-size: 13px; }
.crumb.pos { border-color: var(--posl); } .crumb.neg { border-color: var(--negl); }
.rulegrp { margin: 8px 0; border-left: 4px solid var(--rc); padding-left: 8px; }
.rulegrp .rh { font-weight: 600; color: var(--rc); }
.step { margin: 3px 0; font-size: 13px; }
.step .pr { color: var(--muted); font-size: 12px; }
.sg-pos { color: var(--posi); font-weight: 600; } .sg-neg { color: var(--negi); font-weight: 600; }
.box { border: 1px solid var(--line); border-radius: 8px; padding: 8px 10px; background: var(--bg); }
#panel h3.math { text-transform: none; letter-spacing: 0; font-size: 13px; }
.steps { display: grid; grid-template-columns: auto auto 1fr; gap: 1px 8px; font-size: 13px; align-items: baseline; }
.steps .tag { font-size: 11px; font-weight: 700; color: #fff; border-radius: 3px; padding: 0 4px; }
.steps .m { color: var(--muted); }
.steps > span { cursor: default; }
.verdict { font-size: 12px; font-weight: 600; border-radius: 4px; padding: 0 6px; white-space: nowrap; }
.v-ok { background: color-mix(in srgb, var(--ok) 22%, transparent); color: var(--ok); }
.v-bad { background: color-mix(in srgb, var(--bad) 22%, transparent); color: var(--bad); }
.v-open { background: color-mix(in srgb, var(--open) 22%, transparent); color: var(--open); }
.v-need { background: color-mix(in srgb, var(--open) 12%, transparent); color: var(--muted); }
.srow { display: grid; grid-template-columns: 52px 1fr auto; gap: 8px; align-items: center; margin: 3px 0; font-size: 13px; }
.chips { display: flex; flex-wrap: wrap; gap: 2px; }
.chip { min-width: 14px; height: 14px; border-radius: 3px; font-size: 10px; line-height: 14px; text-align: center;
  color: #fff; font-weight: 700; padding: 0 2px; }
.chip.q { background: transparent; color: var(--muted); border: 1px dashed var(--muted); line-height: 12px; }
.chip.new { outline: 2px solid var(--sel); outline-offset: 1px; }
details summary { cursor: pointer; color: var(--muted); font-size: 12px; text-transform: uppercase; letter-spacing: .04em; margin: 16px 0 6px; }
.legend { display: flex; flex-wrap: wrap; gap: 10px; font-size: 12px; color: var(--muted); }
.legend i { display: inline-block; width: 18px; height: 3px; vertical-align: middle; margin-right: 4px; }
.legend b { display: inline-block; width: 10px; height: 10px; border-radius: 50%; vertical-align: middle; margin-right: 4px; }
.muted { color: var(--muted); }
.lbl { font-size: 11px; fill: var(--ink); pointer-events: none; }
.badge { font-size: 10px; fill: var(--muted); pointer-events: none; }
@media (max-width: 820px) {
  #app { grid-template-columns: 1fr; grid-template-rows: auto 60vh auto; height: auto; }
  header { grid-column: 1; } #panel { border-left: 0; border-top: 1px solid var(--line); }
  #view { height: 60vh; }
}
</style>
</head>
<body>
<div id="app">
  <header>
    <h1 id="ttl"></h1>
    <span id="pat"></span>
    <span class="meta" id="sum"></span>
    <div class="tools">
      <input type="search" id="q" placeholder="find p… e.g. 0257" aria-label="find a Plücker coordinate">
      <button id="fold" title="fold or unfold the selected subtree (double-click)">fold</button>
      <button id="all" title="unfold every subtree, or fold back to a readable depth">unfold all</button>
      <button id="fit" title="fit to window (f)">fit</button>
      <button id="theme" title="light / dark">◐</button>
    </div>
  </header>
  <div id="view"><svg id="svg"><g id="world"></g></svg></div>
  <aside id="panel"></aside>
</div>
<script>
const D = __DATA__;
const C = D.colors;
const root = document.documentElement;
root.style.setProperty('--posl', C.posLine); root.style.setProperty('--negl', C.negLine);
root.style.setProperty('--posi', C.pos); root.style.setProperty('--negi', C.neg);

// ---------------------------------------------------------------- data
const byId = new Map(D.nodes.map(n => [n.id, n]));
for (const n of D.nodes) n.kids = [];
for (const n of D.nodes) if (n.parent !== null && n.parent !== undefined) byId.get(n.parent).kids.push(n);
for (const n of D.nodes) n.kids.sort((x, y) => (y.sign || 0) - (x.sign || 0));   // + child left
const ROOT = D.nodes.find(n => n.parent === null || n.parent === undefined);
const folded = new Set();
let selected = ROOT, query = '';

// ---------------------------------------------------------------- header
document.getElementById('ttl').textContent = D.name;
const pat = document.getElementById('pat');
const sgn = v => v > 0 ? '+' : '−';
let ph = '<span class="muted">sign a</span> ' + D.a.map(v => `<span class="${v > 0 ? 'sg-pos' : 'sg-neg'}">${sgn(v)}</span>`).join(' ');
if (D.B) {
  ph += ' &nbsp;<span class="muted">sign B</span> <span class="mat" style="grid-template-columns:repeat(' + D.n + ',16px)">' +
    D.B.flat().map(v => `<span style="background:${v > 0 ? C.pos : C.neg}">${sgn(v)}</span>`).join('') + '</span>';
}
pat.innerHTML = ph;
document.getElementById('sum').textContent = D.summary + (D.rules ? ' · rules ' + D.rules.join(', ') : '') + ` · ${D.seconds}s`;

// ---------------------------------------------------------------- layout (tidy binary tree)
const DX = 34, DY = 64;
function layout() {
  let next = 0;
  (function place(n) {
    const kids = folded.has(n.id) ? [] : n.kids;
    if (!kids.length) { n.x = next++ * DX; }
    else { kids.forEach(place); n.x = (kids[0].x + kids[kids.length - 1].x) / 2; }
    n.y = n.depth * DY;
    n.visibleKids = kids;
  })(ROOT);
}
function size(n) { return 1 + n.kids.reduce((s, k) => s + size(k), 0); }

// ---------------------------------------------------------------- drawing
const svg = document.getElementById('svg'), world = document.getElementById('world');
const NS = 'http://www.w3.org/2000/svg';
const el = (tag, attrs, parent) => { const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; };
const ruleColor = r => C.rules[r] || '#888';
const hits = () => {
  if (!query) return new Set();
  const s = new Set();
  for (const n of D.nodes)
    if (n.branchS === query || n.choiceS === query || n.steps.some(st => st.S === query)) s.add(n.id);
  return s;
};
function draw() {
  layout();
  world.innerHTML = '';
  const H = hits();
  const gE = el('g', {}, world), gN = el('g', {}, world);
  (function walk(n) {
    for (const k of n.visibleKids) {
      const w = Math.min(1.2 + 0.14 * (k.forced || 0), 3.4);
      el('path', { d: `M${n.x},${n.y} C${n.x},${(n.y + k.y) / 2} ${k.x},${(n.y + k.y) / 2} ${k.x},${k.y}`,
        fill: 'none', stroke: k.sign > 0 ? C.posLine : C.negLine, 'stroke-width': w }, gE);
      walk(k);
    }
    const g = el('g', { transform: `translate(${n.x},${n.y})`, 'data-id': n.id, style: 'cursor:pointer' }, gN);
    if (n.id === selected.id) el('circle', { r: 15, fill: 'none', stroke: 'var(--sel)', 'stroke-width': 3 }, g);
    if (H.has(n.id)) el('circle', { r: 19, fill: 'none', stroke: 'var(--hit)', 'stroke-width': 2, 'stroke-dasharray': '3 2' }, g);
    if (n.status === 'conflict') {
      const rc = n.conflict ? ruleColor(n.conflict.rule) : 'var(--bottom)';
      el('rect', { x: -7, y: -7, width: 14, height: 14, fill: 'var(--bottom)', stroke: rc, 'stroke-width': 3 }, g);
    } else if (n.status === 'completion') {
      el('circle', { r: 8, fill: C.completion }, g);
    } else if (n.status === 'stopped') {
      el('circle', { r: 7, fill: 'var(--node)', stroke: C.stopped, 'stroke-width': 1.5, 'stroke-dasharray': '2 2' }, g);
    } else {
      el('circle', { r: 9, fill: 'var(--node)', stroke: n.status === 'root' ? 'var(--ink)' : 'var(--muted)', 'stroke-width': 1.5 }, g);
      // rule mix of the signs forced here: a small stacked bar under the node
      const cnt = {}; n.steps.forEach(s => cnt[s.rule] = (cnt[s.rule] || 0) + 1);
      const tot = n.steps.length;
      if (tot) {
        let x = -10;
        for (const r of Object.keys(cnt)) { const w = 20 * cnt[r] / tot;
          el('rect', { x, y: 12, width: w, height: 3, fill: ruleColor(r) }, g); x += w; }
      }
    }
    if (n.branchS && !folded.has(n.id)) {
      const t = el('text', { x: 12, y: -10, class: 'lbl' }, g); t.textContent = 'p' + n.branchS;
    }
    if (folded.has(n.id) && n.kids.length) {
      const t = el('text', { x: 0, y: 26, class: 'badge', 'text-anchor': 'middle' }, g);
      t.textContent = '+' + (size(n) - 1);
    }
    const tip = el('title', {}, g);
    tip.textContent = n.status + (n.choiceS ? ` · chose p${n.choiceS} ${n.sign > 0 ? '> 0' : '< 0'}` : '') +
      (n.steps.length ? ` · ${n.steps.length} forced` : '');
  })(ROOT);
}

// ---------------------------------------------------------------- side panel
const panel = document.getElementById('panel');
function pathTo(n) { const p = []; while (n) { p.unshift(n); n = n.parent === null || n.parent === undefined ? null : byId.get(n.parent); } return p; }
const statusName = { root: 'Root', branch: 'Branch node', completion: 'Completion', conflict: 'Contradiction', stopped: 'Not searched' };
function show(n) {
  let h = `<h2>${statusName[n.status] || n.status}${n.status === 'completion' ? ' ' + n.completion : ''}
    <span class="muted" style="font-weight:400">· node ${n.id}, depth ${n.depth}</span></h2>`;
  if (n.choice) h += `<div>chose <span class="${n.sign > 0 ? 'sg-pos' : 'sg-neg'}">${n.choice}</span>
    <span class="muted">(${n.choiceMinor})</span></div>`;
  const p = pathTo(n).slice(1);
  h += '<h3>Path from the root</h3><div class="crumbs">' + (p.length ? p.map(m =>
    `<span class="crumb ${m.sign > 0 ? 'pos' : 'neg'}" data-go="${m.id}">p${m.choiceS} ${m.sign > 0 ? '+' : '−'}</span>`).join('')
    : '<span class="muted">the root</span>') + '</div>';
  const facts = [];
  if (n.forced !== null && n.forced !== undefined) facts.push(`${n.forced} signs forced here`);
  if (n.unknown !== null && n.unknown !== undefined) facts.push(`${n.unknown} still open`);
  if (n.branch) facts.push(`branches on ${n.branch} <span class="muted">(${n.branchMinor})</span>`);
  if (n.kids.length) facts.push(`subtree of ${size(n) - 1} nodes`);
  if (facts.length) h += '<h3>Summary</h3><div>' + facts.join('<br>') + '</div>';
  if (n.conflict) {
    h += `<h3>Contradiction</h3><div class="box rulegrp" style="--rc:${ruleColor(n.conflict.rule)}">
      <div class="rh">⊥ ${n.conflict.lines[0] || ''}</div>${n.conflict.lines.slice(1).map(l => `<div>${l}</div>`).join('')}
      ${n.conflict.premises.length ? `<div class="pr muted">premises: ${n.conflict.premises.join(', ')}</div>` : ''}</div>`;
  }
  h += analysis(n);
  if (n.steps.length) {
    h += `<h3>Signs forced at this node · ${n.steps.length}</h3><div class="steps">`;
    for (const s of n.steps) {
      const tip = (s.head || s.rule) + (s.detail && s.detail !== s.minor ? ' · ' + strip(s.detail) : '') +
        (s.premises.length ? '\nfrom ' + s.premS.map(x => 'p' + x).join(', ') : '');
      h += `<span class="tag" style="background:${ruleColor(s.rule)}" title="${esc(tip)}">${s.rule}</span>
        <span class="${s.sign > 0 ? 'sg-pos' : 'sg-neg'}" title="${esc(tip)}">${s.coord}</span>
        <span class="m" title="${esc(tip)}">${s.minor}</span>`;
    }
    h += '</div>';
  } else if (n.status !== 'stopped') {
    h += '<h3>Signs forced at this node</h3><div class="muted">none</div>';
  }
  h += `<h3>Legend</h3><div class="legend">
    <span><i style="background:${C.posLine}"></i>chose +</span><span><i style="background:${C.negLine}"></i>chose −</span>
    <span><b style="background:${C.completion}"></b>completion</span><span><b style="background:var(--bottom);border-radius:0"></b>contradiction</span>
    ${Object.entries(C.rules).map(([r, c]) => `<span><i style="background:${c}"></i>${r}</span>`).join('')}</div>
    <div class="muted" style="margin-top:8px;font-size:12px">click: inspect · double-click: fold · space: unfold all below · wheel: zoom · drag: pan ·
    f: fit · ↑ ← → walk the tree</div>`;
  panel.innerHTML = h;
  panel.querySelectorAll('[data-go]').forEach(e => e.onclick = () => select(byId.get(+e.dataset.go)));
}

// ---------------------------------------------------------------- what the signs say at a node
const SG = D.signs;
const esc = t => String(t).replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;');
const strip = t => String(t).replace(/<[^>]+>/g, '');
function tableAt(n) {
  const t = Object.assign({}, SG.init), fresh = new Set();
  for (const m of pathTo(n)) {
    if (m.choiceS) { t[m.choiceS] = m.sign; if (m === n) fresh.add(m.choiceS); }
    for (const st of m.steps) { t[st.S] = st.sign; if (m === n) fresh.add(st.S); }
  }
  return [t, fresh];
}
const minorSign = (t, S) => t[S] ? SG.coords[S].c * t[S] : 0;
function chip(sign, label, isNew) {
  const cls = 'chip' + (sign ? '' : ' q') + (isNew ? ' new' : '');
  const bg = sign > 0 ? `background:${C.pos}` : sign < 0 ? `background:${C.neg}` : '';
  return `<span class="${cls}" style="${bg}" title="${esc(strip(label))}">${sign > 0 ? '+' : sign < 0 ? '−' : '?'}</span>`;
}
// a sum of terms with positive unknown weights must have sign `target`
function verdict(signs, target) {
  const want = signs.filter(x => x === target).length, against = signs.filter(x => x === -target).length;
  const open = signs.length - want - against;
  if (want === signs.length) return ['v-ok', 'holds'];
  if (against === signs.length) return ['v-bad', 'violated'];
  if (want) return ['v-open', 'possible'];
  return open ? ['v-need', 'needs an open term'] : ['v-bad', 'violated'];
}
const VHELP = { holds: 'every term has the required sign: true for all magnitudes',
  violated: 'every term has the wrong sign: impossible for all magnitudes',
  possible: 'terms of both signs: decided by the magnitudes',
  'needs an open term': 'no term has the required sign yet: one of the open ones must' };
function vtag(v) { return `<span class="verdict ${v[0]}" title="${VHELP[v[1]]}">${v[1]}</span>`; }
function analysis(n) {
  const [t, fresh] = tableAt(n);
  let h = '';
  // stability: c_k = sum over |I| = k of xi_I det B_II with xi_I > 0 must be positive (k = n: det B > 0)
  h += '<h3 class="math">Stability · c<sub>k</sub> = Σ<sub>|I|=k</sub> ξ<sub>I</sub> det B<sub>I</sub> &gt; 0</h3>';
  for (const row of SG.stab) {
    const sg = row.items.map(S => minorSign(t, S));
    const chips = row.items.map((S, i) => chip(sg[i], SG.coords[S].m + '  (p' + S + ')', fresh.has(S))).join('');
    h += `<div class="srow"><span>c<sub>${row.k}</sub>${row.k === D.n ? ' <span class="muted">det</span>' : ''}</span>
      <span class="chips">${chips}</span>${vtag(verdict(sg, 1))}</div>`;
  }
  const allP = SG.stab.every(r => r.items.every(S => minorSign(t, S) > 0));
  if (allP) h += '<div class="muted" style="font-size:12px">every principal minor is positive: B is a P-matrix</div>';
  // feasibility: x*_j det B = sum_k (-1)^(k+j) a_k det B_{[n]\k,[n]\j}, and x*_j > 0
  const dB = minorSign(t, SG.det.S);
  const target = dB || 1;
  h += `<h3 class="math">Feasibility · x*<sub>j</sub> det B = Σ<sub>k</sub> ± a<sub>k</sub> det B<sub>[n]∖k,[n]∖j</sub></h3>
    <div class="muted" style="font-size:12px;margin-bottom:4px">each sum must have the sign of det B
    (${dB ? (dB > 0 ? '+' : '−') : 'open, + assumed as stability requires'})</div>`;
  for (const row of SG.feas) {
    const sg = row.terms.map(x => t[x.S] ? x.f * t[x.S] : 0);
    const chips = row.terms.map((x, i) => chip(sg[i], `term k=${x.k}: a${x.k} · ${SG.coords[x.S].m}  (p${x.S})`, fresh.has(x.S))).join('');
    h += `<div class="srow"><span>x*<sub>${row.j}</sub></span><span class="chips">${chips}</span>${vtag(verdict(sg, target))}</div>`;
  }
  // the whole sign table, by size of the minor
  const bySize = {};
  for (const [S, info] of Object.entries(SG.coords)) (bySize[info.k] = bySize[info.k] || []).push(S);
  const known = Object.values(t).filter(v => v).length, total = Object.keys(t).length;
  h += `<details><summary>All ${total} minors · ${known} signs known, ${total - known} open</summary>`;
  for (const k of Object.keys(bySize).sort((x, y) => x - y)) {
    if (+k === 0) continue;
    const chips = bySize[k].map(S => chip(minorSign(t, S), SG.coords[S].m + '  (p' + S + ')', fresh.has(S))).join('');
    h += `<div class="srow"><span>${k}×${k}</span><span class="chips">${chips}</span><span></span></div>`;
  }
  h += '<div class="muted" style="font-size:12px">chips show the sign of the minor of B; a gold ring marks signs decided at this node</div></details>';
  return h;
}
function select(n) { selected = n; draw(); show(n); }

// ---------------------------------------------------------------- pan / zoom
let vb = { x: 0, y: 0, w: 100, h: 100 };
const setVB = () => svg.setAttribute('viewBox', `${vb.x} ${vb.y} ${vb.w} ${vb.h}`);
function fit() {
  const b = world.getBBox(), r = svg.getBoundingClientRect(), pad = 40;
  const s = Math.max((b.width + 2 * pad) / r.width, (b.height + 2 * pad) / r.height, 0.35);
  vb = { w: r.width * s, h: r.height * s, x: 0, y: 0 };
  vb.x = b.x + b.width / 2 - vb.w / 2; vb.y = b.y - pad;
  setVB();
}
svg.addEventListener('wheel', e => {
  e.preventDefault();
  const r = svg.getBoundingClientRect(), k = e.deltaY > 0 ? 1.12 : 1 / 1.12;
  const mx = vb.x + (e.clientX - r.left) / r.width * vb.w, my = vb.y + (e.clientY - r.top) / r.height * vb.h;
  vb = { x: mx - (mx - vb.x) * k, y: my - (my - vb.y) * k, w: vb.w * k, h: vb.h * k }; setVB();
}, { passive: false });
let drag = null;
const view = document.getElementById('view');
svg.addEventListener('pointerdown', e => { drag = { x: e.clientX, y: e.clientY, vb: { ...vb }, moved: false }; view.classList.add('drag'); });
window.addEventListener('pointermove', e => {
  if (!drag) return; const r = svg.getBoundingClientRect();
  const dx = (e.clientX - drag.x) / r.width * vb.w, dy = (e.clientY - drag.y) / r.height * vb.h;
  if (Math.abs(e.clientX - drag.x) + Math.abs(e.clientY - drag.y) > 3) drag.moved = true;
  vb.x = drag.vb.x - dx; vb.y = drag.vb.y - dy; setVB();
});
window.addEventListener('pointerup', e => {
  view.classList.remove('drag');
  if (drag && !drag.moved) { const g = e.target.closest('[data-id]'); if (g) select(byId.get(+g.dataset.id)); }
  drag = null;
});
svg.addEventListener('dblclick', e => { const g = e.target.closest('[data-id]'); if (g) toggle(byId.get(+g.dataset.id)); });
function toggle(n) { if (!n.kids.length) return; folded.has(n.id) ? folded.delete(n.id) : folded.add(n.id); draw(); show(n); }
// space: unfold everything below n, however many folded nodes lie in between; if nothing below is folded, fold n
function unfoldBelow(n) {
  if (!n.kids.length) return;
  const below = [], stack = [n];
  while (stack.length) { const m = stack.pop(); if (folded.has(m.id)) below.push(m.id); stack.push(...m.kids); }
  if (below.length) below.forEach(id => folded.delete(id)); else folded.add(n.id);
  document.getElementById('all').textContent = folded.size ? 'unfold all' : 'fold deep';
  draw(); show(n);
}

// ---------------------------------------------------------------- keys, search, theme
document.getElementById('fold').onclick = () => toggle(selected);
document.getElementById('fit').onclick = fit;
document.getElementById('theme').onclick = () => {
  const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
  root.dataset.theme = dark ? 'light' : 'dark';
};
document.getElementById('q').addEventListener('input', e => { query = e.target.value.replace(/[^0-9]/g, ''); draw(); });
window.addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT') return;
  const vis = selected.visibleKids || [];
  if (e.key === 'f') fit();
  else if (e.key === ' ') unfoldBelow(selected);
  else if (e.key === 'ArrowUp' && selected.parent !== null && selected.parent !== undefined) select(byId.get(selected.parent));
  else if (e.key === 'ArrowLeft' && vis[0]) select(vis[0]);
  else if (e.key === 'ArrowRight' && vis.length) select(vis[vis.length - 1]);
  else if (e.key === 'ArrowDown' && vis[0]) select(vis[0]);
  else return;
  e.preventDefault();
});
// large trees start folded at the deepest level that keeps about 400 nodes visible
function foldToBudget(budget) {
  folded.clear();
  if (D.nodes.length <= budget) return;
  const perDepth = [];
  for (const n of D.nodes) perDepth[n.depth] = (perDepth[n.depth] || 0) + 1;
  let total = 0, d = 0;
  while (d < perDepth.length && total + perDepth[d] <= budget) total += perDepth[d++];
  for (const n of D.nodes) if (n.depth === Math.max(d - 1, 0) && n.kids.length) folded.add(n.id);
}
foldToBudget(400);
document.getElementById('all').onclick = () => {
  if (folded.size) folded.clear(); else foldToBudget(400);
  document.getElementById('all').textContent = folded.size ? 'unfold all' : 'fold deep';
  draw(); fit();
};
document.getElementById('all').textContent = folded.size ? 'unfold all' : 'fold deep';
draw(); show(selected); fit();
window.addEventListener('resize', fit);
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
