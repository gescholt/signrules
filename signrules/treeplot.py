"""
The sign-completion search tree drawn without Graphviz (for large trees).

``complete(eco, record_tree=True)`` records every node of the binary search tree in
``result.tree``.  This module lays the tree out as a tidy binary tree (leaves evenly spaced,
each internal node above the middle of its children) and draws it with

* ``plot_search_tree`` — matplotlib, for figures and PNG/PDF files;
* ``search_tree_plotly`` — plotly, interactive HTML with hover details (needs ``plotly``).

Same conventions as ``certificate.search_tree_diagram``: an internal node shows the Plücker
coordinate it branches on, the edge to a child is red when the sign chosen was + and blue when
it was −, leaves are completions (green), contradictions (black ⊥) or unexplored subtrees
(grey, when a node or completion limit stopped the search).  Nothing else is printed on the
tree; the counts (open coordinates, forced signs) are in the plotly hover text.

    result = complete(eco, use=("R1", "R2n"), record_tree=True)
    plot_search_tree(eco, result); plt.savefig("tree.png")
    save_search_tree_html(eco, result, "tree.html")
"""

from typing import Dict, List, Optional, Tuple

import networkx as nx

from .rules import Ecology
from .search import SearchResult
from . import labels as L
from . import style
from .certificate import tree_title, repeated_subtrees, subtree_ids, _edge_width, _open_shade


# --------------------------------------------------------------------------- graph and layout

def tree_graph(result: SearchResult) -> nx.DiGraph:
    """The recorded search tree as a networkx DiGraph; node attributes are the recorded dicts."""
    if not result.tree:
        raise ValueError("no tree: run complete(..., record_tree=True)")
    G = nx.DiGraph()
    for nd in result.tree:
        G.add_node(nd["id"], **nd)
        if nd["parent"] is not None:
            G.add_edge(nd["parent"], nd["id"], sign=nd["sign"], forced=nd.get("forced") or 0)
    return G


def tree_layout(result: SearchResult, stopped_width: float = 0.45) -> Dict[int, Tuple[float, float]]:
    """
    Tidy layout: leaves at increasing x in search order, parents centred over their children,
    y = −depth.  Unexplored leaves (search stopped) take only ``stopped_width`` of a slot, so a
    long run of them does not stretch the picture.
    """
    children: Dict[Optional[int], List[int]] = {}
    status = {}
    for nd in result.tree:
        children.setdefault(nd["parent"], []).append(nd["id"])
        status[nd["id"]] = nd["status"]
    pos: Dict[int, Tuple[float, float]] = {}
    next_x = [0.0]

    def place(i, depth):
        kids = children.get(i, [])
        if not kids:
            width = stopped_width if status[i] == "stopped" else 1.0
            x = next_x[0] + width / 2
            next_x[0] += width
        else:
            xs = [place(k, depth + 1) for k in kids]
            x = sum(xs) / len(xs)
        pos[i] = (x, -float(depth))
        return x

    place(0, 0)
    return pos


def _name(eco: Ecology, S, plucker: bool, fmt: str) -> str:
    if S is None:
        return ""
    return L.plucker(S, fmt=fmt) if plucker else L.minor(eco, S, fmt=fmt)


# --------------------------------------------------------------------------- matplotlib

def plot_search_tree(eco: Ecology, result: SearchResult, ax=None, labels="auto", plucker: bool = True,
                     title=None, use=None, node_size: Optional[float] = None, legend: bool = True,
                     width_by_forced: bool = True, repeats: int = 2, shade_open: bool = False):
    """
    Draw the search tree with matplotlib.  ``labels``: True / False / "auto" (off above 300 nodes).
    Returns the axes.  Figure size is chosen from the tree's width and depth when ``ax`` is None.
    Colour channels as in ``certificate.search_tree_diagram``: edge colour = sign chosen, edge width
    = signs forced, ⊥ leaves in the colour of their rule, pale hulls around repeated subtrees
    (contiguous in this layout), optional grey shading by open coordinates.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch, Patch

    style.use_math_fonts()
    G = tree_graph(result)
    pos = tree_layout(result)
    width = max(x for x, _ in pos.values()) + 1
    depth = max(-y for _, y in pos.values())
    if labels == "auto":
        labels = len(G) <= 300
    if ax is None:
        w = min(60.0, max(6.0, 0.42 * width + 1.5))
        h = min(40.0, max(3.5, 0.7 * depth + 2.0))
        _, ax = plt.subplots(figsize=(w, h))
    size = node_size or max(60.0, min(900.0, 60000.0 / max(len(G), 1)))

    # hulls around repeated subtrees: subtrees are contiguous in x in this layout
    classes = repeated_subtrees(result, max_classes=repeats) if repeats else []
    node_r = 0.5 * (size ** 0.5) / 30.0                       # rough node radius in data units
    for k, (ids, copies, sz, root) in enumerate(classes):
        for r in ids:
            pts = [pos[j] for j in subtree_ids(result, r)]
            x0, x1 = min(x for x, _ in pts) - 0.45, max(x for x, _ in pts) + 0.45
            y0, y1 = min(y for _, y in pts) - 0.4, max(y for _, y in pts) + 0.4
            ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=0.15",
                                        facecolor=style.HULLS[k % len(style.HULLS)], edgecolor="none", zorder=0))
    for u, v, d in G.edges(data=True):
        (xu, yu), (xv, yv) = pos[u], pos[v]
        ax.plot([xu, xv], [yu, yv], color=style.sign_line(d["sign"]), zorder=1,
                lw=_edge_width(d.get("forced"), base=0.9, per_sign=0.12, cap=3.0) if width_by_forced else 1.1)
    root_open = result.tree[0].get("unknown") or 0
    looks = {"root": dict(node_color="white", edgecolors=style.INK, linewidths=1.6, node_shape="o"),
             "branch": dict(node_color="white", edgecolors=style.GREY, linewidths=0.9, node_shape="o"),
             "completion": dict(node_color=style.COMPLETION, edgecolors=style.COMPLETION, linewidths=1.0, node_shape="o"),
             "stopped": dict(node_color="white", edgecolors=style.STOPPED, linewidths=1.0, node_shape="o")}
    for status, look in looks.items():
        nodes = [i for i, d in G.nodes(data=True) if d["status"] == status]
        if nodes:
            s = size * (0.55 if status in ("completion", "stopped") else 1.0)
            if status == "branch" and shade_open:
                look = dict(look, node_color=[_open_shade(G.nodes[i], root_open) for i in nodes])
            nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax, node_size=s, **look)
    rules_seen = {}
    for i, d in G.nodes(data=True):
        if d["status"] == "conflict":
            lab = L.rule(d.get("reason", ""), eco.n)
            rules_seen[lab.rule] = lab.headline()
            nx.draw_networkx_nodes(G, pos, nodelist=[i], ax=ax, node_size=size * 0.55, node_shape="s",
                                   node_color=style.rule_color(lab.rule), edgecolors=style.rule_color(lab.rule), linewidths=1.0)
    if labels:
        fs = 9 if len(G) <= 60 else 7 if len(G) <= 160 else 5.5
        for i, d in G.nodes(data=True):
            x, y = pos[i]
            if d["status"] in ("root", "branch") and d.get("branch_on") is not None:
                ax.annotate(_name(eco, d["branch_on"], plucker, "tex"), (x, y), ha="center", va="center", fontsize=fs)
            elif d["status"] == "completion":
                ax.annotate(str(d.get("completion")), (x, y), ha="center", va="center", fontsize=fs * 0.75, color="white")
            elif d["status"] == "conflict":
                ax.annotate(r"$\bot$", (x, y), ha="center", va="center", fontsize=fs, color="white")
            elif d["status"] == "stopped":
                ax.annotate("…", (x, y), ha="center", va="center", fontsize=fs * 0.8, color=style.GREY)
    ax.set_axis_off()
    ax.margins(x=0.02, y=0.06)
    if title is None:
        title = tree_title(eco, result, use)
    if isinstance(title, (list, tuple)):
        head = title[0] if isinstance(title[0], str) else title[0][0]
        rest = [t if isinstance(t, str) else t[0] for t in title[1:]]
        ax.set_title(head + ("\n" + "\n".join(rest) if rest else ""), fontsize=11)
    else:
        ax.set_title(title, fontsize=11)
    if legend:
        handles = [Line2D([], [], color=style.POS_LINE, lw=2, label="chose +"),
                   Line2D([], [], color=style.NEG_LINE, lw=2, label="chose −")]
        if width_by_forced:
            handles.append(Line2D([], [], color=style.GREY, lw=3, label="thicker: more signs forced"))
        handles.append(Line2D([], [], marker="o", color="none", markerfacecolor=style.COMPLETION, markersize=8, label="completion"))
        for rule, head in sorted(rules_seen.items()):
            handles.append(Line2D([], [], marker="s", color="none", markerfacecolor=style.rule_color(rule), markersize=8,
                                  label=f"⊥ by {head}"))
        if any(d["status"] == "stopped" for _, d in G.nodes(data=True)):
            handles.append(Line2D([], [], marker="o", color="none", markerfacecolor="white", markeredgecolor=style.STOPPED,
                                  markersize=8, label="not searched"))
        for k, (ids, copies, sz, root) in enumerate(classes):
            handles.append(Patch(facecolor=style.HULLS[k % len(style.HULLS)], edgecolor="none",
                                 label=f"same subtree searched {copies}×  (from {L.plucker(root)}, {sz} nodes)"))
        if shade_open:
            handles.append(Line2D([], [], marker="o", color="none", markerfacecolor="#a0a0a0", markeredgecolor=style.GREY,
                                  markersize=8, label="darker: fewer coordinates open"))
        ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False, ncol=min(len(handles), 5))
    return ax


def save_search_tree_png(eco: Ecology, result: SearchResult, path: str, dpi: int = 130, **kw) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax = plot_search_tree(eco, result, **kw)
    ax.figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(ax.figure)


# --------------------------------------------------------------------------- plotly

def search_tree_plotly(eco: Ecology, result: SearchResult, plucker: bool = True, title=None, use=None, labels="auto"):
    """Interactive plotly figure of the search tree; hover shows the node's full record."""
    import plotly.graph_objects as go

    G = tree_graph(result)
    pos = tree_layout(result)
    if labels == "auto":
        labels = len(G) <= 400
    fig = go.Figure()
    for sign, name in ((1, "chose +"), (-1, "chose −")):
        ex, ey = [], []
        for u, v, d in G.edges(data=True):
            if d["sign"] == sign:
                ex += [pos[u][0], pos[v][0], None]
                ey += [pos[u][1], pos[v][1], None]
        fig.add_trace(go.Scatter(x=ex, y=ey, mode="lines", name=name, line=dict(color=style.sign_line(sign), width=1.2),
                                 hoverinfo="skip"))
    fills = {"root": "white", "branch": "white", "completion": style.COMPLETION, "conflict": style.BOTTOM, "stopped": "white"}
    borders = {"root": style.INK, "branch": style.GREY, "completion": style.COMPLETION, "conflict": style.BOTTOM, "stopped": style.STOPPED}
    for status, name in (("branch", "branch node"), ("root", "root"), ("completion", "completion"),
                         ("conflict", "contradiction"), ("stopped", "not searched")):
        nodes = [i for i, d in G.nodes(data=True) if d["status"] == status]
        if not nodes:
            continue
        hover = []
        for i in nodes:
            d = G.nodes[i]
            parts = [f"<b>{name}</b> (node {i}, depth {d['depth']})"]
            if d["basis"] is not None:
                parts.append(f"chose {_name(eco, d['basis'], plucker, 'html')} {'&gt; 0' if d['sign'] > 0 else '&lt; 0'}")
            if d.get("forced"):
                parts.append(f"forced by propagation: {d['forced']}")
            if d.get("unknown") is not None:
                parts.append(f"open coordinates: {d['unknown']}")
            if d.get("branch_on") is not None:
                parts.append(f"branches on {_name(eco, d['branch_on'], plucker, 'html')}")
            if status == "completion":
                parts.append(f"completion {d.get('completion')}")
            if status == "conflict":
                parts.append("⊥ " + L.rule(d.get("reason", ""), eco.n).headline("html"))
            hover.append("<br>".join(parts))
        text = None
        if labels:
            text = [(_name(eco, G.nodes[i]["branch_on"], plucker, "html") if status in ("root", "branch")
                     else str(G.nodes[i].get("completion")) if status == "completion"
                     else "⊥" if status == "conflict" else "…") for i in nodes]
        fig.add_trace(go.Scatter(
            x=[pos[i][0] for i in nodes], y=[pos[i][1] for i in nodes], name=name,
            mode="markers+text" if labels else "markers", text=text, textposition="middle center",
            textfont=dict(size=9 if len(G) <= 100 else 7, color="white" if status in ("completion", "conflict") else style.INK),
            marker=dict(size=(22 if len(G) <= 100 else 14 if len(G) <= 400 else 7) * (0.7 if status in ("completion", "conflict", "stopped") else 1),
                        color=fills[status], line=dict(color=borders[status], width=1.5),
                        symbol="square" if status == "conflict" else "circle"),
            hovertext=hover, hoverinfo="text"))
    if title is None:
        title = tree_title(eco, result, use)
    if isinstance(title, (list, tuple)):
        title = "<br>".join(L.convert(t if isinstance(t, str) else t[0], "html") for t in title)
    fig.update_layout(title=title, showlegend=True, plot_bgcolor="white",
                      xaxis=dict(visible=False), yaxis=dict(visible=False),
                      hovermode="closest", margin=dict(l=20, r=20, t=80, b=20))
    return fig


def save_search_tree_html(eco: Ecology, result: SearchResult, path: str, **kw) -> None:
    fig = search_tree_plotly(eco, result, **kw)
    fig.write_html(path, include_plotlyjs="cdn")
