"""
Graphviz layout, matplotlib rendering — so that graph labels can be mathematics.

Graphviz cannot typeset ``$p_{0257} > 0$``; matplotlib can, but has no layered DAG layout.  A
``Diagram`` takes the best of both: nodes are measured with matplotlib's mathtext engine, passed
to ``dot`` as fixed-size empty boxes, and the returned positions (``dot -Tjson``: nodes,
cluster bounding boxes, edge splines) are drawn with matplotlib.

    d = Diagram(title="$n=4$")
    d.node("a", ["$p_{0257} > 0$"], fill="#f6d3cd")
    d.node("r", ["R4 · Z-block", "$Q=\\{1,2,3\\}$"], shape="box", size=8)
    d.edge("a", "r")
    d.save("figure")            # figure.png + figure.pdf (+ figure.dot with plain-text labels)

Labels are lists of lines; a line is a string or a ``(text, fontsize, colour)`` tuple.
"""

import json
import shutil
import subprocess
from typing import Dict, List, Optional, Sequence, Tuple, Union

from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle
from matplotlib.path import Path as MplPath

from . import labels as L
from . import style

PT = 72.0                                     # points per inch
Line = Union[str, Tuple[str, float, str]]


class _Measure:
    """Text extents in points, from a hidden figure; cached."""

    def __init__(self):
        self._fig = None
        self._cache: Dict[Tuple[str, float], Tuple[float, float]] = {}

    def __call__(self, text: str, size: float) -> Tuple[float, float]:
        key = (text, size)
        if key not in self._cache:
            if self._fig is None:
                import matplotlib.pyplot as plt
                style.use_math_fonts()
                self._fig = plt.figure(figsize=(1, 1), dpi=100)
                self._renderer = self._fig.canvas.get_renderer()
            t = self._fig.text(0, 0, text or " ", fontsize=size)
            bb = t.get_window_extent(renderer=self._renderer)
            t.remove()
            # window extents are in display pixels at dpi=100
            self._cache[key] = (bb.width * PT / 100.0, max(bb.height * PT / 100.0, size * 0.9))
        return self._cache[key]


measure = _Measure()


def _lines(label, size: float, color: str) -> List[Tuple[str, float, str]]:
    if label is None:
        return []
    if isinstance(label, str):
        label = [label]
    out = []
    for ln in label:
        if isinstance(ln, str):
            for piece in ln.split("\n"):
                out.append((piece, size, color))
        else:
            text, s, c = ln
            for piece in text.split("\n"):
                out.append((piece, s or size, c or color))
    return out


class Diagram:
    """A graph of boxed nodes with mathtext labels: laid out by Graphviz (``dot -Tjson``), drawn with matplotlib.

    Build it with ``node``, ``edge``, ``cluster`` and ``legend_item``, then ``save(stem)`` writes ``stem.png``,
    ``stem.pdf`` and a plain-text ``stem.dot``.
    """

    def __init__(self, rankdir: str = "TB", nodesep: float = 0.22, ranksep: float = 0.4,
                 fontsize: float = 10, title: Optional[Line] = None, title_size: float = 11):
        self.rankdir, self.nodesep, self.ranksep, self.fontsize = rankdir, nodesep, ranksep, fontsize
        self.title, self.title_size = title, title_size
        self.nodes: Dict[str, dict] = {}
        self.edges: List[Tuple[str, str, dict]] = []
        self.clusters: Dict[str, dict] = {}
        self.ranks: List[List[str]] = []
        self.legend: List[Tuple[str, dict]] = []
        self.hulls: List[dict] = []
        self._layout = None

    # ------------------------------------------------------------------ building

    def node(self, nid: str, label, *, fill: str = "white", edge: str = style.GREY, lw: float = 1.0,
             ls: str = "-", shape: str = "round", size: Optional[float] = None, color: str = style.INK,
             pad: Tuple[float, float] = (7.0, 4.0), cluster: Optional[str] = None, gap: float = 2.5) -> str:
        """Add a node; ``shape`` is 'round' (rounded box), 'box', 'circle' or 'plain' (text only)."""
        size = size or self.fontsize
        lines = _lines(label, size, color)
        ws, hs = zip(*(measure(t, s) for t, s, _ in lines)) if lines else ((0,), (0,))
        w = max(ws) + 2 * pad[0]
        h = sum(hs) + gap * (len(lines) - 1) + 2 * pad[1]
        if shape == "circle":
            w = h = max(w, h)
        self.nodes[nid] = dict(lines=lines, fill=fill, edge=edge, lw=lw, ls=ls, shape=shape, w=w, h=h,
                               cluster=cluster, heights=hs, gap=gap)
        if cluster is not None:
            self.clusters.setdefault(cluster, dict(label=None, color=style.LIGHT, fill=None, nodes=[]))
            self.clusters[cluster]["nodes"].append(nid)
        return nid

    def edge(self, u: str, v: str, *, color: str = "#777777", lw: float = 0.9, ls: str = "-",
             arrow: bool = True, label: Optional[str] = None, size: Optional[float] = None,
             label_color: str = style.INK, constraint: bool = True, tailport: str = "", headport: str = "",
             ltail: str = "", lhead: str = "") -> None:
        """``ltail`` / ``lhead`` name a cluster: the edge is then clipped at the cluster's border."""
        self.edges.append((u, v, dict(color=color, lw=lw, ls=ls, arrow=arrow, label=label,
                                      size=size or self.fontsize * 0.85, label_color=label_color,
                                      constraint=constraint, tailport=tailport, headport=headport,
                                      ltail=ltail, lhead=lhead)))

    def cluster(self, cid: str, label: Optional[str] = None, color: str = style.LIGHT,
                fill: Optional[str] = None, size: Optional[float] = None, parent: Optional[str] = None,
                label_pos: str = "tl", lw: float = 0.9, margin: float = 8.0) -> None:
        """A Graphviz cluster (kept contiguous by the layout); ``parent`` nests it inside another one."""
        c = self.clusters.setdefault(cid, dict(nodes=[]))
        c.update(label=label, color=color, fill=fill, size=size or self.fontsize * 0.95, parent=parent,
                 label_pos=label_pos, lw=lw, margin=margin)

    def same_rank(self, nids: Sequence[str]) -> None:
        self.ranks.append(list(nids))

    def legend_item(self, text: str, **node_style) -> None:
        self.legend.append((text, node_style))

    # ------------------------------------------------------------------ layout

    def _dot_source(self, plain_labels: bool = False) -> str:
        """DOT for layout (empty fixed-size nodes) or, with ``plain_labels``, a human-readable export."""
        g = ['digraph G {',
             f'  graph [rankdir={self.rankdir}, nodesep={self.nodesep}, ranksep={self.ranksep}, '
             f'fontname="Helvetica", fontsize={self.fontsize}, dpi=72, compound=true];',
             f'  node [fontname="Helvetica", fontsize={self.fontsize}, style="rounded,filled"];',
             '  edge [arrowsize=0.6];']

        def node_line(nid: str, nd: dict, indent: str) -> str:
            shape = {"round": "box", "box": "box", "circle": "circle", "plain": "plaintext"}[nd["shape"]]
            if plain_labels:
                text = "\\n".join(L.convert(t, "text") for t, _, _ in nd["lines"]).replace('"', "'")
                return (f'{indent}"{nid}" [label="{text}", shape={shape}, fillcolor="{nd["fill"]}", '
                        f'color="{nd["edge"]}"];')
            return (f'{indent}"{nid}" [label="", shape={shape}, width={nd["w"] / PT:.3f}, '
                    f'height={nd["h"] / PT:.3f}, fixedsize=true];')

        clustered = {n: cid for cid, c in self.clusters.items() for n in c["nodes"]}

        def rank_line(r, indent):
            return indent + "{ rank=same; " + "; ".join(f'"{n}"' for n in r) + "; }"

        def emit_cluster(cid: str, indent: str):
            c = self.clusters[cid]
            g.append(f'{indent}subgraph "cluster_{cid}" {{')
            attrs = ['style="rounded"', f'color="{c.get("color", style.LIGHT)}"', f'margin={c.get("margin", 8.0)}']
            if c.get("label"):
                w, _ = measure(c["label"], c["size"])
                lab = L.convert(c["label"], "text").replace('"', "'") if plain_labels else "M" * max(1, int(w / (0.62 * c["size"])))
                attrs += [f'label="{lab}"', f'fontsize={c["size"]}', f'labeljust={"r" if c.get("label_pos", "tl").endswith("r") else "l"}']
            else:
                attrs.append('label=""')
            g.append(f'{indent}  ' + "; ".join(attrs) + ";")
            for nid in c["nodes"]:
                g.append(node_line(nid, self.nodes[nid], indent + "  "))
            for r in self.ranks:                       # rank constraints must sit inside their cluster
                if all(clustered.get(n) == cid for n in r):
                    g.append(rank_line(r, indent + "  "))
            for sub, sc in self.clusters.items():
                if sc.get("parent") == cid:
                    emit_cluster(sub, indent + "  ")
            g.append(f"{indent}}}")

        for cid, c in self.clusters.items():
            if c.get("parent") is None:
                emit_cluster(cid, "  ")
        for nid, nd in self.nodes.items():
            if nid not in clustered:
                g.append(node_line(nid, nd, "  "))
        for r in self.ranks:
            if not all(clustered.get(n) == clustered.get(r[0]) and clustered.get(n) is not None for n in r):
                g.append(rank_line(r, "  "))
        for u, v, e in self.edges:
            attrs = []
            if e["label"]:
                w, _ = measure(e["label"], e["size"])
                lab = L.convert(e["label"], "text") if plain_labels else "M" * max(1, int(w / (0.62 * e["size"])))
                attrs.append(f'label="{lab}", fontsize={e["size"]}')
            if not e["constraint"]:
                attrs.append("constraint=false")
            if e["ls"] == "invis":
                attrs.append("style=invis")
            if e["tailport"]:
                attrs.append(f'tailport={e["tailport"]}')
            if e["headport"]:
                attrs.append(f'headport={e["headport"]}')
            if e["ltail"]:
                attrs.append(f'ltail="cluster_{e["ltail"]}"')
            if e["lhead"]:
                attrs.append(f'lhead="cluster_{e["lhead"]}"')
            if plain_labels:
                attrs.append(f'color="{e["color"]}"')
                if not e["arrow"]:
                    attrs.append("arrowhead=none")
            g.append(f'  "{u}" -> "{v}"' + (f' [{", ".join(attrs)}]' if attrs else "") + ";")
        g.append("}")
        return "\n".join(g)

    def to_dot(self) -> str:
        """Graphviz source with plain-text labels."""
        return self._dot_source(plain_labels=True)

    def layout(self) -> dict:
        """Node positions and edge splines from ``dot -Tjson`` (cached)."""
        if self._layout is not None:
            return self._layout
        if shutil.which("dot") is None:
            raise RuntimeError("Graphviz 'dot' is required for the layout")
        out = subprocess.run(["dot", "-Tjson"], input=self._dot_source().encode(), capture_output=True, check=True)
        data = json.loads(out.stdout)
        pos, clusters, gvid = {}, {}, {}
        for obj in data.get("objects", []):
            gvid[obj["_gvid"]] = obj["name"]
            if obj["name"].startswith("cluster_"):
                clusters[obj["name"][8:]] = tuple(float(x) for x in obj["bb"].split(","))
            elif "pos" in obj:
                x, y = (float(v) for v in obj["pos"].split(","))
                pos[obj["name"]] = (x, y)
        edges = []
        for e in data.get("edges", []):
            edges.append((gvid[e["tail"]], gvid[e["head"]], e["pos"],
                          tuple(float(v) for v in e["lp"].split(",")) if "lp" in e else None))
        bb = tuple(float(x) for x in data["bb"].split(","))
        self._layout = dict(pos=pos, clusters=clusters, edges=edges, bb=bb)
        return self._layout

    # ------------------------------------------------------------------ drawing

    @staticmethod
    def _spline(pos: str) -> Tuple[MplPath, Tuple[float, float]]:
        """Graphviz edge 'pos' (spline with optional s,/e, end points) to a matplotlib Path."""
        toks = pos.split()
        start = end = None
        pts = []
        for t in toks:
            if t.startswith("e,"):
                end = tuple(float(v) for v in t[2:].split(","))
            elif t.startswith("s,"):
                start = tuple(float(v) for v in t[2:].split(","))
            else:
                pts.append(tuple(float(v) for v in t.split(",")))
        verts, codes = [], []
        if start is not None:
            verts += [start, pts[0]]
            codes += [MplPath.MOVETO, MplPath.LINETO]
        else:
            verts.append(pts[0])
            codes.append(MplPath.MOVETO)
        for k in range(1, len(pts) - 2, 3):
            verts += pts[k:k + 3]
            codes += [MplPath.CURVE4] * 3
        if end is not None:
            verts.append(end)
            codes.append(MplPath.LINETO)
        return MplPath(verts, codes), (end or verts[-1])

    def draw(self, dpi: int = 170):
        """Draw the laid-out diagram with matplotlib; returns the figure."""
        import matplotlib.pyplot as plt
        style.use_math_fonts()
        lay = self.layout()
        x0, y0, x1, y1 = lay["bb"]
        margin = 6.0
        title_h = 0.0
        if self.title:
            tl = _lines(self.title, self.title_size, style.INK)
            title_h = sum(measure(t, s)[1] for t, s, _ in tl) + 8.0
        legend_h = legend_w = 0.0
        if self.legend:
            fs = self.fontsize * 0.85
            legend_h = fs + 14.0
            legend_w = sum(fs * 1.6 + 4 + measure(t, fs)[0] + 18 for t, _ in self.legend)
        need_w = max(x1 - x0, legend_w,
                     max((measure(t, s)[0] for t, s, _ in _lines(self.title, self.title_size, style.INK)), default=0))
        extra = (need_w - (x1 - x0)) / 2
        W, H = need_w + 2 * margin, y1 - y0 + 2 * margin + title_h + legend_h
        fig = plt.figure(figsize=(W / PT, H / PT), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(x0 - margin - extra, x1 + margin + extra)
        ax.set_ylim(y0 - margin - legend_h, y1 + margin + title_h)
        ax.set_aspect("equal")
        ax.axis("off")

        for h in self.hulls:
            xs0 = [lay["pos"][n][0] - self.nodes[n]["w"] / 2 for n in h["nodes"]]
            xs1 = [lay["pos"][n][0] + self.nodes[n]["w"] / 2 for n in h["nodes"]]
            ys0 = [lay["pos"][n][1] - self.nodes[n]["h"] / 2 for n in h["nodes"]]
            ys1 = [lay["pos"][n][1] + self.nodes[n]["h"] / 2 for n in h["nodes"]]
            hx0, hx1, hy0, hy1 = min(xs0) - h["pad"], max(xs1) + h["pad"], min(ys0) - h["pad"], max(ys1) + h["pad"]
            ax.add_patch(FancyBboxPatch((hx0, hy0), hx1 - hx0, hy1 - hy0, boxstyle="round,pad=0,rounding_size=6",
                                        facecolor=h["fill"], edgecolor=h["edge"] or "none", lw=0.8, zorder=0))
            if h["label"]:
                ax.text(hx1 - 3, hy1 - 2, h["label"], ha="right", va="top", fontsize=h["size"], color=style.GREY, zorder=1)

        def depth(cid):
            k, c = 0, self.clusters[cid]
            while c.get("parent") is not None:
                k, c = k + 1, self.clusters[c["parent"]]
            return k

        for cid, (cx0, cy0, cx1, cy1) in sorted(lay["clusters"].items(), key=lambda kv: depth(kv[0])):
            c = self.clusters[cid]
            ax.add_patch(FancyBboxPatch((cx0, cy0), cx1 - cx0, cy1 - cy0, boxstyle="round,pad=0,rounding_size=6",
                                        facecolor=c.get("fill") or "none", edgecolor=c.get("color") or "none",
                                        lw=c.get("lw", 0.9), zorder=0.1 * depth(cid)))
            if c.get("label"):
                right = c.get("label_pos", "tl").endswith("r")
                ax.text(cx1 - 8 if right else cx0 + 8, cy1 - 4, c["label"], ha="right" if right else "left", va="top",
                        fontsize=c["size"], color=style.GREY if c.get("fill") else style.INK, zorder=1)

        drawn: Dict[Tuple[str, str], int] = {}
        for u, v, pos, lp in lay["edges"]:
            cands = [e for a, b, e in self.edges if a == u and b == v]
            e = cands[min(drawn.get((u, v), 0), len(cands) - 1)]
            drawn[(u, v)] = drawn.get((u, v), 0) + 1
            if e["ls"] == "invis":
                continue
            path, _ = self._spline(pos)
            patch = FancyArrowPatch(path=path, arrowstyle="-|>" if e["arrow"] else "-", mutation_scale=9,
                                    color=e["color"], lw=e["lw"], ls=e["ls"], shrinkA=0, shrinkB=0, zorder=2)
            ax.add_patch(patch)
            if e["label"] and lp:
                ax.text(lp[0], lp[1], e["label"], ha="center", va="center", fontsize=e["size"], color=e["label_color"],
                        zorder=4, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))

        for nid, nd in self.nodes.items():
            x, y = lay["pos"][nid]
            w, h = nd["w"], nd["h"]
            if nd["shape"] == "circle":
                ax.add_patch(Circle((x, y), w / 2, facecolor=nd["fill"], edgecolor=nd["edge"], lw=nd["lw"], ls=nd["ls"], zorder=3))
            elif nd["shape"] == "round":
                ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=5",
                                            facecolor=nd["fill"], edgecolor=nd["edge"], lw=nd["lw"], ls=nd["ls"], zorder=3))
            elif nd["shape"] == "box":
                ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, facecolor=nd["fill"], edgecolor=nd["edge"],
                                       lw=nd["lw"], ls=nd["ls"], zorder=3))
            total = sum(nd["heights"]) + nd["gap"] * (len(nd["lines"]) - 1)
            top = y + total / 2
            for (text, size, color), lh in zip(nd["lines"], nd["heights"]):
                ax.text(x, top - lh / 2, text, ha="center", va="center", fontsize=size, color=color, zorder=5)
                top -= lh + nd["gap"]

        if self.title:
            top = y1 + margin + title_h - 2
            for text, size, color in _lines(self.title, self.title_size, style.INK):
                ax.text((x0 + x1) / 2, top, text, ha="center", va="top", fontsize=size, color=color)
                top -= measure(text, size)[1] + 2
        if self.legend:
            x = x0 - margin - extra + 4
            yl = y0 - margin - legend_h / 2 - 2
            fs = self.fontsize * 0.85
            for text, st in self.legend:
                bw, bh = fs * 1.6, fs * 1.1
                if st.get("shape") == "circle":
                    ax.add_patch(Circle((x + bh / 2, yl), bh / 2, facecolor=st.get("fill", "white"),
                                        edgecolor=st.get("edge", style.GREY), lw=st.get("lw", 1.0), ls=st.get("ls", "-")))
                    x += bh
                else:
                    ax.add_patch(FancyBboxPatch((x, yl - bh / 2), bw, bh, boxstyle="round,pad=0,rounding_size=3",
                                                facecolor=st.get("fill", "white"), edgecolor=st.get("edge", style.GREY),
                                                lw=st.get("lw", 1.0), ls=st.get("ls", "-")))
                    x += bw
                ax.text(x + 4, yl, text, ha="left", va="center", fontsize=fs, color=style.GREY)
                x += measure(text, fs)[0] + 18
        return fig, ax

    def save(self, stem: str, formats: Sequence[str] = ("png", "pdf"), dpi: int = 170, dot: bool = True) -> None:
        """Write stem.png / stem.pdf (and stem.dot with plain-text labels)."""
        import matplotlib.pyplot as plt
        fig, _ = self.draw(dpi=dpi)
        for fmt in formats:
            fig.savefig(f"{stem}.{fmt}", dpi=dpi)
        plt.close(fig)
        if dot:
            with open(f"{stem}.dot", "w") as f:
                f.write(self.to_dot())
