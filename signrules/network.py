"""
The interaction network of a sign pattern, drawn as the article's ``\\econet`` macro.

Conventions (``fig:legend`` of the article in feasible_ecologies): a pink node is free (a_i > 0), a blue
node obligate (a_i < 0); a benefit of i from j (b_ij < 0) is a blue dashed arrow i → j; a pair that loses
both ways (b_ij, b_ji > 0) is a thin red line; a one-way loss is not drawn.  Species 0 is on top and the
others follow clockwise, as in the macro.

* ``econet_edges(a, B)`` — the edges, as the macro's ``i/j/b`` and ``i/j/c``;
* ``econet_tex(a, B)`` — the macro call, ready to paste into the article;
* ``draw_network(ax, a, B)`` / ``save_network(a, B, stem)`` — the matplotlib drawing (png, pdf, and the
  ``.tex`` call);
* ``network_text(a, B, width, height)`` — the same network in a terminal: a braille canvas as a
  ``rich.text.Text`` (used by ``signrules edit``).
"""
import math
from typing import List, Optional, Sequence, Tuple

from . import style

Edge = Tuple[int, int, str]


def econet_edges(a: Sequence[int], B) -> List[Edge]:
    """(i, j, "b") for every benefit b_ij < 0 (an arrow i → j), (i, j, "c") for every pair i < j that loses both ways."""
    n = len(a)
    out: List[Edge] = []
    for i in range(n):
        for j in range(n):
            if i != j and B[i][j] < 0:
                out.append((i, j, "b"))
    for i in range(n):
        for j in range(i + 1, n):
            if B[i][j] > 0 and B[j][i] > 0:
                out.append((i, j, "c"))
    return out


def econet_tex(a: Sequence[int], B) -> str:
    """The article's ``\\econet[n]{signs}{edges}`` call for this pattern."""
    n = len(a)
    signs = ",".join("p" if v > 0 else "m" for v in a)
    edges = ",".join(f"{i}/{j}/{t}" for i, j, t in econet_edges(a, B))
    opt = "" if n == 3 else f"[{n}]"
    return f"\\econet{opt}{{{signs}}}{{{edges}}}"


def layout(n: int) -> List[Tuple[float, float]]:
    """Node positions of ``\\econet``: angle 90° − 360° i/n (plus 90° when n = 2), radius 0.6 for n ≤ 3, else 1."""
    rad = 1.0 if n > 3 else 0.6
    out = []
    for i in range(n):
        ang = math.radians(90 - 360 * i / n + (90 if n == 2 else 0))
        out.append((rad * math.cos(ang), rad * math.sin(ang)))
    return out


# --------------------------------------------------------------------------- matplotlib

def draw_network(ax, a: Sequence[int], B, node_size: float = 0.13, labels: Optional[Sequence[str]] = None) -> None:
    """Draw the network on ``ax`` (equal aspect, no axes)."""
    from matplotlib.patches import Circle, FancyArrowPatch
    n = len(a)
    pos = layout(n)
    r = node_size * (1.0 if n > 3 else 0.8)
    nodes = []
    for i, (x, y) in enumerate(pos):
        c = Circle((x, y), r, facecolor=style.FREE_BG if a[i] > 0 else style.OBLIGATE_BG,
                   edgecolor=style.NODE_EDGE, lw=0.9, zorder=3)
        ax.add_patch(c)
        nodes.append(c)
        ax.text(x, y, labels[i] if labels else str(i), ha="center", va="center", fontsize=11, zorder=4)
    for i, j, t in econet_edges(a, B):                    # bend left 18°, as \econet; clipped at the nodes
        if t == "c":
            look = dict(arrowstyle="-", color=style.LOSS_LINE, lw=0.9, zorder=1)
        else:
            look = dict(arrowstyle="-|>,head_length=0.5,head_width=0.22", color=style.BENEFIT, lw=1.3,
                        linestyle=(0, (4, 2.5)), mutation_scale=12, zorder=2)
        ax.add_patch(FancyArrowPatch(pos[i], pos[j], connectionstyle="arc3,rad=0.16", patchA=nodes[i],
                                     patchB=nodes[j], shrinkA=1.5, shrinkB=1.5, **look))
    pad = 0.35
    ax.set_xlim(-1.0 - pad if n > 3 else -0.6 - pad, 1.0 + pad if n > 3 else 0.6 + pad)
    ax.set_ylim(-1.0 - pad if n > 3 else -0.6 - pad, 1.0 + pad if n > 3 else 0.6 + pad)
    ax.set_aspect("equal")
    ax.axis("off")


def save_network(a: Sequence[int], B, stem: str, title: str = "") -> List[str]:
    """Write ``stem``.png, ``stem``.pdf and ``stem``.tex (the ``\\econet`` call); return the paths."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(2.6, 2.6))
    draw_network(ax, a, B)
    if title:
        ax.set_title(title, fontsize=9)
    paths = []
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=200, bbox_inches="tight")
        paths.append(f"{stem}.{ext}")
    plt.close(fig)
    with open(f"{stem}.tex", "w") as f:
        f.write(econet_tex(a, B) + "\n")
    paths.append(f"{stem}.tex")
    return paths


# --------------------------------------------------------------------------- terminal (braille)

_DOTS = ((0x01, 0x08), (0x02, 0x10), (0x04, 0x20), (0x40, 0x80))   # [row][col] bit of a braille cell
_ARROWS = "→↗↑↖←↙↓↘"


class _Canvas:
    """A braille canvas: ``w × h`` cells of 2 × 4 dots, one colour per cell (the last one drawn)."""

    def __init__(self, w: int, h: int):
        self.w, self.h = w, h
        self.bits = [[0] * w for _ in range(h)]
        self.color = [[None] * w for _ in range(h)]
        self.glyph = {}                                   # (row, col) -> (text, style): overrides the dots

    def dot(self, x: float, y: float, color: str) -> None:
        xi, yi = int(round(x)), int(round(y))
        c, r = xi // 2, yi // 4
        if 0 <= c < self.w and 0 <= r < self.h and (r, c) not in self.glyph:
            self.bits[r][c] |= _DOTS[yi % 4][xi % 2]
            self.color[r][c] = color

    def line(self, p, q, color: str, dash: Optional[Tuple[int, int]] = None) -> None:
        (x0, y0), (x1, y1) = p, q
        steps = max(int(2 * max(abs(x1 - x0), abs(y1 - y0))), 1)
        for k in range(steps + 1):
            if dash and (k // 2) % sum(dash) >= dash[0]:
                continue
            t = k / steps
            self.dot(x0 + t * (x1 - x0), y0 + t * (y1 - y0), color)

    def put(self, row: int, col: int, text: str, st: str) -> None:
        for k, ch in enumerate(text):
            if 0 <= row < self.h and 0 <= col + k < self.w:
                self.glyph[(row, col + k)] = (ch, st)

    def text(self):
        from rich.text import Text
        out = Text()
        for r in range(self.h):
            for c in range(self.w):
                if (r, c) in self.glyph:
                    ch, st = self.glyph[(r, c)]
                    out.append(ch, style=st)
                elif self.bits[r][c]:
                    out.append(chr(0x2800 + self.bits[r][c]), style=self.color[r][c])
                else:
                    out.append(" ")
            if r < self.h - 1:
                out.append("\n")
        return out


def network_text(a: Sequence[int], B, width: int = 40, height: int = 16, legend: bool = True):
    """The network as a ``rich.text.Text`` of ``width × height`` cells (plus two legend lines)."""
    from rich.text import Text
    n = len(a)
    cv = _Canvas(width, height)
    W, H = 2 * width, 4 * height                          # dots; a dot is about square in a terminal
    R = 0.5 * min(W, H) - 6
    cx, cy = W / 2, H / 2
    unit = max(abs(x) for p in layout(max(n, 4)) for x in p) or 1.0
    pos = [(cx + R * x / unit, cy - R * y / unit) for x, y in layout(n)]
    rx, ry = 5.0, 4.5                                     # the node label's half size in dots

    def clip(p, q):
        """The point where the segment p → q leaves the node ellipse around q."""
        dx, dy = p[0] - q[0], p[1] - q[1]
        d = math.hypot(dx / rx, dy / ry) or 1.0
        return q[0] + dx / d, q[1] + dy / d

    edges = econet_edges(a, B)
    both = {(i, j) for i, j, t in edges if t == "b" and (j, i, "b") in edges}
    for i, j, t in edges:
        p, q = pos[i], pos[j]
        if t == "c":
            cv.line(clip(q, p), clip(p, q), style.LOSS)
            continue
        if (i, j) in both:                                # two arrows: shift each to its own side
            dx, dy = q[0] - p[0], q[1] - p[1]
            L = math.hypot(dx, dy) or 1.0
            ox, oy = -dy / L * 1.6, dx / L * 1.6
            p, q = (p[0] + ox, p[1] + oy), (q[0] + ox, q[1] + oy)
        s, e = clip(q, p), clip(p, q)
        cv.line(s, e, style.BENEFIT, dash=(3, 2))
        ang = math.atan2(-(e[1] - s[1]), e[0] - s[0])
        k = int(round(ang / (math.pi / 4))) % 8
        cv.put(int(e[1]) // 4, int(e[0]) // 2, _ARROWS[k], f"bold {style.BENEFIT}")
    for i, (x, y) in enumerate(pos):
        lab = f" {i} "
        bg, fg = (style.FREE_BG, style.FREE_FG) if a[i] > 0 else (style.OBLIGATE_BG, style.OBLIGATE_FG)
        cv.put(int(y) // 4, int(x) // 2 - len(lab) // 2, lab, f"bold {fg} on {bg}")
    out = cv.text()
    if legend:
        out.append("\n")
        out.append(" 0 ", style=f"bold {style.FREE_FG} on {style.FREE_BG}"); out.append(" free  ")
        out.append(" 0 ", style=f"bold {style.OBLIGATE_FG} on {style.OBLIGATE_BG}"); out.append(" obligate\n")
        out.append("⠒⠒→", style=style.BENEFIT); out.append(" i benefits from j  ")
        out.append("⠤⠤⠤", style=style.LOSS); out.append(" both lose")
    return out
