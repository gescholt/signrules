"""
Illustrations of the sub-steps of the sign-rule algorithm (docs/SIGNRULES.md).

    python docs/figures/make_figures.py        # writes docs/figures/*.png / *.pdf

Each figure is drawn from an ``Ecology`` (or, for the dictionary, from a basis S): every set,
signature, permutation and formula printed on it is computed by the package, and all text is
typeset as mathematics through ``signrules.labels`` — nothing is typed by hand, so the same
function draws the figure for any other pattern.

  dictionary.png     p_S, a Plücker coordinate of [diag(a) | B], collapses to a minor of B
  tope.png           R7: the feasible sign vectors of a are cells of a hyperplane arrangement (n = 2)
  zblock.png         R4: a row-sign twist and a column permutation expose a Z-matrix
  balanced.png       R3: a balanced signed graph and the signature that makes S B S a Z-matrix

Colour follows docs/figures/README.md: red = positive entry / coordinate, blue = negative; the two
highlight tints mark subsets (chosen columns, a group of species) and never a sign.
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon

from signrules import Ecology
from signrules import labels as L
from signrules import style

OUT = Path(__file__).resolve().parent
style.use_math_fonts()


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=170, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def sign_cell(ax, x, y, v, outline=False, size=15):
    """One entry of a sign matrix: tint and glyph by sign; ``outline`` for a diagonal entry."""
    ax.add_patch(Rectangle((x, y), 1, 1, facecolor=style.sign_tint(v), edgecolor="#999999", lw=0.8))
    if outline:
        ax.add_patch(Rectangle((x + 0.07, y + 0.07), 0.86, 0.86, facecolor="none", edgecolor=style.INK, lw=1.4))
    ax.text(x + 0.5, y + 0.5, "$+$" if v > 0 else "$-$", ha="center", va="center", fontsize=size, color=style.sign_ink(v))


def sign_matrix(ax, M, title, row_labels=None, col_labels=None, outline_diag=False):
    """Draw sign(M) as a grid, rows top to bottom, with index labels and a math title."""
    n = M.shape[0]
    ax.set_xlim(-0.8, n + 0.3)
    ax.set_ylim(-0.4, n + 0.9)
    ax.set_aspect("equal")
    ax.axis("off")
    for i in range(n):
        for j in range(n):
            sign_cell(ax, j, n - 1 - i, M[i, j], outline=outline_diag and i == j)
    for i in range(n):
        ax.text(-0.4, n - 1 - i + 0.5, (row_labels or [str(k) for k in range(n)])[i], ha="center", va="center",
                fontsize=9, color=style.GREY)
        ax.text(i + 0.5, n + 0.1, (col_labels or [str(k) for k in range(n)])[i], ha="center", fontsize=9, color=style.GREY)
    ax.set_title(title, fontsize=12, pad=14)


# --------------------------------------------------------------------------- 1. dictionary

def fig_dictionary(eco, S, name="dictionary"):
    """The Plücker coordinate p_S of [diag(a) | B]: columns S of the n × 2n matrix, and the minor it collapses to."""
    n = eco.n
    T, P, Q, eps = eco.idx.decompose(tuple(S))
    fig, ax = plt.subplots(figsize=(1.35 * (2 * n + 2) + 3.2, 0.95 * n + 2.4))
    ax.set_xlim(-0.8, 2 * n + 1.2 + n + 3.2)
    ax.set_ylim(-2.2, n + 1.6)
    ax.set_aspect("equal")
    ax.axis("off")

    def cell(x, y, text, bg, fg=style.INK, bold=False, dim=False):
        ax.add_patch(Rectangle((x, y), 1, 1, facecolor=bg, edgecolor="#aaaaaa", lw=0.8))
        ax.text(x + 0.5, y + 0.5, text, ha="center", va="center", fontsize=11,
                color="#bbbbbb" if dim else fg, fontweight="bold" if bold else "normal")

    gap = 0.35
    # the n × 2n matrix, columns indexed 0 … 2n−1: that index set is where S lives
    for i in range(n):
        y = n - 1 - i
        for j in range(n):
            chosen = j in T
            cell(j, y, f"$a_{{{L._i(j)}}}$" if i == j else "$0$", style.HILITE_A if chosen else "white",
                 dim=(i != j), bold=chosen and i == j)
        for j in range(n):
            chosen = j in Q
            keep = chosen and i in P
            cell(n + gap + j, y, f"$b_{{{L._i(i)}{L._i(j)}}}$", style.HILITE_B if keep else "white",
                 dim=chosen and not keep, bold=keep)
    ax.plot([n + gap / 2, n + gap / 2], [-0.15, n + 0.15], color=style.INK, lw=1.0)
    for col in range(2 * n):                     # column indices of [diag(a) | B]: the ground set of Gr(n, 2n)
        x = col + (gap if col >= n else 0) + 0.5
        inS = col in S
        ax.text(x, n + 0.3, f"${L._i(col)}$", ha="center", va="bottom", fontsize=11,
                color=style.INK if inS else "#aaaaaa", fontweight="bold" if inS else "normal")
    ax.text(n / 2, n + 1.05, r"$\mathrm{diag}(a)$", ha="center", fontsize=12)
    ax.text(n + gap + n / 2, n + 1.05, "$B$", ha="center", fontsize=12)
    ax.text(n + gap / 2, -0.75, rf"$S = {L.tex_set(S)} = T \cup (n + Q)$,   $T={L.tex_set(T)}$,   "
            rf"$P = [n] \setminus T = {L.tex_set(P)}$,   $Q={L.tex_set(Q)}$", ha="center", fontsize=11)

    # the collapsed minor: rows P, columns Q of B
    x0 = 2 * n + gap + 1.6
    k = len(P)
    ytop = n / 2 + k / 2
    if k:
        for r, i in enumerate(P):
            for c, j in enumerate(Q):
                cell(x0 + 1.0 + c, ytop - 1 - r, f"$b_{{{L._i(i)}{L._i(j)}}}$", style.HILITE_B, bold=True)
        ax.text(x0 + 0.45, n / 2, r"$\det$", ha="center", va="center", fontsize=13)
        ax.text(x0 + 1.0 + k + 0.25, n / 2, rf"$= \det B_{{{L.tex_set(P)},{L.tex_set(Q)}}}$", ha="left", va="center", fontsize=13)
    else:
        ax.text(x0 + 0.45, n / 2, r"$1$", ha="center", va="center", fontsize=13)
    ax.annotate("", xy=(x0 - 0.2, n / 2), xytext=(2 * n + gap + 0.25, n / 2),
                arrowprops=dict(arrowstyle="-|>", color=style.INK, lw=1.0))
    ax.text(x0 + 0.2 + (k + 1.4) / 2, ytop - k - 0.6, "expand along the $a$-columns; rows $T$ are used up",
            ha="center", va="top", fontsize=9.5, color=style.GREY)
    ax.text((2 * n + gap + x0 + k + 2) / 2, -1.7, L.dictionary(eco, tuple(S)) + rf"   ($\varepsilon(S) = {'+' if eps > 0 else '-'}1$)",
            ha="center", fontsize=13)
    save(fig, name)


# --------------------------------------------------------------------------- 2. tope (R7)

def fig_tope(B, name="tope", L_box=3.0):
    """
    n = 2: the arrangement {x_1 = 0, x_2 = 0, (Bx)_1 = 0, (Bx)_2 = 0}.  Each cell carries the sign
    vector (sign x | sign Bx); the cells inside the positive quadrant are the feasible sign a.
    """
    B = np.asarray(B, dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 6.6))
    ax.set_xlim(-L_box, L_box)
    ax.set_ylim(-L_box, L_box)
    ax.set_aspect("equal")
    ax.add_patch(Rectangle((0, 0), L_box, L_box, facecolor=style.QUADRANT, edgecolor="none", zorder=0))
    ax.axhline(0, color=style.INK, lw=1.1)
    ax.axvline(0, color=style.INK, lw=1.1)
    i0, i1 = L._i(0), L._i(1)
    ax.text(L_box - 0.1, 0.12, f"$x_{{{i1}}} = 0$", ha="right", va="bottom", fontsize=10)
    ax.text(0.12, L_box - 0.1, f"$x_{{{i0}}} = 0$", ha="left", va="top", fontsize=10)
    t = np.linspace(-L_box, L_box, 2)
    for i, (b, ls) in enumerate(zip(B, ("-", "--"))):
        d = np.array([-b[1], b[0]]) / np.linalg.norm(b)
        ax.plot(t * d[0], t * d[1], color=style.GREY, lw=1.6, ls=ls)
        ax.annotate(rf"$(Bx)_{{{L._i(i)}}} = 0$", xy=(2.75 * d[0], 2.75 * d[1]), xytext=(-8, 10),
                    textcoords="offset points", color=style.GREY, fontsize=10, ha="right",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))
    # the cells: contiguous runs of directions with the same sign vector
    angles = np.linspace(0, 2 * np.pi, 1440, endpoint=False) + np.pi / 1440       # avoid the lines exactly
    svs = [tuple(np.sign(np.r_[np.cos(t_), np.sin(t_), B @ np.array([np.cos(t_), np.sin(t_)])]).astype(int)) for t_ in angles]
    runs = []
    for t_, sv in zip(angles, svs):
        if runs and runs[-1][0] == sv:
            runs[-1][1].append(t_)
        else:
            runs.append((sv, [t_]))
    if len(runs) > 1 and runs[0][0] == runs[-1][0]:
        runs[0] = (runs[0][0], [t_ - 2 * np.pi for t_ in runs[-1][1]] + runs[0][1])
        runs.pop()
    feasible = []
    for sv, angs in runs:
        mid = float(np.mean(angs))
        pos = 2.2 * np.array([np.cos(mid), np.sin(mid)])
        inpos = sv[0] > 0 and sv[1] > 0
        label = "$(" + r"\," .join(L.sign_symbol(v) for v in sv[:2]) + r"\mid " + r"\,".join(L.sign_symbol(v) for v in sv[2:]) + ")$"
        ax.text(pos[0], pos[1], label, ha="center", va="center", fontsize=11, color=style.INK if inpos else style.GREY,
                bbox=dict(boxstyle="round,pad=0.25", fc=style.HILITE_A if inpos else "white", ec="#bbbbbb", lw=0.6))
        if inpos:
            feasible.append(tuple(sv[2:]))
    ax.text(0.15, L_box - 0.45, "$x > 0$", fontsize=11, color="#8a6100")
    ax.set_xlabel(f"$x_{{{i0}}}$")
    ax.set_ylabel(f"$x_{{{i1}}}$")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r"cells of the arrangement carry $(\mathrm{sign}\,x \mid \mathrm{sign}\,Bx)$;"
                 "\n" + r"$\mathrm{sign}\,a$ is feasible iff $(+\,+ \mid \mathrm{sign}\,a)$ is a cell with $x>0$",
                 fontsize=10.5)
    feas = r",\ ".join(L.tex_signs(f) for f in feasible)
    offdiag = r",\ ".join(rf"b_{{{L._i(i)}{L._i(j)}}} {'>' if B[i, j] > 0 else '<'} 0"
                          for i in range(2) for j in range(2) if i != j)
    fig.text(0.5, 0.005, rf"${offdiag}$:   feasible $\mathrm{{sign}}\,a \in \{{{feas}\}}$",
             ha="center", fontsize=10.5, color=style.INK)
    save(fig, name)


# --------------------------------------------------------------------------- 3. Z-block twist (R4)

def fig_zblock(eco, name="zblock"):
    """
    R4 on the whole matrix (P = ∅): the first twist (s, σ) that ``Ecology.z_twists`` finds turns
    sign B into a Z-matrix with positive diagonal, shown as B → S B → S B P_σ.
    """
    n = eco.n
    Q = tuple(range(n))
    chi = eco.initial_table()
    _, D, A = eco.schur_block(chi, (), Q)
    twists = eco.z_twists(Q, D)
    if not twists:
        raise ValueError("no Z-twist of the whole matrix for this pattern")
    s_d, sig_d = twists[0]
    s = np.array([s_d[i] for i in Q])
    sigma = [sig_d[i] for i in Q]
    B = eco.B
    SB = np.diag(s) @ B
    SBP = SB[:, sigma]
    prod_s = int(np.prod(s))
    sgn_sigma = L.perm_sign(sigma) if hasattr(L, "perm_sign") else _perm_sign(sigma)
    det_sign = prod_s * sgn_sigma
    twisted_rhs = np.array([s_d[i] * A[i] for i in Q])

    fig, axes = plt.subplots(1, 3, figsize=(3.1 * n + 3.5, 1.0 * n + 1.6))
    sign_matrix(axes[0], B, r"$\mathrm{sign}\,B$")
    sign_matrix(axes[1], SB, rf"$S\,B$,   $S = \mathrm{{diag}}{L.tex_signs(s)}$")
    sign_matrix(axes[2], SBP, rf"$S\,B\,P_\sigma$,   $\sigma = {L.tex_perm(sigma)}$",
                col_labels=[str(L._i(j)) for j in sigma], outline_diag=True)
    for ax in axes[:2]:
        ax.annotate("", xy=(n + 0.55, n / 2), xytext=(n + 0.05, n / 2),
                    arrowprops=dict(arrowstyle="-|>", color=style.INK, lw=1.0), annotation_clip=False)
    axes[2].text(n / 2, -0.3, "Z-matrix, positive diagonal", ha="center", va="top", fontsize=10, color=style.INK)
    if np.all(twisted_rhs > 0):
        line1 = (r"$S\,B\,P_\sigma \cdot (P_\sigma^{-1} x^*) = S\,a > 0$ with $x^* > 0$  ⇒  M-matrix  ⇒  "
                 r"$\det(S\,B\,P_\sigma) > 0$")
    else:
        line1 = r"$S\,a < 0$  ⇒  $S\,B\,P_\sigma$ is not an M-matrix"
    line2 = (rf"$\det(S\,B\,P_\sigma) = \left(\prod_i s_i\right)\,\mathrm{{sgn}}(\sigma)\,\det B "
             rf"= ({'+' if prod_s > 0 else '-'}1)({'+' if sgn_sigma > 0 else '-'}1)\det B = "
             rf"{'' if det_sign > 0 else '-'}\det B$")
    fig.text(0.5, 0.04, line1 + "\n" + line2, ha="center", va="bottom", fontsize=11)
    plt.subplots_adjust(bottom=0.22, top=0.86, wspace=0.25)
    save(fig, name)


def _perm_sign(seq):
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


# --------------------------------------------------------------------------- 4. balanced (R3)

def fig_balanced(eco, name="balanced"):
    """R3: the signed interaction graph of a symmetric pattern, its balancing signature s, and sign(S B S)."""
    n = eco.n
    s = eco.balancing_signature()
    if s is None:
        raise ValueError("pattern is not symmetric and balanced")
    B = eco.B
    SBS = np.diag(s) @ B @ np.diag(s)
    fig, axes = plt.subplots(1, 3, figsize=(3.1 * n + 3.5, 1.0 * n + 1.6), gridspec_kw=dict(width_ratios=[1.15, 1, 1]))
    ax = axes[0]
    ang = np.pi / 2 - 2 * np.pi * np.arange(n) / n               # species 0 on top, clockwise
    pos = {i: (np.cos(ang[i]), np.sin(ang[i])) for i in range(n)}
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")
    for i in range(n):
        for j in range(i + 1, n):
            ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]], color=style.sign_line(B[i, j]), lw=2.2, zorder=1)
    for i in range(n):
        ax.add_patch(Circle(pos[i], 0.19, facecolor=style.sign_tint(s[i]), edgecolor=style.INK, lw=1.2, zorder=3))
        ax.text(pos[i][0], pos[i][1] + 0.04, f"${L._i(i)}$", ha="center", va="center", fontsize=11, zorder=4)
        ax.text(pos[i][0], pos[i][1] - 0.09, f"$s_{{{L._i(i)}}} = {L.sign_symbol(s[i])}$", ha="center", va="center",
                fontsize=7.5, color=style.sign_ink(s[i]), zorder=4)
    ax.set_title(r"edge colour $= \mathrm{sign}\,b_{ij}$,   node colour $= s_i$", fontsize=11, pad=14)
    sign_matrix(axes[1], B, r"$\mathrm{sign}\,B$")
    sign_matrix(axes[2], SBS, rf"$S\,B\,S$,   $S = \mathrm{{diag}}{L.tex_signs(s)}$", outline_diag=True)
    axes[2].text(n / 2, -0.3, "Z-matrix", ha="center", va="top", fontsize=10, color=style.INK)
    fig.text(0.5, 0.06,
             r"$s_i\, s_j\, b_{ij} < 0$ on every edge (balanced)  ⇒  $\mathrm{diag}(x^*)\,S\,B\,S$ is a Hurwitz Z-matrix, "
             r"hence an M-matrix" + "\n"
             r"⇒  every principal minor of $B$ is positive, every almost-principal minor has the sign of its $b_{ij}$",
             ha="center", va="bottom", fontsize=11)
    plt.subplots_adjust(bottom=0.22, top=0.86, wspace=0.25)
    save(fig, name)


# --------------------------------------------------------------------------- examples

def obligate_cycle(n):
    """i benefits only from i+1 (b_{i,i+1} < 0), competes with everyone else, a < 0."""
    B = np.ones((n, n), dtype=int)
    for i in range(n):
        B[i, (i + 1) % n] = -1
    return Ecology([-1] * n, B, f"obligate {n}-cycle")


def two_groups(groups):
    """Mutualism inside each group, competition between groups (symmetric, balanced)."""
    n = sum(len(g) for g in groups)
    B = np.ones((n, n), dtype=int)
    for g in groups:
        for i in g:
            for j in g:
                if i != j:
                    B[i, j] = -1
    return Ecology([1] * n, B, "two mutualistic groups")


if __name__ == "__main__":
    n4 = Ecology([1, -1, 1, -1], np.ones((4, 4), dtype=int), "")
    fig_dictionary(n4, S=(0, 2, 5, 7))
    fig_tope(B=[[1.0, -1.6], [0.7, 1.0]])            # species 2 helps 1 (b_12 < 0), 1 harms 2 (b_21 > 0)
    fig_zblock(obligate_cycle(4))
    fig_balanced(two_groups([[0, 1], [2, 3]]))
