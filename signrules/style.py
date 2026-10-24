"""
One palette for every drawing of the package.

Colour carries exactly three meanings (the interaction networks of ``network.py`` follow the
article's ``\\econet`` instead: ``FREE_BG`` … ``LOSS_LINE``):

* **sign** — red is positive, blue is negative — for entries of B (competition / mutualism in the
  interaction graphs), for Plücker coordinates in derivation graphs, and for the branch taken in a
  search tree;
* **status** of a node that has no sign: a branch choice (dashed border), a contradiction
  (black ⊥), a completion (green), an unexplored subtree (grey);
* **rule** — one hue per rule (``RULE_COLORS``): the headline of a rule box in a derivation
  graph, the ⊥ leaf that a rule closes in a search tree, the rule's box in the overview diagram.
  None of these hues is the red or the blue of the signs.

Everything else (matrices, arrows, hulls) is greyscale or a pale tint, so that the sign colours
are the only thing that jumps out.  Where a figure needs to point at a *subset* (chosen columns
of a matrix, a group of species, copies of a repeated subtree) it uses pale tints, never the
sign colours.
"""

import matplotlib

# sign colours: ink (text / lines) and tint (fills)
POS, NEG = "#b8321e", "#1f4fa3"
POS_TINT, NEG_TINT = "#f6d3cd", "#d2ddf5"
POS_LINE, NEG_LINE = "#d0503a", "#3a6bd0"

# interaction networks, as the article's \econet macro: node fill = growth sign, a benefit of i from j
# (b_ij < 0) a dashed arrow i -> j, a pair that loses both ways a thin line (compet!60), one-way losses not drawn
FREE_BG, OBLIGATE_BG = "#f9d4d1", "#cfdcf7"     # posbg / negbg: a_i > 0 free, a_i < 0 obligate
FREE_FG, OBLIGATE_FG = "#8b1a1a", "#1c3f9e"     # posfg / negfg: the sign text on those fills
BENEFIT = "#1f4fd6"                             # mutual
LOSS = "#d62728"                                # compet
LOSS_LINE = "#e67d7e"                           # compet!60 (60 % compet on white)
NODE_EDGE = "#666666"                           # black!60

# status colours
CHOICE = "#222222"                 # a branch choice: dashed border in this colour
BOTTOM = "#222222"                 # contradiction
COMPLETION = "#2e7d32"
STOPPED = "#b5b5b5"
RULE_INK, RULE_EDGE = "#333333", "#9a9a9a"
INK, GREY, LIGHT = "#222222", "#777777", "#dddddd"

# highlight tints for subsets (not signs)
HILITE_A, HILITE_B = "#fbe7bd", "#d9ead3"
QUADRANT = "#fff4d6"
# pale tints for hulls of repeated subtrees (categorical, in order of use)
HULLS = ["#f3e6c8", "#dcefd8", "#e6dff2", "#d9eef3", "#f6e2ea", "#ececec"]

# one hue per rule; GP = the 3-term Grassmann–Plücker relations
RULE_COLORS = {"GP": "#6f6f6f", "R1": "#c98a00", "R2": "#3aa89a", "R3": "#8c6d1f", "R4": "#7b3fa0",
               "R5": "#d9662b", "R6": "#b8336a", "R7": "#006d77",
               "R8": "#4a6fa5"}


def rule_color(rule: str) -> str:
    """The hue of a rule (``RULE_COLORS``), ink colour for anything else."""
    return RULE_COLORS.get(rule, INK)


def sign_ink(sign) -> str:
    """Text and line colour of a sign: red for +, blue for −."""
    return POS if sign > 0 else NEG


def sign_tint(sign) -> str:
    """Fill colour of a sign: pale red for +, pale blue for −."""
    return POS_TINT if sign > 0 else NEG_TINT


def sign_line(sign) -> str:
    """Edge colour of a sign in a search tree: red for +, blue for −."""
    return POS_LINE if sign > 0 else NEG_LINE


def use_math_fonts() -> None:
    """Computer Modern mathtext (LaTeX look) with a plain sans text font; call once before drawing."""
    matplotlib.rcParams["mathtext.fontset"] = "cm"
    matplotlib.rcParams["mathtext.rm"] = "serif"
    matplotlib.rcParams["font.family"] = "DejaVu Sans"
    matplotlib.rcParams["axes.unicode_minus"] = True
