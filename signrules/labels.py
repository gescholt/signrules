"""
Text for figures and reports, generated from the data — never typed by hand.

Everything the drawings say is produced here, so that the certificate DAGs, search trees,
completion derivations and the rule illustrations all speak the same language:

* a sign is a **Plücker coordinate** ``p_S`` of the point ``[diag(a) | B]`` of ``Gr(n, 2n)``,
  written ``$p_{0257} > 0$``; the minor of ``B`` behind it (``det B_{P,Q}``) is a gloss that is
  only shown when a rule is stated in terms of it (``style="minor"`` / ``"both"``);
* a rule application is a short headline (``R4 · Z-block``) and a line of mathematics giving its
  parameters (``$Q=\\{1,2,3\\},\\ s=(-,-,-),\\ \\sigma=(2\\,3\\,1)$``), parsed from the free-text
  ``why`` strings that ``rules.Ecology`` logs;
* a contradiction is ``⊥`` with the rule that raised it.

Three output formats: ``"tex"`` (matplotlib mathtext, the default), ``"html"`` (Graphviz /
plotly labels) and ``"text"`` (terminal).  Indices are printed as the code stores them (0-based);
set ``ONE_BASED = True`` to print them as in the papers.
"""

import re
from typing import List, Optional, Sequence, Tuple

from .plucker_index import UNKNOWN

ONE_BASED = False          # print species / column indices from 1 instead of 0

MINUS = "−"                # U+2212, for html / text output


# --------------------------------------------------------------------------- small pieces

def _i(k: int) -> str:
    return str(k + 1 if ONE_BASED else k)


def idx(seq: Sequence[int]) -> str:
    """Indices run together: (0, 2, 5, 7) -> '0257'."""
    return "".join(_i(k) for k in seq)


def tex_set(seq: Sequence[int]) -> str:
    return r"\emptyset" if not seq else r"\{" + ",".join(_i(k) for k in seq) + r"\}"


def tex_signs(signs: Sequence[int]) -> str:
    """(+1, -1, -1) -> '(+,-,-)'."""
    return "(" + ",".join("+" if s > 0 else "-" for s in signs) + ")"


def tex_perm(seq: Sequence[int]) -> str:
    """A permutation written as the sequence of images, e.g. (2\\,3\\,1)."""
    return "(" + r"\," .join(_i(k) for k in seq) + ")"


def sign_symbol(sign, fmt: str = "tex") -> str:
    if sign == UNKNOWN or sign is None:
        return "?"
    if fmt == "tex":
        return "+" if sign > 0 else "-"
    return "+" if sign > 0 else MINUS


# --------------------------------------------------------------------------- format conversion

_TEX_TO_HTML = [
    (r"\det", "det"), (r"\emptyset", "∅"), (r"\sigma", "σ"), (r"\lambda", "λ"), (r"\bot", "⊥"),
    (r"\{", "{"), (r"\}", "}"), (r"\,", " "), (r"\ ", " "), (r"\cdot", "·"), (r"\leq", "≤"), (r"\le", "≤"),
    (r"\Rightarrow", "⇒"), (r"\times", "×"), (r"\ne", "≠"), (r"\varepsilon", "ε"), (r"\sum", "Σ"),
    (r"\mathrm", ""), (r"\operatorname", ""),
]


def _strip_math(s: str, sub_open: str, sub_close: str, sup_open: str, sup_close: str) -> str:
    """Turn a mathtext string into plain html/text: braces, subscripts, superscripts, symbols."""
    out = []
    parts = s.split("$")
    for k, part in enumerate(parts):
        if k % 2 == 0:                      # outside math
            out.append(part)
            continue
        t = part.replace(">", " > ").replace("<", " < ").replace("=", " = ")
        t = re.sub(r"\\(mathrm|operatorname)\{([^{}]*)\}", r"\2", t)
        t = t.replace(r"\{", "\u27e8").replace(r"\}", "\u27e9")        # set braces, restored below
        for a, b in _TEX_TO_HTML:
            t = t.replace(a, b)
        # subscripts / superscripts with braces, then single characters
        t = re.sub(r"_\{([^{}]*)\}", lambda m: sub_open + m.group(1) + sub_close, t)
        t = re.sub(r"\^\{([^{}]*)\}", lambda m: sup_open + m.group(1) + sup_close, t)
        t = re.sub(r"_(\w)", lambda m: sub_open + m.group(1) + sub_close, t)
        t = re.sub(r"\^(\w|\*)", lambda m: sup_open + m.group(1) + sup_close, t)
        t = t.replace("\u27e8", "{").replace("\u27e9", "}")
        t = t.replace(" - ", " " + MINUS + " ").replace("(-", "(" + MINUS).replace(",-", "," + MINUS)
        t = re.sub(r"(?<![\w)])-(?=[\w(])", MINUS, t)      # unary minus
        t = re.sub(r"\s+", " ", t).strip()
        out.append(t)
    return "".join(out)


def convert(tex: str, fmt: str) -> str:
    """A tex label in the requested format."""
    if fmt == "tex":
        return tex
    if fmt == "html":
        s = _strip_math(tex, "<sub>", "</sub>", "<sup>", "</sup>")
        return s.replace("&", "&amp;").replace(" < ", " &lt; ").replace(" > ", " &gt; ")
    return _strip_math(tex, "_", "", "^", "")


# --------------------------------------------------------------------------- bases

def plucker(S: Sequence[int], sign=UNKNOWN, fmt: str = "tex") -> str:
    """The Plücker coordinate p_S, with its sign when known: '$p_{0257} > 0$'."""
    core = f"p_{{{idx(S)}}}"
    if sign == UNKNOWN or sign is None:
        return convert(f"${core}$", fmt)
    return convert(f"${core} {'>' if sign > 0 else '<'} 0$", fmt)


def minor_name(eco, S: Sequence[int]) -> str:
    """The minor of B behind basis S, as tex without $: '\\det B_{\\{1,3\\},\\{0,2\\}}'."""
    _, P, Q, _ = eco.idx.decompose(S)
    if not P:
        return "1"
    if len(P) == eco.n:
        return r"\det B"
    if P == Q:
        return r"\det B_{" + tex_set(P) + "}"
    return r"\det B_{" + tex_set(P) + "," + tex_set(Q) + "}"


def minor(eco, S: Sequence[int], sign=UNKNOWN, fmt: str = "tex") -> str:
    """The minor of B behind basis S with the sign that the Plücker sign ``sign`` gives it."""
    name = minor_name(eco, S)
    if sign == UNKNOWN or sign is None:
        return convert(f"${name}$", fmt)
    _, P, Q, _ = eco.idx.decompose(S)
    _, eps = eco.idx.key(P, Q)
    m = eps * eco.idx.a_factor(eco.a, P) * sign
    return convert(f"${name} {'>' if m > 0 else '<'} 0$", fmt)


def dictionary(eco, S: Sequence[int], fmt: str = "tex") -> str:
    """The dictionary entry for S:  p_S = ε a_T · det B_{P,Q}  with the actual sign ε."""
    T, P, Q, eps = eco.idx.decompose(S)
    a_part = "".join(f"a_{{{_i(j)}}}" for j in T)
    sign = "" if eps > 0 else "-"
    parts = [sign + a_part, minor_name(eco, S)] if a_part else [sign + minor_name(eco, S)]
    rhs = r"\," .join(p for p in parts if p not in ("", "1")) or "1"
    return convert(f"$p_{{{idx(S)}}} = {rhs}$", fmt)


def basis(eco, S: Sequence[int], sign=UNKNOWN, style: str = "plucker", fmt: str = "tex") -> List[str]:
    """Label lines for a basis: ['$p_S > 0$'] or with the minor gloss as a second line."""
    if style == "minor":
        return [minor(eco, S, sign, fmt)]
    lines = [plucker(S, sign, fmt)]
    if style == "both":
        lines.append(minor(eco, S, sign, fmt))
    return lines


# --------------------------------------------------------------------------- rules

_TUPLE = r"\(([^()]*)\)"


def _ints(s: str) -> Tuple[int, ...]:
    return tuple(int(x) for x in re.findall(r"-?\d+", s))


class RuleLabel:
    """A parsed rule application: ``rule`` (R1..R8 or GP), ``name``, ``detail`` (tex) and ``clause``."""

    def __init__(self, rule: str, name: str, detail: str = "", clause: str = ""):
        self.rule, self.name, self.detail, self.clause = rule, name, detail, clause

    def headline(self, fmt: str = "tex") -> str:
        return f"{self.rule} · {self.name}" if self.name else self.rule

    def lines(self, fmt: str = "tex", with_clause: bool = False) -> List[str]:
        out = [self.headline(fmt)]
        if self.detail:
            out.append(convert(self.detail, fmt))
        if with_clause and self.clause:
            out.append(convert(self.clause, fmt))
        return out

    def key(self) -> Tuple[str, str, str]:
        """What identifies one application (the clause is dropped: one box, several outputs)."""
        return (self.rule, self.name, self.detail)

    def __repr__(self):
        return f"RuleLabel({self.rule}, {self.name!r}, {self.detail!r}, {self.clause!r})"


def rule(why: str, n: Optional[int] = None) -> RuleLabel:
    """Parse a ``why`` string logged by ``rules.Ecology`` (or the GP marker) into a RuleLabel."""
    if why.startswith("Grassmann") or "Grassmann" in why.split(":")[0]:
        return RuleLabel("GP", "3-term relation")
    head, _, clause = why.partition(":")
    head = head.strip()
    clause = clause.strip()

    if head.startswith("R1"):
        m = re.search(r"λ\^\{?(\d+)\}?.*principal (\d+)-minors", head)
        if m:
            k = m.group(2)
            return RuleLabel("R1", "stability", rf"$c_{{{k}}} = \sum_{{|I|={k}}} \det B_{{I}} > 0$")
        return RuleLabel("R1", "stability", r"$\det B > 0$")

    if head.startswith("R2"):
        m = re.search(r"x\*_(\d+).*J=" + _TUPLE, why)
        if m:
            j, J = int(m.group(1)), _ints(m.group(2))
            return RuleLabel("R2", "Cramer", rf"$x^*_{{{_i(j)}}}\,\det B_{{{tex_set(J)}}}$, $J={tex_set(J)}$")
        return RuleLabel("R2", "Cramer")

    if head.startswith("R3"):
        return RuleLabel("R3", "balanced", "", _clause_tex(clause))

    if head.startswith("R4"):
        m = re.search(r"P=" + _TUPLE + r" Q=" + _TUPLE + r" s=" + _TUPLE + r" σ=" + _TUPLE, head)
        P, Q, s, sigma = (_ints(g) for g in m.groups())
        detail = (rf"$P={tex_set(P)},\ Q={tex_set(Q)}$" if P else rf"$Q={tex_set(Q)}$") + "\n" + \
                 rf"$s={tex_signs(s)},\ \sigma={tex_perm(sigma)}$"
        return RuleLabel("R4", "Z-block", detail, _clause_tex(clause))

    if head.startswith("R5"):
        m = re.search(r"det B_" + _TUPLE + "," + _TUPLE, head)
        P, Q = _ints(m.group(1)), _ints(m.group(2))
        if n is not None and len(P) == n:
            name = r"\det B"
        else:
            name = r"\det B_{" + tex_set(P) + ("" if P == Q else "," + tex_set(Q)) + "}"
        return RuleLabel("R5", "Laplace", f"${name}$")

    if head.startswith("R6"):
        m = re.search(r"P=" + _TUPLE + r" Q=" + _TUPLE, head)
        P, Q = _ints(m.group(1)), _ints(m.group(2))
        return RuleLabel("R6", "Schur block", rf"$P={tex_set(P)},\ Q={tex_set(Q)}$", _clause_tex(clause))

    if head.startswith("R7"):
        if "one-sided" in clause:
            return RuleLabel("R7", "tope", "", r"circuit one-sided against $(+,\mathrm{sign}\,a)$")
        return RuleLabel("R7", "tope", "", clause)
    if head.startswith("R8"):
        m = re.search(r"s=" + _TUPLE + r" σ=" + _TUPLE, head)
        s, sigma = (_ints(g) for g in m.groups())
        return RuleLabel("R8", "real eigenvalue", rf"$s={tex_signs(s)},\ \sigma={tex_perm(sigma)}$",
                         r"$S B P_\sigma$ M-matrix, $\sigma\neq\mathrm{id}$: real eigenvalue $>0$")
    if head.startswith("SYM"):
        m = re.search(r"π=" + _TUPLE, head)
        return RuleLabel("SYM", "symmetry", rf"$\pi={tex_perm(_ints(m.group(1)))}$", _clause_tex(clause))
    return RuleLabel(head.split()[0] if head else "?", head)


def _clause_tex(clause: str) -> str:
    """The sub-clause after the colon, lightly mathematised."""
    c = clause
    c = re.sub(r"ordering " + _TUPLE, lambda m: "ordering $" + tex_perm(_ints(m.group(1))) + "$", c)
    c = c.replace(" > 0", " $>0$").replace(" < 0", " $<0$").replace(" ≤ 0", r" $\leq 0$")
    c = re.sub(r"row (\d+) of D x = A", lambda m: rf"row ${_i(int(m.group(1)))}$ of $D\,x^*_Q = A_Q$", c)
    c = c.replace("not an M-matrix, but ", "no M-matrix: ").replace("not an M-matrix, ", "no M-matrix: ")
    c = c.replace("M-matrix, ", "")
    return c


def conflict(eco, reason: str, fmt: str = "tex") -> List[str]:
    """Label lines of a contradiction node: ⊥, the rule, and what clashed (in Plücker language)."""
    lab = rule(reason, None if eco is None else eco.n)
    lines = [r"$\bot$" if fmt == "tex" else "⊥", lab.headline(fmt)]
    if lab.detail:
        lines.append(convert(lab.detail, fmt))
    tail = reason.rpartition(":")[2].strip() if ":" in reason else ""
    m = re.search(r"det B_" + _TUPLE + "," + _TUPLE + r" must be ([<>]) 0 but is ([<>]) 0", tail)
    if m and eco is not None:
        P, Q = _ints(m.group(1)), _ints(m.group(2))
        S, eps = eco.idx.key(P, Q)
        want = eps * eco.idx.a_factor(eco.a, P) * (1 if m.group(3) == ">" else -1)
        lines.append(convert(rf"$p_{{{idx(S)}}}$ forced ${sign_symbol(want)}$, already ${sign_symbol(-want)}$", fmt))
    elif "every term is" in tail:
        m = re.search(r"every term is ([<>]) 0 but the sum must be ([<>]) 0", tail)
        lines.append(convert(rf"every term ${m.group(1)}0$, sum must be ${m.group(2)}0$", fmt))
    elif tail and tail != reason:
        lines.append(convert(_clause_tex(tail), fmt))
    return lines


# --------------------------------------------------------------------------- patterns

def pattern(eco, fmt: str = "tex") -> str:
    """'$n=4$, $a=(+,-,-,-)$' — the data of a run, without the matrix."""
    return convert(rf"$n={eco.n}$, $\mathrm{{sign}}\,a={tex_signs(eco.a)}$", fmt)

