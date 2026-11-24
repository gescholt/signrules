"""
Live terminal dashboard of one sign-completion search.

    signrules run --a=-1,-1,-1,-1,1 --B=1,-1,1,1,1/-1,1,1,1,1/1,1,1,-1,1/1,1,-1,1,1/-1,-1,-1,-1,1
    signrules run --preset parity_k2 --rules R1,R2n --html tree.html
    signrules run --random 5 --max-nodes 5000 --html tree.html   # a random pattern (--seed S to repeat one)
    signrules run                    # menus

The menus need questionary (the ``live`` extra: ``uv tool install 'ecological_sign_rules[live]'``, or
``uv sync --extra live`` in a checkout).  Write ``--a=-1,...`` with ``=``: a value starting with ``-`` is otherwise read
as an option (or write the signs as ``+--``).  ``signrules.cli`` holds the command line, ``signrules.report``
the verdict printed at the end.

    from signrules.livesearch import LiveSearch
    with LiveSearch(eco, use) as live:
        result = complete(eco, use=use, progress=live)

What it shows, refreshed several times a second:

* **progress** — how far the depth-first sweep has got: every finished leaf at depth d closes 2^-d
  of the binary search space. The estimate of the time left extrapolates that fraction; it is
  reliable for bushy trees and pessimistic early in deep, thin ones (a long first branch sweeps
  almost nothing until it closes). Subtrees cut off by a node limit are shown grey;
* **counters** — nodes, nodes per second, completions, contradictions, current and maximal depth,
  open signs at the current node;
* **rules** — for each rule (GP, R1–R8, in the package's rule colours) how many signs it forced
  and how many branches it closed with a contradiction;
* **current path** — the branch choices from the root to the node being searched (red +, blue −);
* **leaves by depth** — where the tree ends: completions and contradictions per depth;
* **events** — the latest completions, contradictions (with their rule) and new depth records.
"""

import time
from collections import Counter, deque
from typing import Optional, Sequence

from rich.console import Group
from rich.live import Live
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from . import labels as L
from .patterns import TREE_NOTES  # re-exported: tests and scripts import it from here
from . import style
from .rules import Ecology

RULES = ["GP", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
POS, NEG = style.POS_LINE, style.NEG_LINE
BLOCKS = " ▁▂▃▄▅▆▇█"


def _rule_of(why: str) -> str:
    return why[:2] if why[:1] == "R" and why[1:2].isdigit() else ("GP" if "Grassmann" in why else why[:2])


def _milestone(k: int) -> bool:
    """1, 2, 5, 10, 20, 50, 100, ...: which occurrences of a repeated event get logged."""
    while k >= 10 and k % 10 == 0:
        k //= 10
    return k in (1, 2, 5)


def _coord(S) -> str:
    return "".join(map(str, S))


class LiveSearch:
    """Progress callback for ``complete`` that keeps statistics and renders them with rich."""

    def __init__(self, eco: Ecology, use: Optional[Sequence[str]] = None, events: int = 8,
                 refresh_per_second: float = 8.0):
        self.eco, self.use = eco, list(use) if use else list(Ecology.ALL_RULES)
        self.t0 = time.time()
        self.nodes = self.completions = self.conflicts = self.stopped = 0
        self.depth = self.max_depth = 0
        self.unknown = self.root_unknown = None
        self.path = []
        self.covered = 0.0                     # share of the search space closed (completions, ⊥)
        self.abandoned = 0.0                   # share left unsearched because a limit stopped the search
        self.forced = Counter()
        self.closed = Counter()
        self.leaves = {}                       # depth -> Counter(kind)
        self.events = deque(maxlen=events)
        self.finished = False
        self._live = Live(self, refresh_per_second=refresh_per_second, transient=False)

    # ------------------------------------------------------------------ callback
    def __call__(self, ev):
        kind, d = ev["kind"], ev["depth"]
        for _ in ev.get("gp_log") or ():
            self.forced["GP"] += 1
        for entry in ev.get("log") or ():
            self.forced[_rule_of(entry[0])] += 1
        if kind == "root":
            self.root_unknown = self.unknown = ev["unknown"]
            self._event("root", f"root propagation leaves {ev['unknown']} signs open")
            return
        if kind == "forced":
            return
        self.path = ev["path"]
        self.depth = d
        if d > self.max_depth:
            self.max_depth = d
            if kind == "node":
                self._event("depth", f"new depth record {d}")
        if kind == "node":
            self.nodes += 1
            self.unknown = ev["unknown"]
            return
        # a leaf
        if kind == "stopped":
            self.abandoned += 2.0 ** -d
        else:
            self.covered += 2.0 ** -d
        self.leaves.setdefault(d, Counter())[kind] += 1
        if kind == "completion":
            self.nodes += 1
            self.completions += 1
            if _milestone(self.completions):
                self._event("completion", f"completion {self.completions:,} at depth {d}")
        elif kind == "conflict":
            self.conflicts += 1
            lab = L.rule(ev.get("reason", ""), self.eco.n)
            self.closed[lab.rule] += 1
            if _milestone(self.closed[lab.rule]):
                k = self.closed[lab.rule]
                suffix = "th" if 10 <= k % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")
                nth = "first" if k == 1 else f"{k:,}{suffix}"
                self._event("conflict", f"{nth} ⊥ by {lab.headline('plain')} at depth {d}")
        elif kind == "stopped":
            self.stopped += 1
            if self.stopped == 1:
                self._event("stopped", f"limit reached at depth {d}: the rest of the tree is not searched")

    def _event(self, kind, text):
        self.events.appendleft((time.time() - self.t0, kind, text))

    # ------------------------------------------------------------------ rendering
    def _header(self):
        t = Text()
        t.append(self.eco.name or "pattern", style="bold")
        t.append("   sign a ")
        for v in self.eco.a:
            t.append("+ " if v > 0 else "− ", style=style.POS if v > 0 else style.NEG)
        t.append("  rules " + ", ".join(self.use), style="dim")
        B = Text()
        for row in self.eco.B:
            for v in row:
                B.append(" + " if v > 0 else " − ", style=f"bold white on {style.POS if v > 0 else style.NEG}")
            B.append("\n")
        grid = Table.grid(padding=(0, 2))
        grid.add_row(t, B)
        return grid

    def _progress(self):
        el = time.time() - self.t0
        f = min(self.covered, 1.0)
        g = min(self.abandoned, 1.0 - f)
        width = 40
        full, skip = int(f * width), int(g * width)
        bar = (Text("█" * full, style=style.COMPLETION) + Text("▒" * skip, style=style.STOPPED)
               + Text("░" * (width - full - skip), style="dim"))
        if self.finished:
            eta = "stopped by a limit" if self.stopped else "done"
        elif f > 0:
            eta = f"~{el * (1 - f) / f:,.0f}s left"
        else:
            eta = "estimating…"
        tail = f" · {100 * g:.2f}% not searched" if g else ""
        return Group(Text.assemble(bar, f"  {100 * f:6.2f}% swept{tail}   ",
                                   (f"{el:,.1f}s elapsed · {eta}", "dim")))

    def _counters(self):
        el = max(time.time() - self.t0, 1e-9)
        tb = Table.grid(padding=(0, 1))
        tb.add_column(justify="right", style="dim"); tb.add_column(justify="left")
        tb.add_row("nodes", f"{self.nodes:,}  ({self.nodes / el:,.0f}/s)")
        tb.add_row("completions", Text(f"{self.completions:,}", style=style.COMPLETION))
        tb.add_row("contradictions", f"{self.conflicts:,}")
        if self.stopped:
            tb.add_row("not searched", Text(f"{self.stopped:,}", style="dim"))
        tb.add_row("depth", f"{self.depth}  (max {self.max_depth})")
        if self.root_unknown is not None:
            tb.add_row("open signs", f"{self.unknown}  (root {self.root_unknown})")
        return Panel(tb, title="counters", border_style="dim")

    def _rules(self):
        tb = Table(box=None, padding=(0, 1), show_header=True, header_style="dim")
        tb.add_column("rule"); tb.add_column("signs forced", justify="right"); tb.add_column("")
        tb.add_column("⊥ closed", justify="right")
        top = max(self.forced.values(), default=1) or 1
        for r in RULES:
            if not (self.forced[r] or self.closed[r]):
                continue
            col = style.rule_color(r)
            bar = Text("█" * max(1, round(16 * self.forced[r] / top)) if self.forced[r] else "", style=col)
            tb.add_row(Text(r, style=f"bold {col}"), f"{self.forced[r]:,}", bar,
                       f"{self.closed[r]:,}" if self.closed[r] else "")
        return Panel(tb, title="rules", border_style="dim")

    def _path(self):
        t = Text()
        shown = self.path[-14:]
        if len(self.path) > len(shown):
            t.append(f"… {len(self.path) - len(shown)} more  ", style="dim")
        for S, sign in shown:
            t.append(f"p{_coord(S)}{'+' if sign > 0 else '−'} ", style=POS if sign > 0 else NEG)
        if not self.path:
            t.append("the root", style="dim")
        return Panel(t, title=f"current path · depth {len(self.path)}", border_style="dim")

    def _leaves(self):
        if not self.leaves:
            return Panel(Text("no leaves yet", style="dim"), title="leaves by depth", border_style="dim")
        depths = range(0, max(self.leaves) + 1)
        top = max(sum(c.values()) for c in self.leaves.values())
        comp = Text(); conf = Text(); axis = Text()
        for d in depths:
            c = self.leaves.get(d, Counter())
            comp.append(BLOCKS[round(8 * c["completion"] / top)] if c["completion"] else " ", style=style.COMPLETION)
            conf.append(BLOCKS[round(8 * c["conflict"] / top)] if c["conflict"] else " ", style="bold")
            axis.append(str(d % 10), style="dim")
        g = Table.grid(padding=(0, 1))
        g.add_row(Text("completions", style="dim"), comp)
        g.add_row(Text("contradictions", style="dim"), conf)
        g.add_row(Text("depth", style="dim"), axis)
        return Panel(g, title="leaves by depth", border_style="dim")

    def _events(self):
        tb = Table.grid(padding=(0, 1))
        colors = {"completion": style.COMPLETION, "conflict": "bold", "depth": "yellow", "root": "dim",
                  "stopped": style.STOPPED}
        for t, kind, text in self.events:
            tb.add_row(Text(f"{t:7.2f}s", style="dim"), Text(text, style=colors.get(kind, "")))
        return Panel(tb, title="events", border_style="dim")

    GLOSS = ("The search fixes one unknown minor sign per step (a node); a leaf where every sign is fixed is a "
             "completion, one where a rule is violated is a contradiction. Rules on the right count the signs "
             "they forced and the branches they closed.")

    def _gloss(self) -> Text:
        """One orienting sentence under the header: what a node, a completion and a contradiction are."""
        return Text(self.GLOSS, style="dim", overflow="fold")

    def __rich__(self):
        mid = Table.grid(expand=True, padding=(0, 1))
        mid.add_column(ratio=1); mid.add_column(ratio=1)
        mid.add_row(self._counters(), self._rules())
        return Group(self._header(), self._gloss(), self._progress(), mid, self._path(), self._leaves(),
                     self._events())

    # ------------------------------------------------------------------ context manager
    def __enter__(self):
        self.t0 = time.time()
        self._live.__enter__()
        return self

    def __exit__(self, *exc):
        self.finished = True
        self._live.refresh()
        return self._live.__exit__(*exc)


# --------------------------------------------------------------------------- command line

def _presets():
    from .patterns import presets
    return presets()


# ------------------------------------------------------------------ the interactive picker (questionary)

RULE_NOTES = {
    "R1": "stability: det B > 0, some principal k-minor positive for every k",
    "R2": "Cramer on every subsystem: x*_j det B_JJ as a signed sum of minors",
    "R2n": "Cramer on the full system only (with R1: the algorithm of Çelik–Haas–Scholten–Wang–Zucal, Ann. Comb. 2026)",
    "R3": "balanced patterns: S B S a Z-matrix forces its principal minors",
    "R4": "Z-blocks: M-matrix blocks under a signature and a permutation (subsumed by R2,R3,R5,R7,R8: idle by default)",
    "R5": "Laplace expansions of every minor",
    "R6": "Schur blocks: minors of the Schur complement (subsumed by GP and R5,R7: idle by default)",
    "R7": "tope: circuit orthogonality (subsumes R2)",
    "R8": "real eigenvalue of benefactor-permutation patterns",
}
STOP_CHOICES = [("after 5 completions", (5, None)), ("after 20 completions", (20, None)),
                ("after 200 completions", (200, None)), ("after 2,000 nodes", (None, 2000)),
                ("after 20,000 nodes", (None, 20000)), ("after 200,000 nodes", (None, 200000)),
                ("run to the end (can be very long)", (None, None))]


def _signs(a) -> str:
    return " ".join("+" if v > 0 else "−" for v in a)


def _group(name: str) -> str:
    """The menu group of a preset."""
    if name.startswith("tree_"):
        return "search trees worth watching (each with its own rules)"
    if name.startswith("r8_") or name == "parity_k2":
        return "killed by R8 alone at the root, open under R1–R7"
    if name.startswith("n6_hard_"):
        return "six species, possible: huge trees (certificate of possibility in signrules.data.possible(6))"
    if name[:2] in ("S-", "A-"):
        return "four species, proved impossible: killed at the root (R1,R2n gives a tree)"
    return "five species, open patterns near the kills: huge trees"


def questionary_kit():
    """(questionary, Choice, Separator, menu kwargs, prompt kwargs), or None when questionary is missing."""
    try:
        import questionary
        from questionary import Choice, Separator, Style
    except ImportError:                                   # pragma: no cover
        print("the interactive menus need questionary:  uv tool install 'ecological_sign_rules[live]'  (in a checkout: uv sync --extra live)")
        return None
    st = Style([("qmark", "fg:#5f87ff bold"), ("question", "bold"), ("pointer", "fg:#5f87ff bold"),
                ("highlighted", "fg:#5f87ff bold"), ("selected", "fg:#5faf5f"), ("separator", "fg:#8a8a8a"),
                ("instruction", "fg:#8a8a8a italic"), ("answer", "fg:#5faf5f bold")])
    kw = dict(style=st, qmark="◆", pointer="▶", use_shortcuts=False, use_indicator=False,
              instruction="(↑↓ to move, enter to choose)")   # menus
    kt = dict(style=st, qmark="◆")                            # text and yes/no prompts
    return questionary, Choice, Separator, kw, kt


def _pick_pattern(presets: dict):
    """Arrow-key menu for the pattern: {"preset", "a", "B", "random"}, or None when the user quits."""
    kit = questionary_kit()
    if kit is None:
        return None
    questionary, Choice, Separator, kw, kt = kit
    groups: dict = {}
    for name in presets:
        groups.setdefault(_group(name), []).append(name)
    width = max(len(n) for n in presets) if presets else 10
    choices = []
    for g, names in groups.items():
        choices.append(Separator(f"  {g}"))
        for n in names:
            a, _ = presets[n]
            note = TREE_NOTES.get(n, "")
            choices.append(Choice(f"{n:<{width}}   n = {len(a)}   a = {_signs(a)}" + (f"   {note}" if note else ""), value=n))
    choices += [Separator("  other"), Choice(f"{'custom':<{width}}   edit the sign matrix (or type a and B)", value="custom"),
                Choice(f"{'random':<{width}}   a random pattern", value="random")]
    pattern = questionary.select("Pattern", choices=choices, **kw).ask()
    if pattern is None:
        return None
    out = {"preset": None, "a": None, "B": None, "random": None}
    if pattern == "random":
        n = questionary.select("Species", choices=[str(k) for k in range(3, 8)], default="5", **kw).ask()
        if n is None:
            return None
        out["random"] = int(n)
    elif pattern == "custom":
        from .cli import edited_pattern
        picked = edited_pattern()                         # the matrix editor, when textual is installed
        if picked is not None:
            return picked
        out["a"] = questionary.text("Growth-rate signs a (1,-1,... or +-...)", default="-1,-1,-1", **kt).ask()
        out["B"] = questionary.text("Signs of B (rows separated by /)", default="1,1,1/-1,1,1/-1,-1,1", **kt).ask()
        if out["a"] is None or out["B"] is None:
            return None
    else:
        out["preset"] = pattern
    return out


def _pick(presets: dict):
    """Arrow-key menus for the pattern, the rules and the stopping rule; None when the user quits."""
    out = _pick_pattern(presets)
    if out is None:
        return None
    questionary, Choice, Separator, kw, kt = questionary_kit()
    from .patterns import preset_rules
    default = preset_rules(out["preset"] or "") or Ecology.ALL_RULES
    rule_choices = [Choice(f"{r:<4} {RULE_NOTES[r]}", value=r, checked=r in default)
                    for r in ("R1", "R2", "R2n", "R3", "R4", "R5", "R6", "R7", "R8")]
    menu = {k: v for k, v in kw.items() if k not in ("use_shortcuts", "use_indicator", "instruction")}
    rules = questionary.checkbox("Rules", choices=rule_choices, instruction="(space to toggle, enter to confirm)",
                                 **menu).ask()
    if rules is None:
        return None
    out["rules"] = tuple(rules)
    stop = questionary.select("Stop", choices=[Choice(t, value=v) for t, v in STOP_CHOICES], **kw).ask()
    if stop is None:
        return None
    out["max_completions"], out["max_nodes"] = stop
    html = questionary.confirm("Write the interactive tree viewer (HTML) at the end?", default=False, **kt).ask()
    if html is None:
        return None
    out["html"] = questionary.text("File", default="tree.html", **kt).ask() if html else None
    return out

