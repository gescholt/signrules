"""
The full-screen app of ``signrules`` (also ``signrules edit``): every step without flags.

    signrules                         # opens with the pattern library
    signrules edit --preset A-e       # opens with a pattern

The left panel is the n × 2n sign matrix [diag(a) | B]: the diagonal block holds the signs of a (only its
diagonal is editable), the right block holds B (its diagonal is fixed +).  Column numbers are those of
[diag(a) | B], so the Plücker coordinate p_0257 is the minor on the columns 0, 2, 5, 7.  The right panel
draws the interaction network as the article's ``\\econet`` (``network.py``).  Below them, recomputed after
every change: whether the pattern is trivial, which library pattern it is up to relabelling, the verdict
of the root propagation, and the rules and node limit of the search.

Keys: the arrows move over the editable signs, tab jumps between a and B, space or enter flips, + and -
set a sign (a click flips too); l opens the library, < and > (or 2–7, or n for a chooser) change the
number of species, u chooses the rules and the node limit, r runs the search on a screen with the live
dashboard and then the verdict (where h writes the tree viewer and c the certificate drawing), e exports
the network (png, pdf and the \\econet call), s shows the classification, ? shows the help, q quits.

Needs textual (the ``live`` extra).
"""
from typing import Dict, List, Optional, Sequence, Tuple

from rich.text import Text
from textual import events, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import ModalScreen, Screen
from textual.widget import Widget
from textual.widgets import Footer, OptionList, Static
from textual.widgets.option_list import Option

from . import style
from .network import network_text

CELL = 3                                                  # characters per matrix cell
ROW_LABEL = 4                                             # characters before the first cell


# --------------------------------------------------------------------------- the pattern library

def library() -> List[Tuple[str, str, List[int], List[List[int]]]]:
    """(group, id, a, B) for the presets and the minimal impossible five- and six-species patterns."""
    from .data import root_kills
    from .livesearch import _group
    from .patterns import presets, short_id
    out = [(_group(k), short_id(k) if k[:2] in ("S-", "A-") else k, a, B) for k, (a, B) in presets().items()]
    for n in (5, 6):
        for k, r in enumerate(root_kills(n), 1):
            rule = (r.get("info") or "").split()[0] if r.get("info") else ""
            out.append((f"minimal impossible, {n} species (the extension screen)",
                        f"n{n}-kill-{k:02d}" + (f"  {rule}" if rule else ""), list(r["a"]), [list(x) for x in r["B"]]))
    return out


_INDEX: Optional[Dict[tuple, List[str]]] = None


def library_matches(a, B) -> List[str]:
    """The library ids of the patterns equal to (a, B) up to relabelling."""
    global _INDEX
    from .patterns import key
    if _INDEX is None:
        _INDEX = {}
        for _, ident, la, lB in library():
            _INDEX.setdefault((len(la), key(la, lB)), []).append(ident.split()[0])
    return _INDEX.get((len(a), key(a, B)), [])


# --------------------------------------------------------------------------- the matrix

def editable(n: int, row: int, col: int) -> bool:
    """a_i at (i, i); b_ij at (i, n + j) for j ≠ i."""
    return col == row if col < n else col - n != row


def matrix_text(a: Sequence[int], B, cursor: Optional[Tuple[int, int]] = None) -> Text:
    """The sign matrix [diag(a) | B] with column numbers, the cursor cell in brackets."""
    n = len(a)
    t = Text()
    t.append(" " * ROW_LABEL)
    t.append("a".center(n * CELL), style="bold")
    t.append(" ")
    t.append("B".center(n * CELL), style="bold")
    t.append("\n" + " " * ROW_LABEL)
    for c in range(2 * n):
        if c == n:
            t.append(" ")
        t.append(f"{c:^{CELL}}", style="dim")
    for i in range(n):
        t.append(f"\n{i:>{ROW_LABEL - 1}} ", style="dim")
        for c in range(2 * n):
            if c == n:
                t.append("│", style="dim")
            if c < n and c != i:
                cell, st = " · ", "dim"
            else:
                v = a[i] if c < n else B[i][c - n]
                fg, bg = (style.FREE_FG, style.FREE_BG) if v > 0 else (style.OBLIGATE_FG, style.OBLIGATE_BG)
                locked = c >= n and c - n == i
                cell, st = f" {'+' if v > 0 else '−'} ", f"{'' if locked else 'bold '}{fg} on {bg}"
            if cursor == (i, c):
                cell = "[" + cell[1] + "]"
                st += " reverse"
            t.append(cell, style=st)
    return t


class MatrixGrid(Widget, can_focus=True):
    """The editable sign matrix; the pattern lives in the app (``app.a``, ``app.B``)."""

    BINDINGS = [
        Binding("up", "move(-1, 0)", "move", show=False), Binding("down", "move(1, 0)", "move", show=False),
        Binding("left", "move(0, -1)", "move", show=False), Binding("right", "move(0, 1)", "move", show=False),
        Binding("tab,shift+tab", "switch_block", "a ⇄ B", show=False),
        Binding("home", "row_end(0)", "row start", show=False), Binding("end", "row_end(1)", "row end", show=False),
        Binding("space,enter", "flip", "flip"), Binding("plus,equals_sign", "set(1)", "+", show=False),
        Binding("minus", "set(-1)", "−", show=False),
    ]
    DEFAULT_CSS = "MatrixGrid { width: auto; height: auto; padding: 0 1; }"

    def __init__(self) -> None:
        super().__init__()
        self.cursor = (0, 0)

    def render(self) -> Text:
        return matrix_text(self.app.a, self.app.B, self.cursor if self.has_focus else None)

    def get_content_width(self, container, viewport) -> int:
        return ROW_LABEL + 2 * len(self.app.a) * CELL + 1

    def get_content_height(self, container, viewport, width) -> int:
        return 2 + len(self.app.a)

    def action_move(self, dr: int, dc: int) -> None:
        """Move to the next editable cell: the a block holds only its diagonal, the B block skips b_ii."""
        n = len(self.app.a)
        r, c = self.cursor
        if dc:
            for _ in range(2 * n):
                c = (c + dc) % (2 * n)
                if editable(n, r, c):
                    break
        if dr:
            if c < n:                                     # a_i: the diagonal follows the row
                r = (r + dr) % n
                c = r
            else:
                for _ in range(n):
                    r = (r + dr) % n
                    if editable(n, r, c):
                        break
        self.cursor = (r, c)
        self.refresh()

    def action_switch_block(self) -> None:
        """Jump between a_r and the first editable entry of row r of B."""
        n = len(self.app.a)
        r, c = self.cursor
        self.cursor = (r, r) if c >= n else (r, n + (1 if r == 0 else 0))
        self.refresh()

    def action_row_end(self, end: int) -> None:
        """The first (a_r) or last editable cell of the row."""
        n = len(self.app.a)
        r, _ = self.cursor
        self.cursor = (r, r) if not end else (r, 2 * n - 1 if r != n - 1 else 2 * n - 2)
        self.refresh()

    def _write(self, value: Optional[int]) -> None:
        n = len(self.app.a)
        r, c = self.cursor
        if not editable(n, r, c):
            self.app.notify("b_ii is fixed +: every species limits itself" if c >= n
                            else "only the diagonal of diag(a) holds signs: a_i sits at row i, column i",
                            severity="warning", timeout=3)
            return
        if c < n:
            self.app.a[r] = -self.app.a[r] if value is None else value
        else:
            self.app.B[r][c - n] = -self.app.B[r][c - n] if value is None else value
        self.app.changed()

    def action_flip(self) -> None:
        self._write(None)

    def action_set(self, value: int) -> None:
        self._write(value)

    def on_click(self, event: events.Click) -> None:
        n = len(self.app.a)
        x, y = event.x - 1 - ROW_LABEL, event.y - 2         # padding 1, two header lines
        if y < 0 or y >= n or x < 0:
            return
        c = x // CELL if x < n * CELL else (x - 1) // CELL  # one separator column after diag(a)
        if 0 <= c < 2 * n:
            self.cursor = (y, c)
            self.focus()
            self._write(None)

    def on_focus(self) -> None:
        self.refresh()


class NetworkView(Widget):
    """The network of a pattern, sized to the widget."""

    DEFAULT_CSS = "NetworkView { width: 1fr; height: 1fr; min-width: 30; }"

    def __init__(self, pattern=None, **kw) -> None:
        super().__init__(**kw)
        self.pattern = pattern                            # None: the app's pattern

    def render(self) -> Text:
        a, B = self.pattern if self.pattern else (self.app.a, self.app.B)
        w, h = max(self.size.width, 20), max(self.size.height - 3, 8)
        return network_text(a, B, w, h)


# --------------------------------------------------------------------------- screens

class LibraryScreen(ModalScreen):
    """The pattern library, with a preview of the highlighted pattern; enter loads it."""

    BINDINGS = [Binding("escape", "dismiss(None)", "back")]
    DEFAULT_CSS = """
    LibraryScreen { align: center middle; }
    #lib { width: 95%; height: 90%; border: round $accent; background: $surface; }
    #lib OptionList { width: 44; height: 1fr; }
    #preview { width: 1fr; padding: 0 1; }
    #preview-matrix { height: auto; }
    """

    def compose(self) -> ComposeResult:
        self.items = library()
        options, group = [], None
        for k, (g, ident, a, B) in enumerate(self.items):
            if g != group:
                group = g
                options.append(Option(Text(g, style="bold"), disabled=True))
            options.append(Option(f"  {ident:<28} n = {len(a)}", id=str(k)))
        with Horizontal(id="lib"):
            yield OptionList(*options)
            with Vertical(id="preview"):
                yield Static(id="preview-matrix")
                yield NetworkView(pattern=([1, 1], [[1, 1], [1, 1]]), id="preview-net")
        yield Footer()

    def on_mount(self) -> None:
        self.border_title = "library"
        self.query_one("#lib").border_title = "library: ↑↓ preview, enter load, esc back"
        self.query_one(OptionList).focus()

    def on_option_list_option_highlighted(self, event: OptionList.OptionHighlighted) -> None:
        if event.option.id is None:
            return
        g, ident, a, B = self.items[int(event.option.id)]
        self.query_one("#preview-matrix", Static).update(Text(f"{ident}\n", style="bold") + matrix_text(a, B))
        net = self.query_one("#preview-net", NetworkView)
        net.pattern = (a, B)
        net.refresh()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        if event.option.id is not None:
            g, ident, a, B = self.items[int(event.option.id)]
            self.dismiss((ident.split()[0], [int(x) for x in a], [[int(x) for x in r] for r in B]))


class SpeciesScreen(ModalScreen):
    """Choose the number of species."""

    BINDINGS = [Binding("escape", "dismiss(None)", "back")]
    DEFAULT_CSS = """
    SpeciesScreen { align: center middle; }
    SpeciesScreen OptionList { width: 44; height: auto; border: round $accent; }
    """

    def compose(self) -> ComposeResult:
        yield OptionList(*[Option(f"n = {k}" + ("   (2^42 patterns: the search can be slow)" if k == 7 else ""),
                                  id=str(k)) for k in range(2, 8)])

    def on_mount(self) -> None:
        ol = self.query_one(OptionList)
        ol.border_title = "number of species"
        ol.highlighted = max(0, min(5, len(self.app.a) - 2))
        ol.focus()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.dismiss(int(event.option.id))


HELP = """[b]signrules[/b]: the pattern as the sign matrix [diag(a) | B], its network, and the search

[b]matrix[/b]   ←↑↓→ move (only over the signs you can change) · tab jumps a ⇄ B · home/end row ends
         space/enter flip · + / − set · click flips
         a_i is the diagonal of the left block; b_ii is fixed +; b_ij < 0: j benefits i
[b]network[/b]  pink node free (a_i > 0), blue obligate (a_i < 0); a blue dashed arrow i → j:
         i benefits from j; a red line: the pair loses both ways; one-way losses are not drawn
[b]status[/b]   trivial or not, the library pattern it is up to relabelling, and the root verdict
         under all rules (recomputed after every change)

[b]l[/b] library   presets, the eleven n = 4 patterns, the minimal impossible n = 5, 6 patterns
[b]< >[/b]       one species less / more (also [ ]); [b]2[/b]–[b]7[/b] set n directly; [b]n[/b] opens the chooser
         the top-left block is kept; new species are free and compete with everyone
[b]u[/b] rules     the rules to run (space toggles; a all, c R1,R2n) and the node limit
[b]r[/b] run       the search, with its live dashboard; then the verdict (the certificate and its replay for a
         kill at the root). There: h writes the tree viewer (HTML), c the certificate drawing, x stops,
         esc comes back here
[b]e[/b] export    write <name>_network.png / .pdf / .tex (the \\econet call for the article)
[b]s[/b] classification  what the rules decide for n = 4, 5, 6 (the shipped screens)
[b]q[/b] quit      leave (the command line of the pattern is printed)

esc closes this help"""


class HelpScreen(ModalScreen):
    BINDINGS = [Binding("escape,question_mark,q", "dismiss(None)", "back")]
    DEFAULT_CSS = """
    HelpScreen { align: center middle; }
    HelpScreen VerticalScroll { width: 96; height: auto; max-height: 90%; border: round $accent;
                                background: $surface; padding: 1 2; }
    """

    def compose(self) -> ComposeResult:
        with VerticalScroll():
            yield Static(HELP)



LIMITS = [("2,000 nodes", 2000), ("20,000 nodes", 20000), ("200,000 nodes", 200000), ("no limit", None)]


class RulesScreen(ModalScreen):
    """The rules to run and the node limit of the search; esc applies."""

    BINDINGS = [Binding("escape", "done", "done"), Binding("a", "all", "all rules"),
                Binding("c", "classical", "R1,R2n (Annals)")]
    DEFAULT_CSS = """
    RulesScreen { align: center middle; }
    #rules-box { width: 100; height: auto; max-height: 90%; border: round $accent; background: $surface; }
    #rules-box SelectionList { height: auto; border: none; }
    #rules-box OptionList { height: auto; border: none; }
    #rules-box .label { padding: 1 1 0 1; color: $text-muted; }
    """

    def __init__(self, rules, limit) -> None:
        super().__init__()
        self.rules, self.limit = tuple(rules), limit

    def compose(self) -> ComposeResult:
        from textual.widgets import SelectionList
        from .livesearch import RULE_NOTES
        with Vertical(id="rules-box"):
            yield Static("rules (space toggles, tab goes to the limit, esc applies)", classes="label")
            yield SelectionList(*[(f"{r:<4} {RULE_NOTES[r]}", r, r in self.rules)
                                  for r in ("R1", "R2", "R2n", "R3", "R4", "R5", "R6", "R7", "R8")], id="rule-list")
            yield Static("stop the search after", classes="label")
            yield OptionList(*[Option(("● " if v == self.limit else "  ") + t, id=str(k))
                               for k, (t, v) in enumerate(LIMITS)], id="limit-list")
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#rules-box").border_title = "rules and limit"
        self.query_one("#rule-list").focus()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.limit = LIMITS[int(event.option.id)][1]
        ol = self.query_one("#limit-list", OptionList)
        for k, (t, v) in enumerate(LIMITS):
            ol.replace_option_prompt_at_index(k, ("● " if v == self.limit else "  ") + t)

    def _set(self, rules) -> None:
        from textual.widgets import SelectionList
        sl = self.query_one("#rule-list", SelectionList)
        sl.deselect_all()
        for r in rules:
            sl.select(r)

    def action_all(self) -> None:
        from .rules import Ecology
        self._set(Ecology.ALL_RULES)

    def action_classical(self) -> None:
        self._set(("R1", "R2n"))

    def action_done(self) -> None:
        from textual.widgets import SelectionList
        order = ("R1", "R2", "R2n", "R3", "R4", "R5", "R6", "R7", "R8")
        chosen = set(self.query_one("#rule-list", SelectionList).selected)
        self.dismiss((tuple(r for r in order if r in chosen), self.limit))


class _Stopped(Exception):
    pass


class RunScreen(Screen):
    """One search: the live dashboard while it runs, then the verdict; h writes the tree viewer."""

    BINDINGS = [Binding("x", "stop", "stop"), Binding("h", "html", "tree viewer (HTML)"),
                Binding("c", "cert", "certificate drawing"), Binding("escape,b", "back", "back to the pattern")]
    DEFAULT_CSS = """
    #dash { height: auto; padding: 0 1; }
    #verdict { height: 1fr; padding: 0 1; }
    """

    def __init__(self, a, B, name: str, rules, limit) -> None:
        super().__init__()
        from .livesearch import LiveSearch
        from .rules import Ecology
        self.eco = Ecology(a, B, name)
        self.rules, self.limit = tuple(rules), limit
        self.tally = LiveSearch(self.eco, self.rules)
        self.res = None
        self.stop_requested = False

    def compose(self) -> ComposeResult:
        yield Static(id="dash")
        with VerticalScroll(id="verdict"):
            yield Static(Text("searching… x stops", style="dim"), id="verdict-text")
        yield Footer()

    def on_mount(self) -> None:
        import time
        self.tally.t0 = time.time()
        self._timer = self.set_interval(0.15, self._tick)
        self.search()

    def _tick(self) -> None:
        self.query_one("#dash", Static).update(self.tally)

    def _progress(self, ev) -> None:
        if self.stop_requested:
            raise _Stopped()
        self.tally(ev)

    @work(thread=True, exclusive=True)
    def search(self) -> None:
        from .search import complete
        try:
            res = complete(self.eco, use=self.rules, progress=self._progress, record_tree=True,
                           keep_derivations=True, max_nodes=self.limit)
        except _Stopped:
            res = None
        panel = None
        if res is not None:
            from .report import verdict_panel
            panel = verdict_panel(self.eco, res, self.rules, closed=self.tally.closed, covered=self.tally.covered)
        self.app.call_from_thread(self._finished, res, panel)

    def _finished(self, res, panel) -> None:
        self.res = res
        self.tally.finished = True
        self._timer.stop()
        self._tick()
        box = self.query_one("#verdict-text", Static)
        if res is None:
            t = self.tally
            box.update(Text(f"stopped by you after {t.nodes:,} nodes: {t.completions:,} completions, "
                            f"{t.conflicts:,} contradictions, {100 * t.covered:.2f}% swept.  esc goes back.", style="bold"))
        else:
            box.update(panel)

    def action_stop(self) -> None:
        if self.res is None:
            self.stop_requested = True

    def _stem(self) -> str:
        return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in (self.eco.name or "pattern"))

    def action_html(self) -> None:
        if self.res is None:
            self.notify("wait for the search to end (x stops it)", severity="warning")
            return
        from .treeview import save_tree_viewer
        path = f"{self._stem()}_tree.html"
        save_tree_viewer(self.eco, self.res, path, use=self.rules)
        self.notify(f"wrote {path}: open it in a browser", timeout=8)

    def action_cert(self) -> None:
        if self.res is None:
            self.notify("wait for the search to end (x stops it)", severity="warning")
            return
        from .certificate import write_certificate, write_completion
        stem = self._stem()
        try:
            if self.res.root_conflict:
                write_certificate(self.eco, self.res, f"{stem}_certificate")
                self.notify(f"wrote {stem}_certificate.png / .pdf / .dot", timeout=8)
            elif self.res.derivations:
                write_completion(self.eco, self.res, f"{stem}_completion0", index=0)
                self.notify(f"wrote {stem}_completion0.png / .pdf / .dot (the first completion)", timeout=8)
            else:
                self.notify("nothing to draw: no kill at the root, no completion", severity="warning")
        except Exception as e:                            # Graphviz missing, typically
            self.notify(f"drawing failed: {e} (needs Graphviz's dot)", severity="error")

    def action_back(self) -> None:
        self.stop_requested = True
        self.app.pop_screen()


class ClassificationScreen(ModalScreen):
    """What the rules decide: the shipped results of the screens, n = 4, 5, 6."""

    BINDINGS = [Binding("escape,s", "dismiss(None)", "back"), Binding("l", "library", "browse the kills")]
    DEFAULT_CSS = """
    ClassificationScreen { align: center middle; }
    ClassificationScreen VerticalScroll { width: 100; height: auto; max-height: 90%; border: round $accent;
                                          background: $surface; padding: 1 2; }
    """

    def compose(self) -> ComposeResult:
        from collections import Counter
        from rich.table import Table
        from .data import screen_records
        tb = Table(box=None, padding=(0, 2), header_style="bold")
        for col in ("n", "candidates", "killed at the root", "by rule", "the rest"):
            tb.add_column(col)
        tb.add_row("4", "11 (Meng et al.)", "11", "R2, R3, R7, R8", "—")
        for n in (5, 6):
            recs = screen_records(n)
            kills = [r for r in recs if r["status"] == "killed_root"]
            rules = Counter((r.get("info") or "?").split()[0] for r in kills)
            tb.add_row(str(n), str(len(recs)), str(len(kills)), ", ".join(f"{r} {c}" for r, c in sorted(rules.items())),
                       f"{len(recs) - len(kills)} possible")
        intro = Text("Every impossible pattern on n species deletes to an impossible one on n − 1 (the deletion "
                     "lemma), so the candidates at n are the one-species extensions of the kills at n − 1.  The "
                     "rules kill the candidates below at the root; the others have certificates of possibility, so "
                     "these kills are exactly the minimal impossible patterns (the certificates of possibility, tropical witnesses, are in "
                     "the repository feasible_ecologies).\n", style="")
        foot = Text("\nl opens the library, where every kill can be loaded (\"minimal impossible\" groups); "
                    "esc goes back.  Recomputing the screens takes minutes: signrules screen 5 / 6 in a shell.",
                    style="dim")
        with VerticalScroll():
            yield Static(intro)
            yield Static(tb)
            yield Static(foot)
        yield Footer()

    def on_mount(self) -> None:
        self.query_one(VerticalScroll).border_title = "classification"

    def action_library(self) -> None:
        self.dismiss("library")


# --------------------------------------------------------------------------- the app

def free_competitors(n: int) -> Tuple[List[int], List[List[int]]]:
    """n free species, every pair competing: the blank pattern to start from."""
    return [1] * n, [[1] * n for _ in range(n)]


class PatternEditor(App):
    """The editor; ``self.return_value`` is (a, B, name, action) or None."""

    TITLE = "signrules"
    BINDINGS = [
        Binding("l", "library", "library"), Binding("n", "species", "species"),
        Binding("less_than_sign,left_square_bracket", "grow(-1)", "n−1", key_display="<"),
        Binding("greater_than_sign,right_square_bracket", "grow(1)", "n+1", key_display=">"),
        *[Binding(str(k), f"set_species({k})", f"n = {k}", show=False) for k in range(2, 8)],
        Binding("u", "rules", "rules"), Binding("r", "run", "run"),
        Binding("e", "export", "export"), Binding("s", "classification", "classification"),
        Binding("question_mark", "help", "help"), Binding("q", "leave('quit')", "quit"),
    ]
    CSS = """
    #top { height: 1fr; }
    #left { width: auto; }
    #matrix-box { width: auto; height: auto; border: round $accent; }
    #network-box { border: round $accent; padding: 0 1; }
    #hint { padding: 0 1; width: auto; }
    #status { height: auto; padding: 0 1; border-top: solid $accent; }
    """

    def __init__(self, a=None, B=None, name: str = "", n: int = 4, rules=None, limit: Optional[int] = 20000,
                 open_library: bool = False) -> None:
        super().__init__()
        from .patterns import preset_rules
        from .rules import Ecology
        self.rules = tuple(rules) if rules else (preset_rules(name) or tuple(Ecology.ALL_RULES))
        self.limit = limit
        self.open_library = open_library
        if a is None:
            a, B = free_competitors(n)
        self.a = [int(x) for x in a]
        self.B = [[int(x) for x in r] for r in B]
        self.origin = name                                # the name of the pattern as loaded
        self.loaded = (list(self.a), [list(r) for r in self.B])
        self._generation = 0

    EDITOR_ACTIONS = {"library", "species", "grow", "set_species", "rules", "run", "export", "classification",
                      "help", "leave"}

    def check_action(self, action: str, parameters) -> Optional[bool]:
        """The editor's keys act only on the main screen (not under the run screen or a dialog)."""
        if action in self.EDITOR_ACTIONS and len(self.screen_stack) > 1:
            return False
        return True

    @property
    def name_now(self) -> str:
        if self.origin and (self.a, self.B) == self.loaded:
            return self.origin
        return f"edited from {self.origin}" if self.origin else "custom"

    def compose(self) -> ComposeResult:
        with Horizontal(id="top"):
            with Vertical(id="left"):
                with Vertical(id="matrix-box"):
                    yield MatrixGrid()
                yield Static(Text("←↑↓→ move · tab a ⇄ B · home/end · space flip · + − set\n"
                                  "< > or 2–7 change n · l library · u rules · r run · ? help", style="dim"), id="hint")
            with Vertical(id="network-box"):
                yield NetworkView()
        yield Static(id="status")
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#network-box").border_title = "network"
        self.query_one(MatrixGrid).focus()
        self.changed()
        if self.open_library:
            self.action_library()

    # ------------------------------------------------------------------ updates
    def changed(self) -> None:
        """Redraw after an edit and recompute the status in the background."""
        grid = self.query_one(MatrixGrid)
        grid.refresh(layout=True)
        self.query_one(NetworkView).refresh()
        self.query_one("#matrix-box").border_title = f"pattern · n = {len(self.a)}"
        self.query_one("#matrix-box").border_subtitle = self.name_now
        self._generation += 1
        self.query_one("#status", Static).update(self._status_text(None))
        self.analyse(self._generation, list(self.a), [list(r) for r in self.B], self.rules)

    def _status_text(self, info) -> Text:
        from .patterns import trivial
        t = Text()
        n = len(self.a)
        t.append(f"n = {n} · ", style="dim")
        if trivial(self.a, self.B):
            t.append("trivial", style=f"bold {style.POS}")
            t.append(": an obligate species without a benefactor, impossible for that reason alone\n")
        else:
            t.append("not trivial\n", style="dim")
        if info is None:
            t.append("checking…", style="dim")
            return t
        matches, verdict, open_signs, seconds = info
        t.append("library: ", style="dim")
        t.append(", ".join(matches) + " (up to relabelling)" if matches else "not in the library")
        t.append("\nroot (under the rules below): ", style="dim")
        if verdict:
            t.append("IMPOSSIBLE", style="bold")
            t.append(f" · killed by {verdict}")
        else:
            t.append(f"{open_signs} signs open after root propagation", style=style.COMPLETION)
            t.append(" · r runs the search", style="dim")
        t.append(f"  ({seconds:.2f}s)", style="dim")
        t.append("\nsearch: ", style="dim")
        t.append(f"rules {', '.join(self.rules) or 'none'} · ")
        t.append(f"stop after {self.limit:,} nodes" if self.limit else "no node limit")
        t.append("   u changes them · r runs", style="dim")
        return t

    @work(thread=True, exclusive=True)
    def analyse(self, generation: int, a, B, rules) -> None:
        from . import labels as L
        from .rules import Ecology
        from .search import complete
        matches = library_matches(a, B)
        res = complete(Ecology(a, B), use=rules, max_nodes=0, max_completions=1, keep_derivations=False)
        verdict = L.rule(res.root_conflict).headline("plain") if res.root_conflict else None
        info = (matches, verdict, res.root_unknowns, res.seconds)
        if generation == self._generation:
            self.call_from_thread(self._show, generation, info)

    def _show(self, generation, info) -> None:
        if generation == self._generation:
            self.query_one("#status", Static).update(self._status_text(info))

    # ------------------------------------------------------------------ actions
    def action_library(self) -> None:
        def load(picked) -> None:
            if picked:
                ident, a, B = picked
                from .patterns import preset_rules
                from .rules import Ecology
                self.a, self.B, self.origin = a, B, ident
                self.loaded = (list(a), [list(r) for r in B])
                own = preset_rules(ident)
                if own or self.rules == ("R1", "R2n"):
                    self.rules = own or tuple(Ecology.ALL_RULES)
                    self.notify(f"rules {', '.join(self.rules)}" + (": the preset's own" if own else ""), timeout=3)
                self.query_one(MatrixGrid).cursor = (0, 0)
                self.changed()
        self.push_screen(LibraryScreen(), load)

    def action_set_species(self, n: int) -> None:
        """Resize to n species: the top-left block is kept, new species are free and compete with everyone."""
        if not n or n == len(self.a):
            return
        if not 2 <= n <= 7:
            self.notify("between 2 and 7 species", severity="warning", timeout=2)
            return
        fa, fB = free_competitors(n)
        m = min(n, len(self.a))
        for i in range(m):
            fa[i] = self.a[i]
            for j in range(m):
                fB[i][j] = self.B[i][j]
        self.a, self.B = fa, fB
        grid = self.query_one(MatrixGrid)
        r, c = grid.cursor
        grid.cursor = (r, c) if r < n and c < n and c == r else (0, 0)
        self.changed()

    def action_grow(self, step: int) -> None:
        self.action_set_species(len(self.a) + step)

    def action_species(self) -> None:
        self.push_screen(SpeciesScreen(), self.action_set_species)

    def action_rules(self) -> None:
        def apply(picked) -> None:
            if picked:
                self.rules, self.limit = picked
                self.changed()
        self.push_screen(RulesScreen(self.rules, self.limit), apply)

    def action_run(self) -> None:
        self.push_screen(RunScreen(self.a, self.B, self.name_now, self.rules, self.limit))

    def action_classification(self) -> None:
        def then(what) -> None:
            if what == "library":
                self.action_library()
        self.push_screen(ClassificationScreen(), then)

    def action_help(self) -> None:
        self.push_screen(HelpScreen())

    def action_export(self) -> None:
        from .network import save_network
        stem = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in self.name_now) + "_network"
        try:
            paths = save_network(self.a, self.B, stem, title=self.name_now)
            self.notify("wrote " + ", ".join(paths), timeout=6)
        except Exception as e:                            # a read-only directory, typically
            self.notify(f"export failed: {e}", severity="error")

    def action_leave(self, action: str) -> None:
        self.exit((list(self.a), [list(r) for r in self.B], self.name_now, action))


def command_line(a, B) -> str:
    """``signrules run --a=... --B=...`` for the pattern."""
    sa = "".join("+" if v > 0 else "-" for v in a)
    sB = "/".join("".join("+" if v > 0 else "-" for v in r) for r in B)
    return f"signrules run --a={sa} --B={sB}"


def edit(a=None, B=None, name: str = "", n: int = 4, open_library: bool = False):
    """Run the app; (a, B, name, "quit") when left with q, or None (closed otherwise)."""
    return PatternEditor(a, B, name, n, open_library=open_library).run()
