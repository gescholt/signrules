"""
``signrules edit``: edit a sign pattern as the matrix [diag(a) | B] and watch its network.

    signrules edit                    # four free competitors to start from
    signrules edit --n 5
    signrules edit --preset A-e

The left panel is the n × 2n sign matrix: the diagonal block holds the signs of a (only its diagonal is
editable), the right block holds B (its diagonal is fixed +).  Column numbers are those of [diag(a) | B],
so the Plücker coordinate p_0257 is the minor on the columns 0, 2, 5, 7.  The right panel draws the
interaction network as the article's ``\\econet`` (``network.py``).  Below them, recomputed after every
change: whether the pattern is trivial, which library pattern it is up to relabelling, and the verdict of
the root propagation under all rules.

Keys: the arrows move over the editable signs, tab jumps between a and B, space or enter flips, + and - set a sign (a click flips too); l opens the
library, < and > (or 2–7, or n for a chooser) change the number of species, r runs the search, t writes the tree viewer, e exports the
network (png, pdf and the \\econet call), ? shows this help, q quits and prints the command line.

Needs textual (the ``live`` extra).  ``edit()`` returns the edited pattern and the action chosen.
"""
from typing import Dict, List, Optional, Sequence, Tuple

from rich.text import Text
from textual import events, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import ModalScreen
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


HELP = """[b]signrules edit[/b]: the pattern as the sign matrix [diag(a) | B]

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
[b]r[/b] run       leave and run `signrules run` on the pattern (dashboard, verdict)
[b]t[/b] tree      leave and write the tree viewer (`signrules tree`)
[b]e[/b] export    write <name>_network.png / .pdf / .tex (the \\econet call for the article)
[b]q[/b] quit      leave and print the command line of the pattern

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


# --------------------------------------------------------------------------- the app

def free_competitors(n: int) -> Tuple[List[int], List[List[int]]]:
    """n free species, every pair competing: the blank pattern to start from."""
    return [1] * n, [[1] * n for _ in range(n)]


class PatternEditor(App):
    """The editor; ``self.return_value`` is (a, B, name, action) or None."""

    TITLE = "signrules edit"
    BINDINGS = [
        Binding("l", "library", "library"), Binding("n", "species", "species"),
        Binding("less_than_sign,left_square_bracket", "grow(-1)", "n−1", key_display="<"),
        Binding("greater_than_sign,right_square_bracket", "grow(1)", "n+1", key_display=">"),
        *[Binding(str(k), f"set_species({k})", f"n = {k}", show=False) for k in range(2, 8)],
        Binding("r", "leave('run')", "run"), Binding("t", "leave('tree')", "tree"),
        Binding("e", "export", "export"), Binding("question_mark", "help", "help"),
        Binding("q", "leave('quit')", "quit"),
    ]
    CSS = """
    #top { height: 1fr; }
    #left { width: auto; }
    #matrix-box { width: auto; height: auto; border: round $accent; }
    #network-box { border: round $accent; padding: 0 1; }
    #hint { padding: 0 1; width: auto; }
    #status { height: auto; padding: 0 1; border-top: solid $accent; }
    """

    def __init__(self, a=None, B=None, name: str = "", n: int = 4) -> None:
        super().__init__()
        if a is None:
            a, B = free_competitors(n)
        self.a = [int(x) for x in a]
        self.B = [[int(x) for x in r] for r in B]
        self.origin = name                                # the name of the pattern as loaded
        self.loaded = (list(self.a), [list(r) for r in self.B])
        self._generation = 0

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
                                  "< > or 2–7 change n · l library · ? help", style="dim"), id="hint")
            with Vertical(id="network-box"):
                yield NetworkView()
        yield Static(id="status")
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#network-box").border_title = "network"
        self.query_one(MatrixGrid).focus()
        self.changed()

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
        self.analyse(self._generation, list(self.a), [list(r) for r in self.B])

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
        t.append("\nroot: ", style="dim")
        if verdict:
            t.append("IMPOSSIBLE", style="bold")
            t.append(f" · killed by {verdict}")
        else:
            t.append(f"{open_signs} signs open after root propagation", style=style.COMPLETION)
            t.append(" · r runs the search", style="dim")
        t.append(f"  ({seconds:.2f}s)", style="dim")
        t.append("\n" + command_line(self.a, self.B), style="dim")
        return t

    @work(thread=True, exclusive=True)
    def analyse(self, generation: int, a, B) -> None:
        from . import labels as L
        from .rules import Ecology
        from .search import complete
        matches = library_matches(a, B)
        res = complete(Ecology(a, B), max_nodes=0, max_completions=1, keep_derivations=False)
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
                self.a, self.B, self.origin = a, B, ident
                self.loaded = (list(a), [list(r) for r in B])
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


def edit(a=None, B=None, name: str = "", n: int = 4):
    """Run the editor; (a, B, name, action) with action "run", "tree" or "quit", or None (closed otherwise)."""
    return PatternEditor(a, B, name, n).run()
