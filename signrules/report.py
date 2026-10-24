"""
The verdict of one search, printed at the end of ``signrules run`` and ``signrules replay``.

* killed at the root: the certificate, one numbered row per forced sign (the sign, the minor it is,
  the rule, the rows it used), ending in the contradiction, and the result of the independent
  replay (``replay.replay_kill``);
* every branch closed: impossible, with the contradictions per rule;
* completions found: not decided by the rules (a completion is a sign table consistent with every
  rule, not a realisation);
* stopped by a limit before any completion: undecided.

    from signrules.report import report
    report(eco, complete(eco), Ecology.ALL_RULES)
"""
from collections import Counter
from typing import Optional, Sequence

from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from . import labels as L
from . import style
from .certificate import certificate_steps
from .replay import replay_kill
from .rules import Ecology
from .search import SearchResult


def _ordered(steps):
    """The certificate steps in an order where every step follows the steps it uses."""
    produced = {S for S, *_ in steps}
    done, out, pending = set(), [], list(steps)
    while pending:
        rest = [st for st in pending if not all(p in done or p not in produced for p in st[3])]
        ready = [st for st in pending if st not in rest]
        if not ready:                       # cannot happen for a well-formed record; keep the given order
            out += rest
            break
        out += ready
        done |= {S for S, *_ in ready}
        pending = rest
    return out


def certificate_table(eco: Ecology, res: SearchResult, verbose: bool = False) -> Table:
    """The root certificate of ``res`` as a table: #, sign, minor, rule, the rows it used."""
    steps, premises, _ = certificate_steps(eco, res)
    steps = _ordered(steps)
    row_of = {S: k for k, (S, *_) in enumerate(steps, 1)}
    tb = Table(box=None, padding=(0, 1), header_style="dim", show_edge=False)
    for col, just in (("#", "right"), ("sign", "left"), ("minor", "left"), ("rule", "left"), ("uses", "left")):
        tb.add_column(col, justify=just, no_wrap=col in ("#", "sign", "minor"))

    def uses(prem):
        rows = sorted(row_of[p] for p in prem if p in row_of)
        entries = sum(p not in row_of for p in prem)
        parts = [" ".join(f"#{r}" for r in rows)] if rows else []
        if entries:
            parts.append(f"{entries} {'entry' if entries == 1 else 'entries'}")
        return " + ".join(parts) or "—"

    for k, (S, sign, why, prem) in enumerate(steps, 1):
        lab = L.rule(why, eco.n)
        rule = Text(lab.headline("plain"), style=f"bold {style.rule_color(lab.rule)}")
        for extra in lab.lines("plain")[1:]:
            if extra != L.minor(eco, S, sign, "plain"):
                rule.append("\n" + extra, style="dim")
        tb.add_row(str(k), Text(L.plucker(S, sign, "plain"), style=style.POS if sign > 0 else style.NEG),
                   L.minor(eco, S, sign, "plain"), rule, uses(prem))
    lines = L.conflict(eco, res.root_conflict, "plain")
    lab = L.rule(res.root_conflict, eco.n)
    rule = Text(lab.headline("plain"), style=f"bold {style.rule_color(lab.rule)}")
    detail = [t for t in lines[2:] if verbose or len(t) <= 70]
    if len(detail) < len(lines[2:]):
        detail.append("(the full clause: --verbose)")
    for extra in detail:
        rule.append("\n" + extra, style="dim")
    tb.add_row("", Text("⊥", style="bold"), "contradiction", rule, uses(premises))
    return tb


def verdict(res: SearchResult):
    """(word, explanation) for the outcome of ``res``."""
    if res.root_conflict:
        return "IMPOSSIBLE", f"killed at the root by {L.rule(res.root_conflict).headline('plain')}"
    if res.completions:
        more = " (the search stopped early: there may be more)" if res.stopped else ""
        k = len(res.completions)
        return ("NOT DECIDED BY THE RULES",
                f"{k:,} completion{'s' if k != 1 else ''}{more}: sign tables consistent with every rule used, so the "
                "rules cannot kill this pattern; it may still be impossible, or have a realisation")
    if res.stopped:
        return "UNDECIDED", "the search stopped at its limit before finding a completion; raise --max-nodes"
    return "IMPOSSIBLE", f"every branch closed: {res.conflicts:,} contradictions in {res.nodes:,} nodes"


def report(eco: Ecology, res: SearchResult, use: Sequence[str], console: Optional[Console] = None,
           closed: Optional[Counter] = None, covered: Optional[float] = None, verbose: bool = False,
           replay: bool = True) -> None:
    """Print the verdict panel of ``res`` (see the module docstring).  ``closed`` (contradictions per rule)
    and ``covered`` (share of the search space swept) come from a ``livesearch.LiveSearch`` tally."""
    (console or Console()).print(verdict_panel(eco, res, use, closed, covered, verbose, replay))


def verdict_panel(eco: Ecology, res: SearchResult, use: Sequence[str], closed: Optional[Counter] = None,
                  covered: Optional[float] = None, verbose: bool = False, replay: bool = True) -> Panel:
    """The verdict of ``report`` as a rich renderable (the TUI shows it in a widget)."""
    word, why = verdict(res)
    colour = style.COMPLETION if word.startswith("NOT") else ("" if word == "UNDECIDED" else "bold")
    body = [Text.assemble((word, f"bold {colour}".strip()), "  ·  ", why)]
    body.append(Text(f"{eco.name + ' · ' if eco.name else ''}{L.pattern(eco, 'plain')} · rules {', '.join(use) or 'none'}"
                     f" · {res.seconds:.2f}s", style="dim"))
    if res.root_conflict:
        body += [Text(""), certificate_table(eco, res, verbose)]
        if replay:
            rep = replay_kill(eco, res, use)
            body += [Text(""), Text.assemble(("replay ", "dim"), ("OK" if rep.ok else "FAILED",
                                                                   f"bold {style.COMPLETION if rep.ok else style.POS}"),
                                             "  " + (f"{rep.steps} step{'s' if rep.steps != 1 else ''} and the contradiction"
                                                    if rep.steps else "the contradiction") +
                                             (" re-derived by the independent checker (signrules.replay)" if rep.ok
                                              else f": {rep.message}"))]
    else:
        stats = [f"{res.nodes:,} nodes", f"{res.conflicts:,} contradictions"]
        if res.root_unknowns is not None:
            stats.insert(0, f"{res.root_unknowns} signs open after root propagation")
        if covered is not None and res.stopped:
            stats.append(f"{100 * min(covered, 1.0):.2f}% of the tree swept")
        body.append(Text(" · ".join(stats), style="dim"))
        if closed:
            body.append(Text("contradictions by rule: " + ", ".join(f"{r} {c:,}" for r, c in closed.most_common()),
                             style="dim"))
    return Panel(Group(*body), title="verdict", title_align="left", border_style="dim")
