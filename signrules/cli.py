"""
The command ``signrules``: one entry point for the sign completion.

    signrules                         # in a terminal: menus for everything
    signrules run --preset A-e        # dashboard, then the verdict (and the certificate of a kill)
    signrules presets                 # the named patterns and their short ids
    signrules <command> --help

The commands and their one-line descriptions are ``COMMANDS``; the help, the menu and the README take
them from there.  ``signrules-live``, ``signrules-tree`` and ``signrules-screen`` are the same as
``signrules run``, ``signrules tree`` and ``signrules screen``.
"""
import argparse
import sys
from typing import Optional, Sequence

from .rules import Ecology

DESCRIPTION = ("Decide whether a Lotka–Volterra sign pattern is impossible: complete the signs of the Plücker "
               "coordinates of [diag(a) | B] with the Grassmann–Plücker relations and the rules R1–R8.")

# name: (summary, description, examples)
COMMANDS = {
    "run": ("complete one pattern with a live dashboard; prints the verdict and, for a kill, the certificate",
            "Complete the sign table of one pattern: propagation at the root, then a depth-first search over the "
            "open signs.  A terminal shows a live dashboard (progress, rules, current path); at the end comes the "
            "verdict: impossible (with the certificate and its independent replay when the kill is at the root), "
            "not decided by the rules (completions exist), or undecided (a limit was hit).",
            "signrules run --preset A-e                        # a four-species kill, certificate printed\n"
            "signrules run --preset S-d --rules R1,R2n         # 64 completions under the Annals algorithm\n"
            "signrules run --a +--- --B \"+-++/-+++/+++-/++-+\"  # your own pattern"),
    "tree": ("complete one pattern and write an interactive HTML viewer of the search tree",
             "Run the search with the tree recorded and write one self-contained HTML page: the binary search "
             "tree, and for every node the signs it forces, the rule behind each, stability and feasibility.",
             "signrules tree --preset S-d --rules R1,R2n --out s-d.html\n"
             "signrules tree --random 5 --seed 7 --max-nodes 3000"),
    "replay": ("kill one pattern at the root and check its certificate with the independent replayer",
               "Run the root propagation only.  When it ends in a contradiction, print the certificate and replay "
               "every step from the identity or inequality it names (signrules.replay), without the engine.",
               "signrules replay --preset A-d\n"
               "signrules replay --preset parity_k2 --cert out/"),
    "edit": ("the full-screen app (also plain `signrules`): library, matrix editor, network, rules, search, verdict",
             "Open the pattern editor: the n × 2n sign matrix [diag(a) | B] on the left (arrows move, space flips), "
             "the interaction network on the right in the style of the article (pink free, blue obligate, dashed "
             "blue arrow i → j where i benefits from j, red line where both lose), and below them the root verdict "
             "and the library pattern it matches, recomputed after every change.  l opens the library, n changes "
             "the number of species, r runs the search, e exports the network.  Needs the live extra (textual).",
             "signrules edit\n"
             "signrules edit --n 5\n"
             "signrules edit --preset A-e"),
    "presets": ("list the named patterns (the eleven n = 4 patterns, the n = 5, 6 presets) with short ids",
                "List the patterns that --preset accepts.  A preset is named by its short id (S-a … S-e, A-a … A-f, "
                "parity_k2, ...) or any unambiguous prefix of its name.",
                "signrules presets"),
    "screen": ("list the candidates for impossible n-species patterns (deletion lemma) and run the rules on each",
               None,                                          # the docstring of signrules.screen
               "signrules screen 6 --candidates-only\n"
               "signrules screen 5 --workers 12 --out out/n5.jsonl"),
}

KNOWN_RULES = tuple(Ecology.ALL_RULES) + ("R2n",)


# --------------------------------------------------------------------------- shared options

def add_pattern_args(ap: argparse.ArgumentParser) -> None:
    """--preset, --a/--B, --random/--seed, --name."""
    g = ap.add_argument_group("pattern (one of --preset, --a with --B, --random)")
    g.add_argument("--preset", help="a named pattern: short id (S-d, A-e, parity_k2) or a prefix of its name; "
                                    "'signrules presets' lists them")
    g.add_argument("--a", help="signs of the growth rates: 1,-1,-1 or +,-,- or +--  "
                               "(write --a=-1,... when the first sign is negative, as for --B)")
    g.add_argument("--B", help="signs of B, rows separated by /: 1,1,1/-1,1,1/-1,-1,1 or +++/-++/--+ "
                               "(diagonal +; B[i,j] < 0: j benefits i)")
    g.add_argument("--random", type=int, metavar="N", help="a uniformly random pattern on N species")
    g.add_argument("--seed", type=int, default=None, help="seed for --random (printed when omitted)")
    g.add_argument("--name", default="", help="a label for the pattern in the output")


def add_rule_args(ap: argparse.ArgumentParser) -> None:
    """--rules, --idle, --root-only."""
    g = ap.add_argument_group("rules")
    g.add_argument("--rules", default=None,
                   help=f"'all' (R1–R8, the default; a tree_* preset brings its own) or a comma list of "
                        f"{', '.join(KNOWN_RULES)}; R1,R2n is the algorithm of the Annals paper")
    g.add_argument("--idle", default="auto",
                   help="rules among --rules that are registered but never run: 'auto' (default: R4 and R6 "
                        "when the rules that subsume them are on), a list such as R4,R6, or '' to run everything")
    g.add_argument("--root-only", default="", help="rules among --rules to run at the root only, such as R4,R6")


def _rule_list(text: str, ap: argparse.ArgumentParser, what: str = "--rules") -> tuple:
    names = tuple(r.strip() for r in text.split(",") if r.strip())
    bad = [r for r in names if r not in KNOWN_RULES]
    if bad:
        ap.error(f"{what}: unknown rule {', '.join(bad)}; the rules are {', '.join(KNOWN_RULES)}")
    return names


def rules_from_args(args, ap: argparse.ArgumentParser):
    """(use, idle, root_only) from the parsed rule options; a bad name is an argparse error."""
    rules = args.rules
    if rules is None:                                     # not given: the preset's own rules, else all
        from .patterns import preset_rules, resolve_preset
        own = None
        if getattr(args, "preset", None):
            try:
                own = preset_rules(resolve_preset(args.preset))
            except ValueError:
                pass
        if own:
            print(f"rules {','.join(own)}: the preset's own (--rules all to override)")
            rules = ",".join(own)
        else:
            rules = "all"
    use = tuple(Ecology.ALL_RULES) if rules == "all" else _rule_list(rules, ap)
    idle = "auto" if args.idle == "auto" else _rule_list(args.idle, ap, "--idle")
    return use, idle, _rule_list(args.root_only, ap, "--root-only")


def pattern_from_args(args, ap: argparse.ArgumentParser) -> Ecology:
    """The Ecology named by the parsed pattern options; bad input is an argparse error (exit code 2)."""
    from .patterns import parse_pattern, presets, random_pattern, resolve_preset
    given = [bool(args.preset), bool(args.a or args.B), bool(args.random)]
    if sum(given) != 1:
        ap.error("give exactly one pattern: --preset NAME, --a with --B, or --random N")
    if args.random:
        if args.seed is None:
            import random
            args.seed = random.randrange(10 ** 6)
        a, B = random_pattern(args.random, args.seed)
        print(f"random pattern (seed {args.seed}):  --a={','.join(map(str, a))} "
              f"--B={'/'.join(','.join(map(str, r)) for r in B)}")
        return Ecology(a, B, args.name or f"random n={args.random}, seed {args.seed}")
    if args.preset:
        try:
            full = resolve_preset(args.preset)
        except ValueError as e:
            ap.error(str(e))
        a, B = presets()[full]
        return Ecology(a, B, args.name or full)
    if not (args.a and args.B):
        ap.error("--a and --B go together")
    try:
        a, B = parse_pattern(args.a, args.B)
    except ValueError as e:
        ap.error(str(e))
    return Ecology(a, B, args.name)


def _interactive() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def _fill_from_menu(args, picked) -> None:
    args.preset, args.a, args.B, args.random = picked["preset"], picked["a"], picked["B"], picked["random"]


# --------------------------------------------------------------------------- the commands

def cmd_run(args, ap) -> int:
    from rich.console import Console
    from .livesearch import LiveSearch, _pick, _presets
    from .report import report
    from .search import complete
    picked_rules = None
    if not (args.preset or args.a or args.B or args.random):
        if not _interactive():
            ap.error("give a pattern: --preset NAME, --a with --B, or --random N")
        picked = _pick(_presets())
        if picked is None:
            return 1
        _fill_from_menu(args, picked)
        picked_rules = picked["rules"]
        args.max_completions = args.max_completions or picked["max_completions"]
        args.max_nodes = args.max_nodes or picked["max_nodes"]
        args.html = args.html or picked["html"]
    eco = pattern_from_args(args, ap)
    use, idle, root_only = rules_from_args(args, ap)
    if picked_rules is not None:
        use = picked_rules
    kw = dict(use=use, record_tree=bool(args.html), keep_derivations=bool(args.cert), max_nodes=args.max_nodes,
              max_completions=args.max_completions, root_only=root_only, idle=idle)
    console = Console()
    tally = LiveSearch(eco, use)
    if args.plain or not sys.stdout.isatty():
        res = complete(eco, progress=tally, **kw)
    else:
        try:
            with tally:
                res = complete(eco, progress=tally, **kw)
        except KeyboardInterrupt:
            console.print(f"\ninterrupted after {tally.nodes:,} nodes: {tally.completions:,} completions, "
                          f"{tally.conflicts:,} contradictions, {100 * tally.covered:.2f}% swept")
            return 130
        console.print()
    report(eco, res, use, console, closed=tally.closed, covered=tally.covered, verbose=args.verbose)
    _outputs(eco, res, use, args, console)
    return 0


def _outputs(eco, res, use, args, console) -> None:
    """--html and --cert."""
    if getattr(args, "html", None):
        from .treeview import save_tree_viewer
        save_tree_viewer(eco, res, args.html, use=use)
        console.print(f"tree viewer written to {args.html}")
    if getattr(args, "cert", None):
        from pathlib import Path
        from .certificate import write_certificate, write_completion
        out = Path(args.cert)
        out.mkdir(parents=True, exist_ok=True)
        stem = out / ("".join(c if c.isalnum() or c in "-_" else "_" for c in (eco.name or "pattern").split()[0]))
        try:
            if res.root_conflict:
                write_certificate(eco, res, f"{stem}_certificate")
                console.print(f"certificate written to {stem}_certificate.png / .pdf / .dot")
            elif res.derivations:
                write_completion(eco, res, f"{stem}_completion0", index=0)
                console.print(f"first completion written to {stem}_completion0.png / .pdf / .dot")
            else:
                console.print("--cert: nothing to draw (no kill at the root, no completion)")
        except Exception as e:                            # Graphviz missing, typically
            console.print(f"--cert failed: {e} (the drawings need Graphviz's dot on the PATH)")


def cmd_tree(args, ap) -> int:
    from .livesearch import _pick_pattern, _presets
    from .search import complete
    from .treeview import save_tree_viewer
    if not (args.preset or args.a or args.B or args.random):
        if not _interactive():
            ap.error("give a pattern: --preset NAME, --a with --B, or --random N")
        picked = _pick_pattern(_presets())
        if picked is None:
            return 1
        _fill_from_menu(args, picked)
    eco = pattern_from_args(args, ap)
    use, idle, root_only = rules_from_args(args, ap)
    res = complete(eco, use=use, record_tree=True, max_nodes=args.max_nodes, idle=idle, root_only=root_only)
    save_tree_viewer(eco, res, args.out, use=use)
    print(f"{res.summary()}  ->  {args.out}")
    return 0


def cmd_replay(args, ap) -> int:
    from rich.console import Console
    from .livesearch import _pick_pattern, _presets
    from .report import report
    from .search import complete
    if not (args.preset or args.a or args.B or args.random):
        if not _interactive():
            ap.error("give a pattern: --preset NAME, --a with --B, or --random N")
        picked = _pick_pattern(_presets())
        if picked is None:
            return 1
        _fill_from_menu(args, picked)
    eco = pattern_from_args(args, ap)
    use, idle, root_only = rules_from_args(args, ap)
    res = complete(eco, use=use, max_nodes=0, max_completions=1, keep_derivations=False, idle=idle)
    console = Console()
    if res.root_conflict is None:
        console.print(f"{eco.name or 'pattern'}: not killed at the root by {', '.join(use) or 'propagation alone'} "
                      f"({res.root_unknowns} signs open after root propagation); nothing to replay. "
                      "'signrules run' searches the branches.")
        return 1
    report(eco, res, use, console, verbose=args.verbose)
    _outputs(eco, res, use, args, console)
    return 0


def cmd_presets(args, ap) -> int:
    from rich.console import Console
    from rich.table import Table
    from .livesearch import _group, _signs
    from .patterns import presets, short_id
    table = presets()
    groups: dict = {}
    for name in table:
        groups.setdefault(_group(name), []).append(name)
    width = max(len(short_id(k)) for k in table)
    console = Console()
    for group, names in groups.items():
        console.print(f"[bold]{group}[/bold]")
        tb = Table(box=None, padding=(0, 2), show_header=False)
        tb.add_column(min_width=width, no_wrap=True); tb.add_column(no_wrap=True); tb.add_column(no_wrap=True)
        tb.add_column(style="dim")
        for name in names:
            a, _ = table[name]
            from .patterns import TREE_NOTES, preset_rules
            note = TREE_NOTES.get(name, name[len(short_id(name)):].strip())
            if preset_rules(name):
                note += f"  [rules {','.join(preset_rules(name))}]"
            tb.add_row(short_id(name), f"n = {len(a)}", f"a = {_signs(a)}", note)
        console.print(tb)
    console.print("[dim]use the id (or a prefix of the name):  signrules run --preset tree_S-d_annals[/dim]")
    return 0


def cmd_edit(args, ap) -> int:
    if not _interactive():
        ap.error("the editor needs a terminal (stdin and stdout must be a tty)")
    try:
        from .editor import command_line, edit
    except ImportError:
        ap.error("the editor needs textual:  pip install -e \".[live]\"  (or pip install textual)")
    if args.preset or args.a or args.B or args.random:
        from .patterns import short_id
        eco = pattern_from_args(args, ap)
        name = short_id(eco.name) if args.preset and not args.name else eco.name
        out = edit(list(eco.a), [list(map(int, r)) for r in eco.B], name)
    else:
        if not 2 <= args.n <= 7:
            ap.error("--n: between 2 and 7 species")
        out = edit(n=args.n)
    if out is None:
        return 1
    a, B, name, action = out
    print(command_line(a, B))
    return 0


def edited_pattern() -> Optional[dict]:
    """The pattern of the editor as the menus' {"preset", "a", "B", "random"}; None without textual or on quit."""
    try:
        from .editor import edit
    except ImportError:
        return None
    out = edit()
    if out is None:
        return None
    a, B, _, _ = out
    return {"preset": None, "random": None, "a": ",".join(map(str, a)), "B": "/".join(",".join(map(str, r)) for r in B)}


def cmd_screen(args, ap) -> int:
    from .screen import run_args
    run_args(args)
    return 0


HANDLERS = {"run": cmd_run, "tree": cmd_tree, "replay": cmd_replay, "edit": cmd_edit, "presets": cmd_presets,
            "screen": cmd_screen}


# --------------------------------------------------------------------------- the parser

def build_parser() -> argparse.ArgumentParser:
    """The parser of ``signrules`` with one subparser per entry of ``COMMANDS``."""
    from . import screen
    ap = argparse.ArgumentParser(
        prog="signrules", description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="With no command, in a terminal: the full-screen app (pattern library, matrix editor, rules, "
               "search, verdict; every action is a key).  Run 'signrules <command> --help' for the options of a command.")
    sub = ap.add_subparsers(dest="command", title="commands", metavar="<command>")
    for name, (summary, description, examples) in COMMANDS.items():
        sp = sub.add_parser(name, help=summary, description=description or screen.__doc__,
                            epilog="examples:\n  " + examples.replace("\n", "\n  "),
                            formatter_class=argparse.RawDescriptionHelpFormatter)
        sp.set_defaults(_parser=sp)
        if name in ("run", "tree", "replay"):
            add_pattern_args(sp)
            add_rule_args(sp)
        if name == "edit":
            add_pattern_args(sp)
            sp.add_argument("--n", type=int, default=4, help="number of species to start from, without a pattern "
                                                           "(default 4: free competitors)")
        if name in ("run", "replay"):
            out = sp.add_argument_group("output")
            out.add_argument("--cert", metavar="DIR",
                             help="write the certificate of a kill (or the derivation of the first completion) "
                                  "as .png, .pdf, .dot to DIR (needs Graphviz)")
            out.add_argument("--verbose", action="store_true", help="print the full clause of the contradiction")
        if name == "run":
            out.add_argument("--html", metavar="FILE", help="also write the interactive tree viewer to FILE")
            out.add_argument("--plain", action="store_true",
                             help="no dashboard, only the verdict (automatic when the output is not a terminal)")
            g = sp.add_argument_group("search")
            g.add_argument("--max-nodes", type=int, default=None, help="stop after this many nodes (default: no limit)")
            g.add_argument("--max-completions", type=int, default=None,
                           help="stop after this many completions (default: no limit; 1 is enough to show the "
                                "rules cannot kill the pattern)")
        if name == "tree":
            g = sp.add_argument_group("output")
            g.add_argument("--max-nodes", type=int, default=20000,
                           help="stop after this many nodes (default 20000, so the page stays small)")
            g.add_argument("--out", default="tree.html", help="the HTML file (default tree.html)")
        if name == "screen":
            screen.add_arguments(sp)
    return ap


def _menu(ap) -> int:
    """'What do you want to do?' with the commands of ``COMMANDS``, then that command's own menus."""
    from .livesearch import questionary_kit
    kit = questionary_kit()
    if kit is None:
        ap.print_help()
        return 0
    questionary, Choice, Separator, kw, kt = kit
    width = max(map(len, COMMANDS))
    choice = questionary.select("What do you want to do?",
                                choices=[Choice(f"{n:<{width}}  {s}", value=n) for n, (s, _, _) in COMMANDS.items()],
                                **kw).ask()
    if choice is None:
        return 1
    if choice == "screen":
        n = questionary.select("Species", choices=["5", "6"], default="6", **kw).ask()
        what = questionary.select("Run", choices=[Choice("count the candidates only", value=True),
                                                  Choice("the rules on every candidate (minutes)", value=False)],
                                  **kw).ask()
        if n is None or what is None:
            return 1
        argv = ["screen", n] + (["--candidates-only"] if what else ["--out", f"out/n{n}_extension_screen.jsonl"])
    else:
        argv = [choice]
    print(f"$ signrules {' '.join(argv)}")
    return main(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point of ``signrules`` (and ``python -m signrules``)."""
    ap = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        if _interactive():
            try:
                from .editor import command_line, edit
            except ImportError:                           # no textual: the question-and-answer menus
                return _menu(ap)
            out = edit(open_library=True)
            if out:
                print(command_line(out[0], out[1]))
            return 0
        ap.print_help()
        return 0
    args = ap.parse_args(argv)
    return HANDLERS[args.command](args, args._parser) or 0


def run_alias(argv: Optional[Sequence[str]] = None) -> int:
    """``signrules-live``."""
    return main(["run", *(sys.argv[1:] if argv is None else argv)])


def tree_alias(argv: Optional[Sequence[str]] = None) -> int:
    """``signrules-tree``."""
    return main(["tree", *(sys.argv[1:] if argv is None else argv)])
