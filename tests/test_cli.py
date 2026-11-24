"""The command ``signrules``: help, pattern input, errors and the verdicts it prints."""
import argparse
import subprocess
import sys

import pytest

from signrules import cli
from signrules.patterns import ELEVEN, PATTERNS, parse_pattern, resolve_preset


def run(capsys, *argv):
    code = cli.main(list(argv))
    return code, capsys.readouterr()


def test_top_level_help_lists_every_command(capsys):
    with pytest.raises(SystemExit):
        cli.main(["--help"])
    out = capsys.readouterr().out
    for name, (summary, _, _) in cli.COMMANDS.items():
        assert name in out and summary.split(";")[0][:40] in " ".join(out.split())


def test_no_command_without_terminal_prints_help(capsys):
    code, cap = run(capsys)
    assert code == 0 and "commands:" in cap.out


def test_every_subcommand_and_option_is_described():
    ap = cli.build_parser(full=True)
    sub = next(a for a in ap._actions if isinstance(a, argparse._SubParsersAction))
    assert set(sub.choices) == set(cli.COMMANDS) == set(cli.HANDLERS)
    for name, sp in sub.choices.items():
        assert sp.description, name
        assert "examples:" in sp.epilog, name
        for action in sp._actions:
            assert action.help and action.help != argparse.SUPPRESS, \
                f"signrules {name} {action.option_strings or action.dest}: no help"


def test_ablation_options_hidden_from_help_but_accepted(capsys):
    """--idle, --root-only, --name are listed by --help-all only; they still parse."""
    for flag, hidden in (("--help", False), ("--help-all", True)):
        with pytest.raises(SystemExit) as e:
            cli.main(["run", flag])
        out = capsys.readouterr().out
        assert e.value.code == 0
        assert ("--idle IDLE" in out) == hidden and ("--root-only ROOT_ONLY" in out) == hidden
    code, cap = run(capsys, "run", "--preset", "A-e", "--idle", "", "--root-only", "R4", "--name", "x", "--plain")
    assert code == 0 and "IMPOSSIBLE" in cap.out


def test_default_node_limit_and_lifting_it(capsys):
    """run stops at 20,000 nodes by default and says so; --max-nodes 0 lifts the limit."""
    code, cap = run(capsys, "run", "--preset", "5cycle", "--plain")
    assert code == 0 and "limit 20,000 nodes" in cap.out and "stopped early" in cap.out
    code, cap = run(capsys, "run", "--preset", "S-d", "--rules", "annals", "--max-nodes", "0", "--plain")
    assert code == 0 and "no node limit" in cap.out and "64 completions" in cap.out and "stopped" not in cap.out


def test_rules_annals_is_r1_r2n(capsys):
    _, annals = run(capsys, "run", "--preset", "S-d", "--rules", "annals", "--plain")
    _, listed = run(capsys, "run", "--preset", "S-d", "--rules", "R1,R2n", "--plain")
    assert "64 completions" in annals.out and "131 nodes" in annals.out
    assert annals.out.replace("annals", "R1,R2n") == listed.out or "64 completions" in listed.out


def test_leading_minus_in_a_or_b_needs_no_equals(capsys):
    code, cap = run(capsys, "run", "--a", "-+-", "--B", "+++/-++/--+", "--plain")
    assert code == 0 and "sign a = (−,+,−)" in cap.out
    with pytest.raises(SystemExit):                       # the value reaches the parser: a pattern error, not argparse's
        cli.main(["run", "--a", "+-", "--B", "-+/++", "--plain"])
    assert "diagonal of B is positive" in capsys.readouterr().err
    assert cli._join_sign_values(["--a", "-1,-1", "--B", "--plain"]) == ["--a=-1,-1", "--B", "--plain"]


def test_presets_one_line_each(capsys):
    _, short = run(capsys, "presets")
    _, long = run(capsys, "presets", "--long")
    assert all(len(line) <= 100 for line in short.out.splitlines())
    assert len(long.out.splitlines()) > len(short.out.splitlines())
    assert "64 completions, 131 nodes" in " ".join(long.out.split())


def test_run_kill_prints_certificate_and_replay(capsys):
    code, cap = run(capsys, "run", "--preset", "A-e", "--plain")
    assert code == 0
    assert "IMPOSSIBLE" in cap.out and "killed at the root by R7" in cap.out
    assert "contradiction" in cap.out and "p_4567 > 0" in cap.out
    assert "replay OK" in cap.out
    assert "independent checker" in cap.out          # the one-sentence meaning of the verdict


def test_run_completions(capsys):
    code, cap = run(capsys, "run", "--preset", "S-d", "--rules", "R1,R2n", "--plain")
    assert code == 0 and "NOT DECIDED BY THE RULES" in cap.out and "64 completions" in cap.out
    assert "may have a realisation" in cap.out


def test_run_exhausted_and_undecided(capsys):
    # S-a under R1,R2n survives the root and dies in the search (every branch closed)
    _, cap = run(capsys, "run", "--preset", "S-a", "--rules", "R1,R2n", "--plain")
    assert "IMPOSSIBLE" in cap.out and "every branch closed" in cap.out and "contradictions by rule" in cap.out
    _, cap = run(capsys, "run", "--preset", "S-a", "--rules", "R1,R2n", "--max-nodes", "1", "--plain")
    assert "UNDECIDED" in cap.out and "node limit" in cap.out


@pytest.mark.parametrize("argv, message", [
    (["run", "--preset", "S"], "ambiguous: S-a, S-b, S-c, S-d, S-e"),
    (["run", "--preset", "zzz"], "no preset 'zzz'"),
    (["run", "--preset", "A-e", "--rules", "R9"], "unknown rule R9"),
    (["run", "--a", "+--", "--B", "++/++"], "B needs 3 rows of 3 signs"),
    (["run", "--a", "+-", "--B=-+/++"], "diagonal of B is positive"),
    (["run", "--a", "+2", "--B", "++/++"], "'2' is none of these"),
    (["run"], "give a pattern"),
    (["run", "--preset", "A-e", "--random", "4"], "exactly one pattern"),
    (["replay", "--preset", "A-e", "--idle", "R0"], "unknown rule R0"),
])
def test_bad_input_is_a_clean_error(capsys, argv, message):
    with pytest.raises(SystemExit) as e:
        cli.main(argv)
    assert e.value.code == 2
    err = capsys.readouterr().err
    assert message in err and "Traceback" not in err


def test_compact_and_comma_signs_agree():
    a1, B1 = parse_pattern("+--", "+++/-++/--+")
    a2, B2 = parse_pattern("1,-1,-1", "1,1,1/-1,1,1/-1,-1,1")
    assert a1 == a2 and (B1 == B2).all()
    assert parse_pattern("+,−,−", "+,+,+/−,+,+/−,−,+")[0] == a2


def test_presets_resolve_by_short_id_and_prefix(capsys):
    assert resolve_preset("a-e") == "A-e competitor + 3-cycle (Thm 3.5)"
    assert resolve_preset("parity") == "parity_k2"
    code, cap = run(capsys, "presets")
    assert code == 0
    for name in list(ELEVEN) + list(PATTERNS):
        assert name.split()[0] in cap.out


def test_replay(capsys):
    code, cap = run(capsys, "replay", "--preset", "parity_k2")
    assert code == 0 and "R8" in cap.out and "replay OK" in cap.out
    code, cap = run(capsys, "replay", "--preset", "S-d", "--rules", "R1,R2n")
    assert code == 1 and "not killed at the root" in cap.out


def test_tree_writes_the_viewer(capsys, tmp_path):
    out = tmp_path / "t.html"
    code, cap = run(capsys, "tree", "--preset", "S-d", "--rules", "R1,R2n", "--out", str(out))
    assert code == 0 and "64 completions" in cap.out and out.read_text().startswith("<!doctype html>")


def test_screen_candidates(capsys):
    code, cap = run(capsys, "screen", "5", "--candidates-only")
    assert code == 0 and "83 non-trivial candidates" in cap.out


def test_module_and_aliases_run():
    for cmd in ([sys.executable, "-m", "signrules", "--help"],
                [sys.executable, "-c", "from signrules.livesearch import main; main(['--preset', 'A-e', '--plain'])"],
                [sys.executable, "-c", "from signrules.treeview import main; main(['--help'])"]):
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        assert r.returncode == 0, r.stderr


def test_readme_command_table_matches_commands():
    from pathlib import Path
    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text()
    for name, (summary, _, _) in cli.COMMANDS.items():
        assert f"| `{name}` | {summary} |" in readme, f"README: the row of {name} differs from cli.COMMANDS"


def test_tree_presets_grow_trees(capsys):
    """Every tree_* preset survives the root under its rules and the search ends (the claim of its note)."""
    from signrules import Ecology, complete
    from signrules.patterns import TREES, preset_rules
    for name, (a, B) in TREES.items():
        res = complete(Ecology(a, B), use=preset_rules(name) or Ecology.ALL_RULES, keep_derivations=False)
        assert res.root_conflict is None and not res.stopped and res.conflicts > 0, name
    code, cap = run(capsys, "run", "--preset", "tree_S-a", "--plain")
    assert "the preset's own" in cap.out and "every branch closed" in cap.out
    code, cap = run(capsys, "run", "--preset", "tree_S-d", "--rules", "all", "--plain")
    assert "the preset's own" not in cap.out and "killed at the root" in cap.out


def test_keyboard_interrupt_keeps_the_partial_result():
    from signrules import Ecology, complete
    a, B = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]

    def interrupt(event, seen=[]):
        seen.append(event)
        if len(seen) == 40:
            raise KeyboardInterrupt
    res = complete(Ecology(a, B), use=("R1", "R2n"), progress=interrupt)
    assert res.interrupted and res.stopped and 0 < len(res.completions) < 64
    assert "stopped early" in res.summary()


def test_version(capsys):
    import signrules
    from importlib.metadata import version
    assert signrules.__version__ == version("signrules")
    with pytest.raises(SystemExit):
        cli.main(["--version"])
    assert capsys.readouterr().out.strip() == f"signrules {signrules.__version__}"
