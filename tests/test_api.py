"""The Python API a user meets first: Ecology.parse / repr / str, SearchResult.verdict, minor_signs, complete()."""
import inspect

import numpy as np
import pytest

from signrules import (Ecology, complete, ELEVEN, parse_pattern, parse_signs, presets, resolve_preset,
                       patterns)
from signrules import report

S_D = "S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"
A_E = "A-e competitor + 3-cycle (Thm 3.5)"


def test_parse_accepts_the_three_string_forms_and_sequences():
    a, B = ELEVEN[A_E]
    a_str, B_str = patterns.compact_pattern(a, B)                  # "+---", "+-++/... "
    by_strings = Ecology.parse(a_str, B_str)
    assert by_strings.a == list(a) and (by_strings.B == np.asarray(B)).all()
    assert Ecology.parse("+,-,-,-", "1,-1,1,1/-1,1,1,1/1,1,1,-1/1,1,-1,1").a == by_strings.a
    assert (Ecology.parse("1,-1,-1,-1", "+-++/-+++/+++-/++-+").B ==
            Ecology.parse("+---", "+-++/-+++/+++-/++-+").B).all()
    assert Ecology.parse(a, B).a == list(a)                       # lists and arrays pass through
    assert (Ecology.parse(a, B_str).B == np.asarray(B)).all()     # mixed
    assert Ecology.parse("+---", B, name="A-e").name == "A-e"
    with pytest.raises(ValueError):
        Ecology.parse("+-", "++/+")
    with pytest.raises(ValueError):
        Ecology.parse("+-", "-+/++")                               # diagonal not +


def test_repr_round_trips_and_str_is_the_sign_matrix():
    eco = Ecology.parse("+---", "+-++/-+++/+++-/++-+", name="A-e")
    assert repr(eco) == "Ecology.parse('+---', '+-++/-+++/+++-/++-+', name='A-e')"
    back = eval(repr(eco))
    assert back.a == eco.a and (back.B == eco.B).all() and back.name == "A-e"
    assert repr(Ecology.parse("+--", "+++/-++/--+")) == "Ecology.parse('+--', '+++/-++/--+')"
    text = str(eco)
    lines = text.split("\n")
    assert lines[0] == "A-e"
    assert "│" in lines[3] and lines[3].rstrip().endswith("+ − + +")      # row 0 of B
    assert lines[-1].rstrip().endswith("+ + − +")
    assert patterns.compact_pattern(eco.a, eco.B) == ("+---", "+-++/-+++/+++-/++-+")


def test_verdict_three_values_and_agrees_with_the_report():
    a, B = ELEVEN[A_E]
    killed = complete(Ecology(a, B))
    assert killed.verdict == "impossible" and killed.root_conflict
    assert repr(killed).startswith("SearchResult(impossible: impossible at the root")
    a, B = ELEVEN[S_D]
    eco = Ecology(a, B)
    open_ = complete(eco, use=("R1", "R2n"))
    assert open_.verdict == "open" and len(open_.completions) == 64
    assert repr(open_) == "SearchResult(open: 64 completions, 131 nodes, 4 pruned)"
    cut = complete(eco, use=("R1", "R2n"), max_nodes=1)
    assert cut.verdict == "undecided" and cut.stopped and not cut.completions
    closed = complete(Ecology(*ELEVEN["S-a obligate mutualism"]), use=("R1", "R2n"))
    assert closed.verdict == "impossible" and not closed.root_conflict and closed.conflicts
    words = {"impossible": "IMPOSSIBLE", "open": "NOT DECIDED BY THE RULES", "undecided": "UNDECIDED"}
    for res in (killed, open_, cut, closed):
        assert report.verdict(res)[0] == words[res.verdict]
    for name, (a, B) in presets().items():
        if len(a) == 4:
            res = complete(Ecology(a, B), max_nodes=200)
            assert report.verdict(res)[0] == words[res.verdict]


def test_minor_signs_names_a_completion():
    a, B = ELEVEN[S_D]
    eco = Ecology(a, B)
    res = complete(eco, use=("R1", "R2n"), max_completions=1)
    m = eco.minor_signs(res.completions[0])
    assert m["det B"] == 1                                         # R1
    assert set(m.values()) <= {1, -1}
    assert len(m) == len(eco.idx.bases) - 1                        # every minor but the empty one
    assert m["det B_{3},{0}"] == int(B[3][0])                     # 1×1 minors are the entries
    partial = eco.minor_signs(eco.initial_table())
    assert 0 < len(partial) < len(m)


def test_complete_documents_the_two_groups():
    params = inspect.signature(complete).parameters
    assert "seed" in params                     # random branch order (random completions or coarse cells)
    doc = inspect.cleandoc(complete.__doc__)
    split = doc.index("For the ablations")
    assert doc.index("Parameters") < split
    for k in ("use", "max_nodes", "max_completions", "record_tree", "progress", "keep_derivations"):
        assert doc.index(f"\n{k} :") < split
    for k in ("idle", "root_only", "stats", "engine", "branching", "probe", "symmetry", "decide", "assume", "seed"):
        assert doc.index(f"\n{k} :") > split


def test_top_level_exports_the_parsers():
    assert parse_signs("+--") == [1, -1, -1]
    assert parse_pattern("+-", "++/-+")[0] == [1, -1]
    assert resolve_preset("A-e")
    assert A_E in presets()
