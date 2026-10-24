"""signrules edit: the keys change the right entries, the library loads, the status follows."""
import asyncio

import pytest

pytest.importorskip("textual")

from signrules.editor import PatternEditor, command_line, editable, library, library_matches  # noqa: E402
from signrules.patterns import ELEVEN  # noqa: E402

AE = ELEVEN["A-e competitor + 3-cycle (Thm 3.5)"]


def drive(app, *steps, size=(130, 34)):
    """Run ``app`` headless; each step is a key (str), a pause (float) or a callable(app)."""
    async def go():
        async with app.run_test(size=size) as pilot:
            await pilot.pause(0.5)
            for s in steps:
                if isinstance(s, float):
                    await pilot.pause(s)
                elif callable(s):
                    s(app)
                else:
                    await pilot.press(s)
    asyncio.run(go())
    return app.return_value


def status(app) -> str:
    return str(app.query_one("#status").render())


def cursor(app):
    from signrules.editor import MatrixGrid
    return app.query_one(MatrixGrid).cursor


def put(r, c):
    """A step that places the cursor (as a click would, without flipping)."""
    def go(app):
        from signrules.editor import MatrixGrid
        app.query_one(MatrixGrid).cursor = (r, c)
    return go


def test_navigation_visits_only_editable_cells():
    seen = []
    note = lambda app: seen.append(cursor(app))                            # noqa: E731
    drive(PatternEditor(*AE, "A-e"),
          "right", note, "right", note, "right", note, "right", note,        # a_0, then b_01 b_02 b_03, wrap to a_0
          "down", note, "down", note,                                        # a_1, a_2: the diagonal follows
          "tab", note, "tab", note,                                          # row 2 of B: b_20 (col 4); back to a_2
          "end", note, "home", note,                                         # b_23 (col 7); a_2
          put(0, 5), "down", note,                                           # b_01 -> b_11 is fixed -> b_21
          "q")
    assert seen == [(0, 5), (0, 6), (0, 7), (0, 0), (1, 1), (2, 2), (2, 4), (2, 2), (2, 7), (2, 2), (2, 5)]


def test_flip_set_and_locked_cells():
    seen = {}
    out = drive(PatternEditor(*AE, "A-e"),
                "right", "space",                                            # (0, 5): b_01 + -> -
                "down", "minus",                                             # (2, 5): b_21 + -> -
                put(1, 5), "space",                                          # b_11: fixed, a warning
                put(0, 1), "space",                                          # not a sign
                put(0, 0), "space",                                          # a_0 + -> -
                lambda app: seen.update(name=app.name_now), "q")
    a, B, name, action = out
    assert action == "quit" and name == "edited from A-e" == seen["name"]
    assert a == [-1, -1, -1, -1] and B[0] == [1, -1, 1, 1] and B[2][1] == -1
    assert all(B[i][i] == 1 for i in range(4))


def test_status_follows_the_pattern():
    seen = {}
    drive(PatternEditor(*AE, "A-e"), 1.0, lambda app: seen.update(before=status(app)),
          put(1, 6), "space", 1.0, lambda app: seen.update(after=status(app)), "q")   # b_12 - -> +
    assert "A-e" in seen["before"] and "killed by R7" in seen["before"]
    assert "A-e" not in seen["after"].split("root:")[0]


def test_change_n_with_keys():
    seen = []
    note = lambda app: seen.append((len(app.a), app.query_one("#matrix-box").border_title))   # noqa: E731
    out = drive(PatternEditor(*AE, "A-e"), ">", note, "greater_than_sign", note, "<", note, "3", note,
                "right_square_bracket", note, "2", "<", note, "7", ">", note, "q")
    assert [k for k, _ in seen] == [5, 6, 5, 3, 4, 2, 7]                     # stays within 2 … 7
    assert "n = 5" in seen[0][1]
    a, B, _, _ = out
    assert len(a) == 7 and a[:2] == list(AE[0][:2]) and B[0][:2] == list(map(int, AE[1][0][:2]))


def test_library_loads_and_species_resize():
    from signrules.patterns import TREES
    downs = ["down"] * (len(TREES) + 3)                    # the first entry is highlighted; S-d is the 4th of ELEVEN
    out = drive(PatternEditor(), "l", 0.3, *downs, "enter", 0.5, "q")
    sd = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]
    assert out[2] == "S-d" and out[0] == list(sd[0]) and out[1] == [list(map(int, r)) for r in sd[1]]
    out = drive(PatternEditor(*AE, "A-e"), "n", 0.3, "down", "enter", 0.3, "q")      # n = 4 -> 5
    a, B, _, _ = out
    assert len(a) == 5 and a[:4] == list(AE[0]) and B[4] == [1] * 5 and B[0][:4] == list(map(int, AE[1][0]))


def test_leave_actions_and_helpers():
    assert drive(PatternEditor(*AE, "A-e"), "r")[3] == "run"
    assert drive(PatternEditor(*AE, "A-e"), "t")[3] == "tree"
    assert drive(PatternEditor(*AE, "A-e"), "question_mark", 0.2, "escape", "q")[3] == "quit"
    assert command_line(*AE) == "signrules run --a=+--- --B=++++/++-+/+++-/+-++"
    assert editable(3, 0, 0) and not editable(3, 0, 1) and not editable(3, 1, 4) and editable(3, 1, 3)
    from signrules.patterns import PATTERNS, TREES
    assert len(library()) == len(TREES) + len(ELEVEN) + len(PATTERNS) + 24 + 41
    assert "A-e" in library_matches(*AE)
