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


def test_flip_set_and_locked_cells():
    seen = {}
    out = drive(PatternEditor(*AE, "A-e"),
                "right", "right", "right", "right", "right", "space",     # (0, 5): b_01 + -> -
                "down", "left", "left", "left", "left", "minus",          # (1, 1): a_1 stays -
                "up", "space",                                            # (0, 1): not a sign
                "right", "right", "right", "space",                       # (0, 4): b_00 fixed
                "left", "left", "left", "left", "space",                  # (0, 0): a_0 + -> -
                lambda app: seen.update(name=app.name_now), "q")
    a, B, name, action = out
    assert action == "quit" and name == "edited from A-e" == seen["name"]
    assert a == [-1, -1, -1, -1] and B[0] == [1, -1, 1, 1] and all(B[i][i] == 1 for i in range(4))


def test_status_follows_the_pattern():
    seen = {}
    drive(PatternEditor(*AE, "A-e"), 1.0, lambda app: seen.update(before=status(app)),
          "right", "right", "right", "right", "down", "right", "space", 1.0,   # (1, 5): b_11 is fixed, nothing
          lambda app: seen.update(same=status(app)),
          "right", "space", 1.0, lambda app: seen.update(after=status(app)), "q")  # (1, 6): b_12 - -> +
    assert "A-e" in seen["before"] and "killed by R7" in seen["before"]
    assert seen["same"] == seen["before"]
    assert "A-e" not in seen["after"].split("root:")[0]


def test_library_loads_and_species_resize():
    out = drive(PatternEditor(), "l", 0.3, "down", "down", "down", "enter", 0.5, "q")   # S-a highlighted first
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
    from signrules.patterns import PATTERNS
    assert len(library()) == len(ELEVEN) + len(PATTERNS) + 24 + 41
    assert "A-e" in library_matches(*AE)
