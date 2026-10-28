"""The live search dashboard: its counters agree with the search result."""
import io


from rich.console import Console

from signrules import Ecology, complete
from signrules.livesearch import LiveSearch, _milestone


def test_counters_match_result_and_sweep_is_complete():
    from signrules.patterns import ELEVEN
    a, B = list(ELEVEN.values())[3]
    eco = Ecology(a, B, "d")
    live = LiveSearch(eco, ("R1", "R2n"))
    res = complete(eco, use=("R1", "R2n"), progress=live, keep_derivations=False)
    assert live.completions == len(res.completions)
    assert live.conflicts == res.conflicts
    assert abs(live.covered - 1.0) < 1e-12               # an exhaustive search sweeps everything
    assert sum(live.closed.values()) == res.conflicts
    out = io.StringIO()
    Console(file=out, width=120).print(live)            # renders without error, with the gloss under the header
    assert "a node" in out.getvalue() and "completion" in out.getvalue()


def test_root_conflict_and_limits():
    from signrules.patterns import ELEVEN
    a, B = list(ELEVEN.values())[0]
    eco = Ecology(a, B)
    live = LiveSearch(eco)
    res = complete(eco, progress=live)
    assert res.root_conflict and live.conflicts == 1 and live.covered == 1.0
    eco = Ecology([1, -1, -1], [[1, 1, 1], [-1, 1, 1], [-1, -1, 1]])
    live = LiveSearch(eco, ("R1", "R2n"))
    res = complete(eco, use=("R1", "R2n"), progress=live, max_nodes=3)
    assert res.stopped and live.stopped > 0 and live.covered + live.abandoned <= 1.0 + 1e-12


def test_milestones():
    assert [k for k in range(1, 120) if _milestone(k)] == [1, 2, 5, 10, 20, 50, 100]


def test_random_pattern_is_reproducible_with_positive_diagonal():
    from signrules.patterns import random_pattern
    a, B = random_pattern(5, seed=3)
    assert (a, B) == random_pattern(5, seed=3)
    assert len(a) == 5 and all(B[i][i] == 1 for i in range(5))
    assert {x for r in B for x in r} <= {-1, 1}
