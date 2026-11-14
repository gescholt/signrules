"""The interaction network drawn as the \\econet of the papers: edges, the macro call, the drawings."""
import pytest

from signrules.network import econet_edges, econet_tex, network_text, save_network
from signrules.patterns import ELEVEN, random_pattern


def from_econet(signs, edges):
    """The pattern of an \\econet call: every listed i/j/b is b_ij < 0, every other entry a loss (+)."""
    n = len(signs)
    a = [1 if s == "p" else -1 for s in signs]
    B = [[1] * n for _ in range(n)]
    for e in edges:
        i, j, t = e.split("/")
        if t == "b":
            B[int(i)][int(j)] = -1
    return a, B


# four three-species networks: two obligate cycles and two patterns with a free species
EX_THREE = [("m,m,m", "0/1/b,1/2/b,2/0/b"), ("m,m,m", "0/1/b,1/0/b,1/2/b,2/1/b,2/0/b,0/2/b"),
            ("p,m,m", "1/2/b,2/1/b,0/1/c,0/2/c"), ("p,m,m", "1/2/b,2/1/b,0/1/b,0/2/b")]


@pytest.mark.parametrize("signs, edges", EX_THREE)
def test_econet_round_trip(signs, edges):
    a, B = from_econet(signs.split(","), edges.split(","))
    out = econet_tex(a, B)
    assert out.startswith(f"\\econet{{{signs}}}{{")             # n = 3: the macro's default, no [n]
    assert set(out.split("}{")[1].rstrip("}").split(",")) == set(edges.split(","))


def test_edges_follow_the_papers():
    a, B = ELEVEN["A-e competitor + 3-cycle (Thm 3.5)"]          # free competitor 0 over the obligate cycle 1 -> 2 -> 3
    assert set(econet_edges(a, B)) == {(1, 2, "b"), (2, 3, "b"), (3, 1, "b"), (0, 1, "c"), (0, 2, "c"), (0, 3, "c")}
    assert econet_tex(a, B).startswith("\\econet[4]{p,m,m,m}")
    a = [1, 1]
    assert econet_edges(a, [[1, -1], [1, 1]]) == [(0, 1, "b")]  # predation: one arrow, the loss not drawn


@pytest.mark.parametrize("n", range(2, 8))
def test_network_text_renders(n):
    a, B = random_pattern(n, seed=n)
    t = network_text(a, B, 40, 14)
    assert len(t.plain.splitlines()) == 16 and all(str(i) in t.plain for i in range(n))


def test_save_network(tmp_path):
    a, B = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]
    paths = save_network(a, B, str(tmp_path / "sd"))
    assert [p.rsplit(".", 1)[1] for p in paths] == ["png", "pdf", "tex"]
    assert (tmp_path / "sd.tex").read_text().strip() == econet_tex(a, B)
