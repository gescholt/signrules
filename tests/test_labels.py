"""The figure text is generated from the run: parser and formats of signrules/labels.py."""

import sys
from pathlib import Path

import numpy as np


from signrules.patterns import ELEVEN                                    # noqa: E402
from signrules import Ecology, complete, labels as L                 # noqa: E402
from signrules.certificate import certificate_diagram, search_tree_diagram   # noqa: E402


def test_plucker_and_minor_labels():
    eco = Ecology([1, -1, 1, -1], np.ones((4, 4), dtype=int))
    S = (0, 2, 5, 7)                       # T = {0, 2}, P = Q = {1, 3}, eps = -1
    assert L.plucker(S, 1) == "$p_{0257} > 0$"
    assert L.plucker(S, -1, "html") == "p<sub>0257</sub> &lt; 0"
    assert L.plucker(S, 1, "text") == "p_0257 > 0"
    assert L.dictionary(eco, S) == r"$p_{0257} = -a_{0}a_{2}\,\det B_{\{1,3\}}$"
    assert L.dictionary(eco, (4, 5, 6, 7)) == r"$p_{4567} = \det B$"
    # p_S = -a_0 a_2 det B_{13} with a_0 a_2 > 0: a positive coordinate is a negative minor
    assert L.minor(eco, S, 1) == r"$\det B_{\{1,3\}} < 0$"


def test_every_logged_rule_string_parses():
    seen = set()
    for name, (a, B) in ELEVEN.items():
        eco = Ecology(a, B, name)
        res = complete(eco, engine="sweep", idle=())   # the sweep reaches R4 on these (run it: it is idle by default)
        for why, *_ in res.root_log:
            lab = L.rule(why, eco.n)
            assert lab.rule in {"R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"}, why
            seen.add(lab.rule)
            for fmt in ("tex", "html", "text"):
                assert all(isinstance(t, str) and t for t in lab.lines(fmt, with_clause=True))
        lines = L.conflict(eco, res.root_conflict)
        assert lines[0] == r"$\bot$" and lines[1].split(" · ")[0].startswith("R")
    assert {"R1", "R4", "R5"} <= seen


def test_rule_parameters_are_typeset():
    lab = L.rule("R4 Z-block P=(0,) Q=(1, 2, 3) s=(1, -1, -1) σ=(2, 3, 1): M-matrix, principal minor > 0", 4)
    assert lab.headline() == "R4 · Z-block"
    assert lab.detail == "$P=\\{0\\},\\ Q=\\{1,2,3\\}$\n$s=(+,-,-),\\ \\sigma=(2\\,3\\,1)$"
    assert lab.lines("text")[1] == "P = {0}, Q = {1,2,3}\ns = (+,−,−), σ = (2 3 1)"
    assert L.rule("R5 Laplace expansion of det B_(0, 1, 2, 3),(0, 1, 2, 3)", 4).detail == r"$\det B$"
    assert L.rule("R2 Cramer: x*_2 det B_JJ, J=(0, 1, 2)").detail == r"$x^*_{2}\,\det B_{\{0,1,2\}}$, $J=\{0,1,2\}$"


def test_diagrams_render(tmp_path):
    a, B = ELEVEN["A-e competitor + 3-cycle (Thm 3.5)"]
    eco = Ecology(a, B, "Theorem 3.5")
    d = certificate_diagram(eco, complete(eco), label_style="both")
    d.save(str(tmp_path / "cert"), formats=("png",))
    assert (tmp_path / "cert.png").exists() and (tmp_path / "cert.dot").exists()
    dot = (tmp_path / "cert.dot").read_text()
    assert "det B" in dot and "p_" in dot and "⊥" in dot
    eco = Ecology([1, -1, 1], np.array([[1, -1, 1], [1, 1, -1], [-1, 1, 1]]))
    res = complete(eco, record_tree=True)
    search_tree_diagram(eco, res, use=Ecology.ALL_RULES).save(str(tmp_path / "tree"), formats=("pdf",))
    assert (tmp_path / "tree.pdf").exists()
