"""The interactive tree viewer: per-node records and the generated HTML."""
import json
import re


from signrules import Ecology, complete
from signrules.treeview import tree_data, tree_viewer_html


def _eleven0():
    from signrules.patterns import ELEVEN
    name, (a, B) = next(iter(ELEVEN.items()))
    return Ecology(a, B, name)


def test_nodes_carry_their_forced_signs_and_conflicts():
    eco = _eleven0()
    res = complete(eco, use=("R1", "R2n"), record_tree=True)
    data = tree_data(eco, res, ("R1", "R2n"))
    assert len(data["nodes"]) == len(res.tree)
    # the forced counts agree with the per-node step lists for branch nodes
    for nd, rec in zip(res.tree, data["nodes"]):
        if nd["status"] == "branch":
            assert nd["forced"] == len(rec["steps"])
    conflicts = [r for r in data["nodes"] if r["status"] == "conflict"]
    assert conflicts and all(r["conflict"]["rule"] in ("GP", "R1", "R2") for r in conflicts)
    assert all(st["rule"] in ("GP", "R1", "R2") for r in data["nodes"] for st in r["steps"])


def test_html_embeds_valid_json_and_no_external_scripts():
    eco = Ecology([1, -1, -1], [[1, 1, 1], [-1, 1, 1], [-1, -1, 1]], "omnivory")
    res = complete(eco, use=("R1", "R2n"), record_tree=True)
    page = tree_viewer_html(eco, res, ("R1", "R2n"))
    assert "<script src" not in page
    payload = re.search(r"const D = (.*?);\nconst C", page, re.S).group(1)
    data = json.loads(payload.replace("<\\/", "</"))
    assert data["n"] == 3 and data["B"][1][0] == -1
    assert sum(1 for n in data["nodes"] if n["status"] == "completion") == len(res.completions)


def test_root_conflict_tree():
    eco = _eleven0()
    res = complete(eco, record_tree=True)          # all rules: killed at the root
    assert res.root_conflict
    data = tree_data(eco, res)
    assert len(data["nodes"]) == 1 and data["nodes"][0]["conflict"]["premises"]
