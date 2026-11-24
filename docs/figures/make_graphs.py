"""
Derivation graphs and search trees for docs/SIGNRULES.md.

    python docs/figures/make_graphs.py         # writes docs/figures/certificate_*, search_tree_*, completion_*

The examples are the eleven four-species patterns of Meng et al. (2025), proved impossible in
Çelik–Isoldi–Portakal–Zucal (``signrules.patterns.ELEVEN``),
and one open three-species pattern; every label on the drawings is generated from the run
(signrules/certificate.py, labels.py).
"""

from pathlib import Path


from signrules.patterns import ELEVEN
from signrules import Ecology, complete, write_certificate, write_search_tree, write_completion

OUT = Path(__file__).resolve().parent
CLASSICAL = ("R1", "R2n")


def pattern(key_start, name):
    key = next(k for k in ELEVEN if k.startswith(key_start))
    a, B = ELEVEN[key]
    return Ecology(a, B, name)


# the figures were drawn with the sweep engine and first-open branching (the certificates follow the
# paper's proofs rule by rule, the trees branch in table order); the default engine finds shorter kills.
# R4 and R6 are idle (subsumed): the Z-block step of the paper's proofs appears as
# three circuits (R7) and Cramer's rule (R2) in certificate_thm35.
ORIGINAL = dict(engine="sweep", branching="first")


def main():
    # certificates: kills at the root under all rules
    for key, name, stem in [("A-e", "competitor + obligate 3-cycle (Theorem 3.5)", "certificate_thm35")]:
        eco = pattern(key, name)
        res = complete(eco, **ORIGINAL)
        assert res.root_conflict, name
        write_certificate(eco, res, str(OUT / stem))
        print("wrote", stem)

    # the classical algorithm (feasibility + weak stability) on Figure 2(d): a full search tree
    eco = pattern("S-d", "Figure 2(d)")
    res = complete(eco, use=CLASSICAL, record_tree=True, **ORIGINAL)
    write_search_tree(eco, res, str(OUT / "search_tree_sd"), use=CLASSICAL, completion_numbers=False)
    write_completion(eco, res, str(OUT / "completion_sd"), index=0)
    print("wrote search_tree_sd, completion_sd:", res.summary())


if __name__ == "__main__":
    main()
