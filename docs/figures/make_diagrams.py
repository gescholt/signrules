"""
The two overview diagrams of docs/SIGNRULES.md, drawn with the package's own renderer so that
their labels are mathematics:

  algorithm.*    one run of complete(): initial signs → propagation to a fixed point → ⊥ or branch
  sign_rules.*   the rules grouped by the hypothesis they use, feeding the sign table

    python docs/figures/make_diagrams.py
"""

from pathlib import Path

from signrules import Diagram, style

OUT = Path(__file__).resolve().parent
BLUE, GREEN, AMBER = "#e3eaf8", "#e6f3df", style.HILITE_A
OLD, NEW, RETIRED = "#e6e6e6", "#fff0c8", "#f7f7f7"


def algorithm():
    d = Diagram(nodesep=0.35, ranksep=0.45, fontsize=11)
    small = lambda t: (t, 9, style.GREY)
    d.node("input", [r"sign pattern  $(\mathrm{sign}\,a,\ \mathrm{sign}\,B)$"], fill=BLUE, edge="#9db0d6")
    d.node("dict", [r"minors of $B$  =  Plücker coordinates $p_S$ of $[\,\mathrm{diag}(a) \mid B\,] \in \mathrm{Gr}(n, 2n)$"],
           fill=BLUE, edge="#9db0d6")
    d.node("init", ["initial signs", small("the $p_S$ determined by the entries")], fill="#f4f4f4")
    d.cluster("fix", label="propagate to a fixed point", color="#aaaaaa")
    d.node("gp", ["Grassmann–Plücker", small("3-term relations")], fill=GREEN, edge="#8fbf7a", cluster="fix")
    d.node("rules", ["sign rules",
                     small("R1 stability · R2 Cramer · R3 balanced · R5 Laplace"),
                     small("R7 tope · R8 real eigenvalue   (R4, R6: retired, subsumed)")], fill=GREEN, edge="#8fbf7a", cluster="fix")
    d.node("conflict", ["contradiction?"], fill="#f4f4f4", shape="round")
    d.node("bottom", [(r"$\bot$", 14, "white"), ("impossible ecology", 10, "white"), ("certificate: the forced signs and their rules", 8.5, "#dddddd")],
           fill=style.BOTTOM, edge=style.BOTTOM)
    d.node("done", ["all $p_S$ known?"], fill="#f4f4f4")
    d.node("completion", ["completion", small("a consistent sign table, not a proof of possibility")], fill=AMBER, edge="#d8b96a")
    d.node("branch", ["branch on one open $p_S$:  $+$ / $-$"], fill="#f4f4f4")
    d.edge("input", "dict")
    d.edge("dict", "init")
    d.edge("init", "gp")
    d.edge("gp", "rules")
    d.edge("rules", "gp", constraint=False)
    d.edge("rules", "conflict")
    d.edge("conflict", "bottom", label="yes")
    d.edge("conflict", "done", label="no")
    d.edge("done", "completion", label="yes")
    d.edge("done", "branch", label="no")
    d.edge("branch", "gp", ls="--", constraint=False)
    d.save(str(OUT / "algorithm"))
    print("wrote algorithm")


def sign_rules():
    d = Diagram(rankdir="LR", nodesep=0.25, ranksep=0.9, fontsize=11)
    small = lambda t: (t, 9, style.GREY)

    def rule(nid, cluster, head, detail, new=True, retired=False):
        d.node(nid, [(head, 11, style.rule_color(nid)), small(detail)], fill=RETIRED if retired else (NEW if new else OLD),
               edge=style.rule_color(nid), lw=1.3, cluster=cluster, **({"ls": ":"} if retired else {}))

    d.cluster("id", label="determinant identities", color="#bbbbbb")
    rule("GP", "id", "Grassmann–Plücker", "3-term relations between the $p_S$", new=False)
    rule("R5", "id", "R5 · Laplace", "expand a minor along a row or column")
    rule("R6", "id", "R6 · Schur blocks", r"minors of $B / B_{PP}$ from their entries   (retired: GP + R5 + R7)", retired=True)
    d.cluster("feas", label=r"feasibility   $x^* > 0$", color="#bbbbbb")
    rule("R7", "feas", "R7 · tope", r"$(+, \mathrm{sign}\,a)$ is a tope of $[\,I \mid B^{T}\,]$   (oriented-matroid Farkas)")
    rule("R2", "feas", "R2 · Cramer on subsystems", r"$x^*_j \det B_{JJ}$ is a signed sum of minors   ($\subset$ R7)")
    d.cluster("stab", label=r"stability   $-\mathrm{diag}(x^*)\,B$ Hurwitz", color="#bbbbbb")
    rule("R1", "stab", "R1 · stability", r"$\det B > 0$;  some principal $k$-minor $> 0$ for each $k$", new=False)
    d.cluster("z", label=r"feasibility + stability + a Z-twist  $\Rightarrow$  M-matrix", color="#bbbbbb")
    rule("R3", "z", "R3 · balanced patterns", r"$S B S$ is a Z-matrix: principal minors $> 0$, almost-principal $\sim b_{ij}$")
    rule("R4", "z", "R4 · Z-blocks", r"a Schur block twisted to a Z-matrix: M-matrix minor signs, or a leading minor $\leq 0$   (retired)", retired=True)
    d.cluster("eig", label=r"feasibility + stability + a Z-twist of the whole matrix", color="#bbbbbb")
    rule("R8", "eig", "R8 · real eigenvalue", r"$S B P_\sigma$ a Z-matrix, $\sigma \neq \mathrm{id}$: $-\mathrm{diag}(x^*) B$ has a positive real eigenvalue   $\Rightarrow \bot$")
    d.node("table", ["sign table", small(r"signs of all minors of $B$ = all $p_S$")], fill=GREEN, edge="#8fbf7a", pad=(12, 8))
    d.node("out", [(r"$\bot$", 14, "white"), ("impossible ecology + certificate", 10, "white")], fill=style.BOTTOM, edge=style.BOTTOM)
    for nid, cl in (("R5", "id"), ("R7", "feas"), ("R1", "stab"), ("R3", "z")):
        d.edge(nid, "table", ltail=cl)
    d.edge("table", "out")
    d.edge("R8", "out", ltail="eig")
    d.legend_item("in the Annals paper (Çelik–Haas–Scholten–Wang–Zucal)", fill=OLD, edge="#999999")
    d.legend_item("new", fill=NEW, edge="#c99a2e")
    d.legend_item("retired: subsumed by the other rules, idle by default", fill=RETIRED, edge="#999999", ls=":")
    d.legend_item("border and title: the rule's colour, used for its ⊥ in every other figure", fill="white", edge=style.rule_color("R4"), lw=1.3)
    d.save(str(OUT / "sign_rules"))
    print("wrote sign_rules")


if __name__ == "__main__":
    algorithm()
    sign_rules()
