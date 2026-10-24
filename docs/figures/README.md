# Figures for docs/SIGNRULES.md

The Graphviz and matplotlib figures of the sign-rule notes, in `docs/figures/`.

Every figure is generated; no label is typed by hand. Regenerate all of them with

```bash
python docs/figures/make_figures.py     # dictionary, tope, zblock, balanced
python docs/figures/make_graphs.py      # certificate_*, search_tree_*, completion_sd
python docs/figures/make_diagrams.py    # algorithm, sign_rules
```

| file | shows | made by |
|---|---|---|
| `algorithm.*` | one run of `complete()`: initial signs → propagation to a fixed point → ⊥ or branch | `make_diagrams.py` |
| `sign_rules.*` | the rules grouped by the hypothesis they use, new ones highlighted, the retired R4 and R6 dotted | `make_diagrams.py` |
| `dictionary.*` | a Plücker coordinate $p_S$ of $[\operatorname{diag}(a)\,\|\,B]$ collapsing to a minor of $B$ | `make_figures.py` — `fig_dictionary(eco, S)` |
| `tope.*` | R7: the feasible sign vectors of $a$ as cells of a hyperplane arrangement ($n=2$) | `make_figures.py` — `fig_tope(B)` |
| `zblock.*` | R4 (retired, subsumed by the other rules): a row-sign twist and a column permutation exposing a Z-matrix | `make_figures.py` — `fig_zblock(eco)` |
| `balanced.*` | R3: a balanced signed graph and the signature making $SBS$ a Z-matrix | `make_figures.py` — `fig_balanced(eco)` |
| `certificate_thm35.*` | machine proof of Theorem 3.5 (competitor + obligate 3-cycle) | `make_graphs.py` — `write_certificate` |
| `certificate_lemma38.*` | machine proof of Figure 3(a) / Lemma 3.8 | `make_graphs.py` — `write_certificate` |
| `search_tree_sd.*` | the binary search tree of Figure 2(d) under `R1,R2n` (`engine="sweep", branching="first"`): 137 nodes, 64 completions, 10 contradictions | `make_graphs.py` — `write_search_tree` |
| `search_tree_n3.*` | the full 8-completion tree of an open $n=3$ pattern | `make_graphs.py` — `write_search_tree` |
| `completion_sd.*` | derivation of one weak completion of Figure 2(d) under `R1,R2n`: branch choices and every forced sign with its rule | `make_graphs.py` — `write_completion` |

Each graph comes as `.png`, `.pdf` and a `.dot` with plain-text labels (the drawing itself is laid
out by Graphviz and typeset by matplotlib, `signrules/gvrender.py`, so that labels can be
mathematics).

## Conventions

* **Language.** A sign is a Plücker coordinate $p_S$, $S\subset\{0,\dots,2n-1\}$, $|S|=n$, of the
  point $[\operatorname{diag}(a)\,|\,B]$ of $\mathrm{Gr}(n,2n)$; the dictionary figure shows how
  $p_S$ collapses to a minor of $B$. Minors of $B$ are named only inside the rules that are stated
  in terms of them (Laplace expansions, Cramer, Z-blocks). Indices are 0-based as in the code
  (`labels.ONE_BASED = True` prints them 1-based as in the papers).
* **Colour** carries three meanings. *Sign*: red = positive, blue = negative — for entries of
  $B$ (competition / mutualism in the interaction graphs), for coordinates in derivation graphs,
  and for the branch taken in a search tree. *Rule*: one hue per rule (`style.RULE_COLORS`, none
  of them the red or the blue of the signs) — the border and title of a rule's box in `sign_rules`, the headline of
  its box and the border of the $\bot$ it raises in a certificate, and the $\bot$ leaves it closes
  in a search tree, so the same rule looks the same in every figure. *Status* of a node that has
  no sign: a completion is green, an unexplored subtree grey, a contradiction black. Pale tints
  mark subsets — chosen columns of a matrix, feasible cells, the copies of a repeated subtree —
  and never a sign. Everything else is greyscale.
* **Provenance** of a sign is the border: thin grey = read off the entries, black = forced by a
  rule, dashed = a branch choice.
* **Search trees** use two more channels: the width of an edge grows with the number of signs the
  propagation forced after that choice (a thick edge is a choice that decided a lot), and
  subtrees that the search visits several times — same coordinates branched on, same choices,
  same propagation, same outcome — are boxed in the same pale tint with their multiplicity (in
  Figure 2(d) of Çelik–Isoldi–Portakal–Zucal the classical search, drawn with `engine="sweep",
  branching="first"`, does the same 15-node subtree eight times: the last three
  coordinates are independent of the earlier choices). `shade_open=True` additionally greys a node
  by the number of coordinates still open there.
* **Not shown**: node counts, numbers of open coordinates, "hypothesis" nodes, and prose captions
  — they belong in the text (or in the plotly hover), not on the figure.
