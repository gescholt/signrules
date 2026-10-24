# signrules

Sign completion on the ecological Grassmannian: a Python package that proves Lotka–Volterra sign
patterns **impossible**, with a machine certificate for every proof and an independent checker that
replays it.

A Lotka–Volterra community `ẋ = diag(x)(a − Bx)` is described qualitatively by the signs of the growth
rates `a` and of the interaction matrix `B` (`b_ii > 0`; `B[i, j] < 0` means species `j` benefits
species `i`). The pattern is *possible* if some magnitudes with these signs have a positive equilibrium
that is locally asymptotically stable, and *impossible* otherwise. Encoding `[diag(a) | B]` as a point of
the Grassmannian `Gr(n, 2n)`, every Plücker coordinate is a minor of `B` up to a known sign, so the
feasibility, stability, M-matrix and Schur-complement arguments of the two papers below all become
constraints on one set of sign variables, the same variables on which the oriented-matroid
(Grassmann–Plücker) propagation acts. The package runs

* **chirotope propagation**, the 3-term Grassmann–Plücker relations (complete for uniform chirotopes), and
* the **sign rules R1–R8**: stability (`det B > 0`, positive coefficients of the characteristic polynomial),
  Cramer's rule on every subsystem, M-matrix minor signs of balanced patterns, Laplace expansions, the
  oriented-matroid Farkas lemma (tope / circuit orthogonality, complete for feasibility), and the real
  eigenvalue of benefactor-permutation patterns,

to a fixed point at every node of a depth-first sign-completion search. A pattern that is *killed* is a
proved impossible ecology, and the backward slice of the log from the contradiction is its certificate.
The rules are checked never to reject a genuine feasible, stable system.

It extends the sign-completion algorithm of

* T. Ö. Çelik, P. A. Haas, G. Scholten, K. Wang, G. Zucal, *Strata of ecological coexistence via
  Grassmannians*, Ann. Comb. (2026), doi:10.1007/s00026-026-00837-7, arXiv:2509.00165,

with the arguments of

* T. Ö. Çelik, V. A. Isoldi, I. Portakal, G. Zucal, *Stable coexistence in ecologies and games*,
  arXiv:2609.09286 (2026),

read as propagation rules. The converse tool, certificates that a pattern *is* possible (tropical
witnesses), lives in the companion repository
[`feasible_ecologies`](https://github.com/gescholt/feasible_ecologies), which depends on this package.

## Install

```bash
pip install -e ".[all,dev]"        # numpy, matplotlib, networkx; rich + questionary (live dashboard), plotly (HTML trees); pytest
pytest -q                          # about a minute
```

Graphviz (`dot`) is needed only for the certificate and completion drawings.

## Use

```python
from signrules import Ecology, complete, write_certificate

# obligate cyclic predation 0 -> 1 -> 2 -> 3 -> 0, competition elsewhere
B = [[ 1, -1,  1,  1],
     [ 1,  1, -1,  1],
     [ 1,  1,  1, -1],
     [-1,  1,  1,  1]]
eco = Ecology(a_signs=[-1, -1, -1, -1], B_signs=B)
result = complete(eco)
print(result.summary())          # impossible at the root (...)
result.root_log                  # every sign forced at the root, with the rule that forced it
write_certificate(eco, result, "cert")      # cert.png + cert.pdf + cert.dot: the derivation DAG of the kill

# feasibility + weak stability only: the algorithm of the Grassmannian paper
print(len(complete(eco, use=("R1", "R2n")).completions))

# a search tree, recorded and drawn
result = complete(eco, use=("R1", "R2n"), record_tree=True)
from signrules import search_tree_text, save_search_tree_png
print(search_tree_text(eco, result, max_depth=3))
save_search_tree_png(eco, result, "tree.png")
```

`complete(eco, use=...)` takes the rule names `R1` stability, `R2` subsystem Cramer (`R2n`: full system
only), `R3` balanced patterns, `R5` Laplace, `R7` tope (subsumes `R2`), `R8` real eigenvalue, and the
retired `R4` Z-blocks and `R6` Schur blocks, which the other rules subsume (`docs/SIGNRULES.md`, "R4 and R6
are subsumed"); `complete()` idles them by default. Every drawing accepts
`label_style="plucker" | "minor" | "both"`.

Patterns to start from are in `signrules.patterns`: `ELEVEN` (the eleven four-species candidates of Meng
et al. 2025), `parity_pattern` (the benefactor-permutation family of R8), `PATTERNS` (five- and
six-species presets), `random_pattern`, `canonical_classes(n)`.

### Exploring one pattern

```bash
signrules-live --preset parity_k2 --max-nodes 3000 --html parity_k2.html    # live terminal dashboard, then the tree viewer
signrules-live --a -1,-1,-1,-1,1 --B "1,-1,1,1,1/-1,1,1,1,1/1,1,1,-1,1/1,1,-1,1,1/-1,-1,-1,-1,1" --max-nodes 3000 --html mine.html
signrules-live --random 5 --max-nodes 3000 --html random.html              # a random pattern (--seed S repeats it)
signrules-live                                                             # arrow-key menus for everything
signrules-tree --a 1,-1,-1 --B "1,1,1/-1,1,1/-1,-1,1" --out tree.html      # only the HTML viewer
```

The dashboard shows the share of the tree already swept, the signs forced and the branches closed by each
rule, the current path, completions and contradictions by depth. The tree viewer is one self-contained
HTML file: click a node to see the signs it forced with their rules, what the known signs say about
stability and feasibility, and every minor of `B`; at a contradiction it names the violated rule and its
premises.

## Results and how to reproduce them

| | | command |
|---|---|---|
| **n = 4** | all eleven candidate impossible ecologies of Meng et al. are killed at the root; with `--rules R1,R2n` the search reproduces the completion counts of the Maple algorithm (64, 64, 1996, 1528) | `python scripts/impossible_eleven.py [--rules R1,R2n]` |
| **n = 5** | exhaustive root screen of all 248,436 non-trivial ecologies up to relabelling: 21 kills, 10 of them new; with R8, 24 | `python scripts/screen_all.py` (100 CPU-minutes); `signrules/data/n5/killed_patterns.json` |
| **n = 5, 6** | extension screens by the deletion lemma: **83** candidates at n = 5, **24** killed; **143** at n = 6, **41** killed; the rest are certified possible in `feasible_ecologies`, so these are exactly the minimal impossible patterns | `signrules-screen 5`, `signrules-screen 6 --base known`; `signrules/data/n{5,6}/extension_screen.jsonl` |
| **all n** | benefactor-permutation families: every species obligate with a single benefactor is impossible for any benefactor permutation; with free species benefiting from everyone, impossible for any non-trivial permutation | `python scripts/families.py` |
| **certificates** | the 11 + 24 + 41 root kills replay in the independent checker; flipping any single step makes the replay fail | `python scripts/replay_certificates.py {4,5,6} --tamper` |
| **rule ladder** | what each rule buys on the minimal impossible patterns | `python scripts/rule_ladder.py 4 5 6` |
| **ablations** | R4 and R6 never fire under the other rules; R4 subsumed up to one branching, exhaustively per block configuration | `python scripts/rule_ablation.py`, `python scripts/r4_subsumption.py` |

`docs/SIGNRULES.md` is the full account: the Plücker/minor dictionary, each rule and the theorem it
encodes, certificates, the search, soundness checks, the screens and the families.

## Layout

```
signrules/
  plucker_index.py    minors of B  <->  Plücker coordinates of [diag(a) | B]; the initial sign table
  chirotope.py        3-term Grassmann–Plücker propagation
  rules.py            rules R1–R8 (Ecology, Conflict); every forced sign is logged with its rule
  engine.py           trigger-driven fixpoint: one propagator per relation, re-run only when a coordinate changes
  search.py           complete(): propagation + rules at every node, conflict-driven branching order, probing
  certificate.py      derivation DAG of a kill or a completion; the search tree as text or Graphviz
  replay.py           independent replay of a certificate: every step re-derived from the identity it names
  symmetry.py         the pattern's relabelling stabilizer, used inside the search
  labels.py, style.py, gvrender.py, treeplot.py   the drawings (one vocabulary, one palette)
  treeview.py         self-contained HTML viewer of a search tree          (signrules-tree)
  livesearch.py       live terminal dashboard of one search               (signrules-live)
  screen.py           extension screen: candidates by the deletion lemma  (signrules-screen)
  patterns.py         ELEVEN, parity_pattern, PATTERNS, random_pattern, canonical_classes, encode/decode
  data/               candidate lists and kills (n = 5, 6), screen summary, ablation records
scripts/              the drivers of the table above
tests/                conventions, soundness against random feasible-stable systems, golden search trees, replay
docs/                 SIGNRULES.md and its figures
```

## Conventions

Species are 0-based. `b_ii > 0`. `B[i, j] < 0`: species `j` benefits species `i` (food, a mutualist);
`B[i, j] > 0`: `i` loses to `j` (competition, predation). `a_i > 0`: species `i` grows when alone (free);
`a_i < 0`: it does not (obligate). A sign table is a flat list indexed by basis id with `+1`, `-1`,
`UNKNOWN = 0`; `eco._minor_key[(P, Q)] = (S, c)` gives `sign(det B_PQ) = c · chi[S]`.

## Licence and citation

MIT. Please cite the two papers above (`CITATION.cff`).
