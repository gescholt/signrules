# signrules

Sign completion on the ecological Grassmannian: a Python package that proves Lotka–Volterra sign
patterns **impossible**, with a machine certificate for every proof and an independent checker that
replays it.

A Lotka–Volterra community `ẋ = diag(x)(a − Bx)` is described qualitatively by the signs of the growth
rates `a` and of the interaction matrix `B`. The pattern is *possible* if some magnitudes with these signs
have a positive equilibrium that is locally asymptotically stable, and *impossible* otherwise. Encoding
`[diag(a) | B]` as a point of the Grassmannian `Gr(n, 2n)`, every Plücker coordinate is a minor of `B` up to
a known sign, so the feasibility, stability, M-matrix and Schur-complement arguments of the two papers below
all become constraints on one set of sign variables, the same variables on which the oriented-matroid
(Grassmann–Plücker) propagation acts. The package runs

* **chirotope propagation**, the 3-term Grassmann–Plücker relations (complete for uniform chirotopes), and
* the **sign rules R1–R8**: stability (`det B > 0`, positive coefficients of the characteristic polynomial),
  Cramer's rule on every subsystem, M-matrix minor signs of balanced patterns, Laplace expansions, the
  oriented-matroid Farkas lemma (tope / circuit orthogonality, complete for feasibility), and the real
  eigenvalue of benefactor-permutation patterns,

to a fixed point at every node of a depth-first sign-completion search. A pattern that is *killed* is a
proved impossible ecology, and the backward slice of the log from the contradiction is its certificate.
The rules are checked never to reject a genuine feasible, stable system.

It extends the sign-completion algorithm (the Maple implementation accompanying the paper, "the Maple
algorithm" below) of

* T. Ö. Çelik, P. A. Haas, G. Scholten, K. Wang, G. Zucal, *Strata of ecological coexistence via
  Grassmannians*, Ann. Comb. (2026), [doi:10.1007/s00026-026-00837-7](https://doi.org/10.1007/s00026-026-00837-7),
  [arXiv:2509.00165](https://arxiv.org/abs/2509.00165) ("the Annals paper"),

with the arguments of

* T. Ö. Çelik, V. A. Isoldi, I. Portakal, G. Zucal, *Stable coexistence in ecologies and games* (2026),
  [arXiv:2609.09286](https://arxiv.org/abs/2609.09286) ("CIPZ"),

read as propagation rules. The patterns both papers start from are those of

* Y. Meng, S. Horvát, C. D. Modes, P. A. Haas, *Impossible ecologies: interaction networks and stability
  of coexistence in ecological communities*, Cell Systems 16 (2025) 101297,
  [doi:10.1016/j.cels.2025.101297](https://doi.org/10.1016/j.cels.2025.101297).

The converse tool, certificates that a pattern *is* possible (tropical witnesses), lives in the companion
repository [`feasible_ecologies`](https://github.com/gescholt/feasible_ecologies), which will depend on
this package.

## Install

Not on PyPI; Python 3.10 or later.

```bash
git clone git@github.com:gescholt/signrules.git && cd signrules
pip install -e .                 # the package: numpy, matplotlib, networkx
pip install -e ".[all,dev]"      # + rich and questionary (signrules-live), plotly (HTML trees), pytest
pytest -q                        # about a minute
```

Graphviz (`dot`; `apt install graphviz` or `brew install graphviz`) is needed only for the Graphviz
drawings: certificates, completions and `write_search_tree`.

## Conventions

Species are 0-based. `b_ii > 0`. `B[i, j] < 0`: species `j` benefits species `i` (food, a mutualist);
`B[i, j] > 0`: `i` loses to `j` (competition, predation). `a_i > 0`: species `i` grows when alone (free);
`a_i < 0`: it does not (obligate). A pattern is *trivial* when some obligate species has only harmful
interactions (it cannot have a positive equilibrium at all). In the matrix `[diag(a) | B]` the columns
`0..n-1` carry `diag(a)` and `n..2n-1` carry `B`; the Plücker coordinate `p_0257` is the minor on columns
`{0, 2, 5, 7}` ([docs/SIGNRULES.md](docs/SIGNRULES.md) §1 gives the dictionary to the minors of `B`). A
sign table is a flat list indexed by basis id with `+1`, `-1`, `UNKNOWN = 0`; `eco._minor_key[(P, Q)] =
(S, c)` gives `sign(det B_PQ) = c · chi[S]`.

## Use

```python
from signrules import Ecology, complete, ELEVEN, write_certificate, replay_kill
from signrules.certificate import certificate_steps

# a free competitor over an obligate 3-cycle (1 eats 2 eats 3 eats 1): CIPZ Theorem 3.5
a, B = ELEVEN["A-e competitor + 3-cycle (Thm 3.5)"]
eco = Ecology(a, B)
result = complete(eco)
print(result.summary())                      # impossible at the root (R5 laplace ...)
result.root_log                              # every sign a rule forced at the root, with the rule and its premises
result.root_gp_log                           # the signs forced by the 3-term relations
for step in certificate_steps(eco, result)[0]:
    print(step)                              # (basis, sign, rule, premises): the backward slice of the log
print(replay_kill(eco, result))              # the independent checker re-derives every step
write_certificate(eco, result, "cert")       # cert.png + cert.pdf + cert.dot: the derivation DAG

# feasibility + weak stability only, the Maple algorithm: 64 weak completions of CIPZ Figure 2(d)
a, B = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]
eco = Ecology(a, B)
result = complete(eco, use=("R1", "R2n"), record_tree=True)
print(result.summary())                      # 64 completions, 131 nodes, 4 pruned

from signrules import search_tree_text, save_search_tree_png
print(search_tree_text(eco, result, max_depth=3))
save_search_tree_png(eco, result, "tree.png")
```

`complete(eco, use=...)` takes the rule names `R1` stability, `R2` subsystem Cramer (`R2n`: full system
only), `R3` balanced patterns, `R5` Laplace, `R7` tope (subsumes `R2`), `R8` real eigenvalue, and the
retired `R4` Z-blocks and `R6` Schur blocks, which the other rules subsume (`docs/SIGNRULES.md`, "R4 and R6
are subsumed"); `complete()` idles them by default (`idle=()` runs them). The Graphviz drawings
(`write_certificate`, `write_completion`, `write_search_tree`) accept `label_style="plucker" | "minor" |
"both"`; the matplotlib and plotly trees and `search_tree_text` take `plucker=False` to label by minors.

Patterns to start from are in `signrules.patterns`: `ELEVEN` (the eleven four-species patterns of Meng et
al., proved impossible in CIPZ), `parity_pattern` (the benefactor-permutation family of R8), `PATTERNS`
(five- and six-species presets), `random_pattern`, `canonical_classes(n)`. The data of the classifications is
in `signrules.data`: `root_kills(5)` returns the 24 minimal impossible five-species patterns,
`screen_records(6)` the 143 six-species candidates with their verdicts (the `open` rows are certified
possible only in `feasible_ecologies`).

### Exploring one pattern

```bash
signrules-live --preset parity_k2 --max-nodes 3000 --html parity_k2.html    # live terminal dashboard, then the tree viewer
signrules-live --a=-1,-1,-1,-1,1 --B=1,-1,1,1,1/-1,1,1,1,1/1,1,1,-1,1/1,1,-1,1,1/-1,-1,-1,-1,1 --max-nodes 3000 --html mine.html
signrules-live --random 5 --max-nodes 3000 --html random.html              # a random pattern (--seed S repeats it)
signrules-live                                                             # arrow-key menus for everything
signrules-tree --a=1,-1,-1 --B=1,1,1/-1,1,1/-1,-1,1 --rules R1,R2n --out tree.html   # only the HTML viewer
```

Write `--a=...` with an equals sign: a value starting with `-` is otherwise read as an option. The dashboard
(needs the `live` extra) shows the share of the tree already swept, the signs forced and the branches closed
by each rule, the current path, completions and contradictions by depth. The tree viewer is one
self-contained HTML file: click a node to see the signs it forced with their rules, what the known signs say
about stability and feasibility, and every minor of `B`; at a contradiction it names the violated rule and
its premises.

## Results and how to reproduce them

| case | result | command |
|---|---|---|
| **n = 4** | all eleven candidate impossible ecologies of Meng et al. are killed at the root; with `--rules R1,R2n` the search reproduces the completion counts of the Maple algorithm quoted in CIPZ (S-d, S-e, A-a, A-c: 64, 64, 1996, 1528) | `python scripts/impossible_eleven.py [--rules R1,R2n]` |
| **n = 5** | exhaustive root screen of all 248,436 non-trivial ecologies up to relabelling: 21 kills under R1–R6 (10 of them new), 24 with R8 | `python scripts/screen_all.py --rules R1,R2,R3,R4,R5,R6` (100 CPU-minutes); `signrules/data/n5/killed_patterns.json` |
| **n = 5, 6** | extension screens by the deletion lemma: **83** candidates at n = 5, **24** killed; **143** at n = 6, **41** killed; the rest are certified possible in `feasible_ecologies`, so these are exactly the minimal impossible patterns | `signrules-screen 5 --out out/n5.jsonl`, `signrules-screen 6 --out out/n6.jsonl`; the shipped lists are `signrules/data/n{5,6}/extension_screen.jsonl` |
| **all n** | benefactor-permutation families: every species obligate with a single benefactor is impossible for any benefactor permutation; with free species benefiting from everyone, impossible for any non-trivial permutation | `python scripts/families.py` |
| **certificates** | the 11 + 24 + 41 root kills replay in the independent checker; flipping any single step makes the replay fail | `for n in 4 5 6; do python scripts/replay_certificates.py $n --tamper; done` |
| **rule ladder** | what each rule buys on the minimal impossible patterns | `python scripts/rule_ladder.py 4 5 6 --out out/ladder.jsonl --tex out/ladder.tex` |
| **ablations** | R4 and R6 never fire under the other rules; R4 subsumed up to one branching, exhaustively per block configuration | `python scripts/rule_ablation.py 3 4 --census`, `python scripts/r4_subsumption.py 3 4 --rules R1,R7 --others` |

Without `--out`, `signrules-screen` and `rule_ladder.py` overwrite the shipped data files.
[docs/SIGNRULES.md](docs/SIGNRULES.md) is the full account: the Plücker/minor dictionary, each rule and the
theorem it encodes, certificates and their replay, the search, soundness checks, the screens and the families.

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
  data/               n5/, n6/: candidate lists and kills; n3/, n4/: sign-correlation census;
                      rule_ablation/, r4_subsumption/, rule_ladder.jsonl: records of the ablations
scripts/              the drivers of the table above
tests/                conventions, soundness against random feasible-stable systems, golden search trees, replay
docs/                 SIGNRULES.md, its figures (docs/figures, with their generators), rule_ladder_table.tex
.github/workflows/    pytest on Python 3.10 and 3.13
CITATION.cff
```

## Licence and citation

MIT. Please cite the Annals paper and CIPZ ([CITATION.cff](CITATION.cff)).
