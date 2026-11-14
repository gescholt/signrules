# signrules

[![tests](https://github.com/gescholt/signrules/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/gescholt/signrules/actions/workflows/ci.yml)

A Python package that proves Lotka–Volterra sign patterns **impossible**, with a machine certificate for
every proof and an independent checker that replays it.

A Lotka–Volterra community `ẋ = diag(x)(a − Bx)` is described qualitatively by a *sign pattern*: the signs
of the growth rates `a` and of the interaction matrix `B`. The pattern is *possible* if some magnitudes
with these signs have a positive equilibrium that is locally asymptotically stable, and *impossible*
otherwise. `signrules` takes a pattern and searches for a proof that it is impossible.

<img src="docs/figures/network_A-e.png" width="200" alt="the interaction network of the pattern A-e">

*The pattern A-e. A dashed arrow i → j means that species i benefits from j; a thin red line joins two
species that each lose to the other. Pink species grow when alone, blue ones do not.*

For this pattern, `signrules run --preset A-e` prints:

```
╭─ verdict ────────────────────────────────────────────────────────────────────────────────────╮
│ IMPOSSIBLE  ·  killed at the root by R7 · tope                                               │
│ A-e competitor + 3-cycle (Thm 3.5) · n = 4, sign a = (+,−,−,−) · rules R1–R8                 │
│ No magnitudes with these signs have a stable positive equilibrium: the certificate below is  │
│ the proof, replayed by an independent checker.                                               │
│                                                                                              │
│   #  sign        minor                      rule                     uses                    │
│   1  p_4567 > 0  det B > 0                  R1 · stability           —                       │
│   2  p_0367 < 0  det B_{1,2},{2,3} > 0      R7 · tope                4 entries               │
│    …                                                                                         │
│  14  p_0256 > 0  det B_{1,3},{1,2} > 0      GP · 3-term relation     #12 #13 + 3 entries     │
│      ⊥           contradiction              R7 · tope                #14 + 4 entries         │
│                                                                                              │
│ replay OK  14 steps and the contradiction re-derived by the independent checker              │
╰──────────────────────────────────────────────────────────────────────────────────────────────╯
```


## Install

Python 3.10 or later; not on PyPI yet.

```bash
git clone https://github.com/gescholt/signrules.git && cd signrules
pip install -e ".[all]"          # the package, the app and the HTML tree viewer
pip install -e .                 # the core alone: numpy, matplotlib, networkx, rich
pip install -e ".[all,dev]"      # + pytest; then `pytest -q` (under a minute)
```

With [uv](https://docs.astral.sh/uv/), `uv sync --all-extras` creates the environment from the lock file
(`uv.lock`, the versions the tests last ran with) and `uv run pytest -q` runs the tests in it.

Graphviz (`dot`; `apt install graphviz` or `brew install graphviz`) is needed only for the drawings of
certificates, completions and `write_search_tree`.

## First run

```bash
signrules run --preset A-e                         # a four-species pattern proved impossible: verdict + certificate
signrules run --a +--- --B "+-++/-+++/+++-/++-+"   # a pattern of your own: signs of a, then the rows of B
signrules                                          # the app: pick a pattern in the library, press r
```

Signs are written `+-...`, `+,-,...` or `1,-1,...`; the rows of `B` are separated by `/` and its diagonal
is `+`. Each run ends in one of three verdicts:

* **Impossible**: no magnitudes with these signs have a stable positive equilibrium. The certificate
  printed below the verdict is the proof, and the independent checker has replayed it.
* **Not decided by the rules**: some sign tables survive every rule, so the pattern may be possible; the
  rules cannot say.
* **Undecided**: the search stopped at its limit, 20,000 nodes by default (`--max-nodes 0` lifts it;
  Ctrl-C stops early and keeps what was found).

## Use from Python

```python
from signrules import Ecology, complete

eco = Ecology.parse("+---", "+-++/-+++/+++-/++-+")   # the strings of the command line; lists of ±1 work too
print(eco)                                           # the sign matrix [diag(a) | B]
result = complete(eco)
print(result.verdict, result.summary())              # impossible  impossible at the root (R8 real eigenvalue ...)
```

For a kill, `result.root_log` holds every sign a rule forced at the root with the rule and its premises;
`certificate_steps(eco, result)` (from `signrules.certificate`) is the backward slice of that log, the
certificate; `replay_kill(eco, result)` re-derives every step with the independent checker;
`write_certificate(eco, result, "cert")` draws it as cert.png / .pdf / .dot (needs Graphviz). For an open
pattern, `result.completions` are the surviving sign tables and `eco.minor_signs(table)` reads one as
`{minor name: sign}`:

```python
from signrules import ELEVEN, search_tree_text, save_search_tree_png

a, B = ELEVEN["S-d {0,1}|{2,3}, {2,3} obligate, a01=(+,-)"]
eco = Ecology(a, B)
result = complete(eco, use=("R1", "R2n"), record_tree=True)   # the first algorithm: rules R1 and R2n only
print(result.verdict, result.summary())                       # open  64 completions, 131 nodes, 4 pruned
print(eco.minor_signs(result.completions[0]))                 # {'det B_{0},{1}': 1, ...}
print(search_tree_text(eco, result, max_depth=3))
save_search_tree_png(eco, result, "tree.png")
```

`complete(eco, use=..., max_nodes=..., max_completions=..., record_tree=..., progress=...)` are the
keywords a user sets; the docstring lists the others under "for the ablations". Patterns to start from
are in `signrules.patterns`: `ELEVEN` (the eleven four-species patterns of Meng et al., proved impossible
in CIPZ), `parity_pattern` (the benefactor-permutation family of R8), `PATTERNS` (five- and six-species
presets), `random_pattern`, `canonical_classes(n)`. The data of the classifications is in
`signrules.data`: `root_kills(5)` returns the 24 minimal impossible five-species patterns,
`screen_records(6)` the 143 six-species candidates with their verdicts (the `open` rows are certified
possible only in `feasible_ecologies`).

## The command line and the app

`signrules --help` and `signrules <command> --help` describe every option; `python -m signrules` is the
same command.

| command | what it does |
|---|---|
| `run` | complete one pattern with a live dashboard; prints the verdict and, for a kill, the certificate |
| `tree` | complete one pattern and write an interactive HTML viewer of the search tree |
| `replay` | kill one pattern at the root and check its certificate with the independent replayer |
| `edit` | the full-screen app (also plain `signrules`): library, matrix editor, network, rules, search, verdict |
| `presets` | list the named patterns (the eleven n = 4 patterns, the n = 5, 6 presets) with short ids |
| `screen` | list the candidates for impossible n-species patterns (deletion lemma) and run the rules on each |

```bash
signrules presets                                  # the named patterns, one line each (--long for the notes)
signrules run --preset tree_S-d                    # a search tree to watch: 64 completions, 4 contradictions
signrules run --preset S-d --rules annals          # the first algorithm (rules R1, R2n only): 64 completions
signrules run --random 5 --html random.html        # a random pattern, then the interactive tree viewer
signrules replay --preset A-d --cert out/          # root kill, independent replay, certificate as png/pdf/dot
```

A preset is named by its short id or any unambiguous prefix of its name. The `tree_*` presets are patterns
the root does not decide, so the search grows a tree with completions and contradictions; the other
four-species presets die at the root, and the five- and six-species ones grow trees too large to finish.
`--rules` takes `all` (the default), `annals` (R1 and R2n, the algorithm of the first paper below) or a
comma list of rule names. `--plain`, or output to a file or pipe, prints only the verdict.

During a search the terminal shows a live dashboard: the share of the tree already swept, the signs forced
and the branches closed by each rule, the current path, completions and contradictions by depth. The tree
viewer (`--html`, or `signrules tree`) is one self-contained HTML file: click a node to see the signs it
forced with their rules and every minor of `B`; at a contradiction it names the violated rule and its
premises.

Plain `signrules` opens a full-screen app on the pattern library (it needs the `all` extra). Its keys are
listed in its footer:

| key | what it does |
|---|---|
| ↑↓ enter (in the library) | preview a pattern's matrix and network, load it; the `tree_*` patterns come first |
| ←↑↓→, tab, space, click | move over the sign matrix `[diag(a) | B]`, jump between a and B, flip a sign |
| `<` `>` or `2`–`7` | one species less / more, or n directly |
| `u` | the rules (checkboxes; `a` all, `c` the two rules of the first algorithm, R1 and R2n) and the node limit |
| `r` | run the search: live dashboard, then the verdict (certificate and replay for a root kill); there `h` writes the tree viewer, `c` the certificate drawing, `x` stops, esc comes back |
| `e` | export the network as png, pdf and the `\econet` call of the papers |
| `s` | the classification: what the rules decide for n = 4, 5, 6 |
| `?`, `q` | help, quit |

## Conventions

Species are numbered from 0 (the papers number them from 1, so species 0 here is species 1 there). `b_ii > 0`. `B[i, j] < 0`: species `j` benefits species `i` (food, a mutualist);
`B[i, j] > 0`: `i` loses to `j` (competition, predation). `a_i > 0`: species `i` grows when alone (free);
`a_i < 0`: it does not (obligate). A pattern is *trivial* when some obligate species has only harmful
interactions (it cannot have a positive equilibrium at all). In the matrix `[diag(a) | B]` the columns
`0..n-1` carry `diag(a)` and `n..2n-1` carry `B`, which is how the app and the certificates number them: `p_0367` is the Plücker coordinate on columns 0, 3, 6, 7.
The dictionary between the minors of `B` and the Plücker coordinates, and the representation of a sign
table, are in [docs/SIGNRULES.md](docs/SIGNRULES.md) §1 and §3.

## How it works

Encoding `[diag(a) | B]` as a point of the Grassmannian `Gr(n, 2n)`, every Plücker coordinate is a minor of
`B` up to a known sign. The conditions for feasibility and stability then become constraints on one set of
sign variables, the same variables on which oriented-matroid propagation acts. The package runs

* **chirotope propagation**: the 3-term Grassmann–Plücker relations, and
* the **sign rules**: `R1` stability (`det B > 0`, positive coefficients of the characteristic polynomial),
  `R2` Cramer's rule on every subsystem (`R2n`: the full system only), `R3` M-matrix minor signs of
  balanced patterns, `R5` Laplace expansions, `R7` the oriented-matroid Farkas lemma (it subsumes `R2`),
  and `R8` the real eigenvalue of benefactor-permutation patterns,

to a fixed point at every node of a depth-first search over the open signs. A pattern that is *killed* is
proved impossible, and the backward slice of the log from the contradiction is its certificate. The tests
check that the rules never reject a genuine feasible, stable system. Two earlier rules, `R4` and `R6`, are
retired because the others subsume them; `complete()` leaves them idle by default.

<img src="docs/figures/algorithm.png" width="520" alt="one run of the search: initial signs, propagation to a fixed point, then a contradiction, a completion or a branch">

[docs/SIGNRULES.md](docs/SIGNRULES.md) is the full account: the dictionary between minors and Plücker
coordinates, each rule and the theorem it encodes, certificates and their replay, the search, the soundness
checks, the screens and the families. [docs/CHANGES.md](docs/CHANGES.md) lists what is new since the Maple
algorithm of the first paper below.

## Results and how to reproduce them

The scripts under `scripts/` are not needed to run the method: they reproduce the numbers of the papers.
"CIPZ" is the second paper below.

| case | result | command |
|---|---|---|
| **n = 4** | all eleven candidate impossible ecologies of Meng et al. are killed at the root; with `--rules annals` the search reproduces the completion counts of the first paper's algorithm quoted in CIPZ (S-d, S-e, A-a, A-c: 64, 64, 1996, 1528) | `python scripts/impossible_eleven.py [--rules R1,R2n]` |
| **n = 5** | exhaustive root screen of all 248,436 non-trivial ecologies up to relabelling: 21 kills by the rules before R8 (the command shown; R4 and R6 have since been retired as subsumed), 24 with R8 | `python scripts/screen_all.py --rules R1,R2,R3,R4,R5,R6` (100 CPU-minutes); `signrules/data/n5/killed_patterns.json` |
| **n = 5, 6** | extension screens by the deletion lemma: **83** candidates at n = 5, **24** killed; **143** at n = 6, **41** killed; the rest are certified possible in `feasible_ecologies`, so these are exactly the minimal impossible patterns | `signrules screen 5 --out out/n5.jsonl`, `signrules screen 6 --out out/n6.jsonl`; the shipped lists are `signrules/data/n{5,6}/extension_screen.jsonl` |
| **all n** | benefactor-permutation families: every species obligate with a single benefactor is impossible for any benefactor permutation; with free species benefiting from everyone, impossible for any non-trivial permutation | `python scripts/families.py` |
| **certificates** | the 11 + 24 + 41 root kills replay in the independent checker; flipping any single step makes the replay fail | `for n in 4 5 6; do python scripts/replay_certificates.py $n --tamper; done` |
| **rule ladder** | what each rule buys on the minimal impossible patterns | `python scripts/rule_ladder.py 4 5 6 --out out/ladder.jsonl --tex out/ladder.tex` |
| **ablations** | R4 and R6 never fire under the other rules; R4 subsumed up to one branching, exhaustively per block configuration | `python scripts/rule_ablation.py 3 4 --census`, `python scripts/r4_subsumption.py 3 4 --rules R1,R7 --others` |

Without `--out`, `signrules screen` and `rule_ladder.py` overwrite the shipped data files. The data is also
available from Python: `signrules.data.root_kills(5)` returns the 24 minimal impossible five-species
patterns, `screen_records(6)` the 143 six-species candidates with their verdicts.

## Papers

The package extends the sign-completion algorithm, implemented in Maple, of

* T. Ö. Çelik, P. A. Haas, G. Scholten, K. Wang, G. Zucal, *Strata of ecological coexistence via
  Grassmannians*, Ann. Comb. (2026), [doi:10.1007/s00026-026-00837-7](https://doi.org/10.1007/s00026-026-00837-7),
  [arXiv:2509.00165](https://arxiv.org/abs/2509.00165),

with the arguments of

* T. Ö. Çelik, V. A. Isoldi, I. Portakal, G. Zucal, *Stable coexistence in ecologies and games* (2026),
  [arXiv:2609.09286](https://arxiv.org/abs/2609.09286),

read as propagation rules. The patterns these papers start from are those of

* Y. Meng, S. Horvát, C. D. Modes, P. A. Haas, *Impossible ecologies: interaction networks and stability
  of coexistence in ecological communities*, Cell Systems 16 (2025) 101297,
  [doi:10.1016/j.cels.2025.101297](https://doi.org/10.1016/j.cels.2025.101297).

Certificates that a pattern the rules leave open *is* possible (tropical witnesses) are computed in the
companion repository [`feasible_ecologies`](https://github.com/gescholt/feasible_ecologies), which depends on
this package. Nothing here needs it.

## Layout

```
signrules/
  plucker_index.py, chirotope.py     minors of B <-> Plücker coordinates; 3-term Grassmann–Plücker propagation
  rules.py, engine.py, search.py     rules R1–R8 (Ecology, Conflict); trigger-driven fixpoint; complete()
  certificate.py, replay.py          derivation DAG of a kill or completion; independent replay of a certificate
  symmetry.py                        the pattern's relabelling stabilizer, used inside the search
  labels.py, style.py, gvrender.py, treeplot.py, treeview.py   the drawings and the HTML tree viewer
  cli.py, report.py, livesearch.py   the signrules command, the verdict, the live dashboard
  editor.py, network.py              the app (textual) and the interaction network as the \econet macro of the papers
  screen.py, patterns.py, data/      extension screen; ELEVEN, PATTERNS, parsers; candidate lists and kills
scripts/              reproduce the numbers of the table above
tests/                conventions, soundness against random feasible-stable systems, golden search trees, replay
docs/                 SIGNRULES.md, CHANGES.md, figures (with their generators), rule_ladder_table.tex
.github/workflows/    pytest on Python 3.10 and 3.13
```

## Licence, citation, contact

MIT. To cite the package, see [CITATION.cff](CITATION.cff). Questions and bug reports:
[open an issue](https://github.com/gescholt/signrules/issues).
