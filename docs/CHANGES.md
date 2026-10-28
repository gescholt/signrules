# What changed since the sign-completion algorithm of our first paper

Internal note. `signrules` starts from the sign-completion algorithm of our first paper (its Maple
implementation, "the Maple algorithm") and adds the arguments of the stable-coexistence paper and new ones
as propagation rules on the same sign variables. Every number below is from [SIGNRULES.md](SIGNRULES.md),
where the measurements are described. A one-page PDF of this note: [changes.pdf](changes.pdf).

## 1. The starting point

The Maple algorithm completes the signs of the Plücker coordinates of `[diag(a) | B]` with the 3-term
Grassmann–Plücker relations, **R1** (weak stability: `det B > 0` and a positive principal k-minor for every
k) and **R2n** (Cramer's rule on the full system), branching on the open signs. In the package these are the
rules `R1,R2n`:

- On the eleven four-species candidates it kills S-a, S-b and S-c only after branching,
  and leaves the other eight with 64 to 15,512 completions.
- On the ten new five-species kills it leaves 158–173 signs open at the root and never closes the search.

## 2. Methods added

- **Rules.**
  - R2: Cramer on every subsystem.
  - R3: M-matrix minors of balanced patterns.
  - R5: Laplace expansions of every minor.
  - R7: the tope rule, oriented-matroid Farkas (complete for feasibility; subsumes R2).
  - R8: the real eigenvalue of benefactor-permutation patterns.
  - R4 (Z-blocks) and R6 (Schur blocks), the arguments of the stable-coexistence paper, were implemented and then shown to be subsumed
    by the others (they never force a sign under them), so they are idle by default.
- **Search.**
  - Failed-literal probing.
  - Symmetry breaking with the pattern's stabilizer (SBDS).
  - Coarse cells (`decide`), exact sign correlations (`assume`), random completions (`seed`).
- **Trust.**
  - Every kill comes with a certificate: the backward slice of the log, from the entries to the
    contradiction.
  - An independent replayer re-derives each step from the identity or inequality it names. Flipping any
    single step makes the replay fail (tamper tests on all 11 + 24 + 41 kills).
  - Soundness tests: the rules never reject a random feasible, stable system.
- **Screens.**
  - The exhaustive root screen of all 248,436 non-trivial five-species classes.
  - Extension screens by the deletion lemma at n = 5 and 6.
  - The infinite families of benefactor-permutation patterns, impossible for every n.

## 3. Optimisations (measured in SIGNRULES.md §3)

| change | effect |
|---|---|
| trigger-driven engine: small monotone propagators, re-run only when a coordinate they read changes | 5–9× faster than the full sweep on n = 4 searches, same fixed points |
| flat sign tables, propagators compiled to ids once per pattern | n = 5, 6 extension screens 97 s → 51 s; n = 7 search 10 → 5 ms per node |
| most-constrained branching (the open sign read by the most propagators) | about 7× fewer pruned nodes than first-open on the 32 known n = 4, 5 kills |
| conflict-activity branching (as VSIDS), the default | n = 6 screen 39,062 → 16,034 nodes; slowest candidate 290 s → 3.4 s |
| R4 and R6 idle (subsumed), propagators kept for the branching order | n = 7 searches about 2× faster (2.4 → 1.1 s, …), identical trees |
| probing at every node (optional) | hardest n = 6 candidate 20,000 → 117 nodes, at a higher cost per node |
| symmetry breaking (optional) | nodes halve on exhaustive searches (38,542 → 18,813) |

Verdicts and completion sets are unchanged by every optimisation; only the trees and the time change.

## 4. What it decides now

Root kills of the minimal impossible patterns (n = 4: the eleven; n = 5, 6: the extension screens):

| rules | n = 4 (11) | n = 5 (24) | n = 6 (41) |
|---|---|---|---|
| R1, R2n (Maple algorithm) | 0 | 0 | 0 |
| R1, R2 | 3 | 2 | 2 |
| R1, R7 | 8 | 13 | 21 |
| all but R8 | 11 | 21 | 34 |
| all rules | **11** | **24** | **41** |

At n = 5 these are 24 minimal impossible patterns. Ten of them are covered by no theorem of the
stable-coexistence paper, and they sharpen the earlier upper bound of 30 candidates to 24. At n = 6 there are 41. Every other candidate has a
certificate of possibility (tropical witness) in the companion repository `feasible_ecologies`, so the
lists are complete.

## 5. Software

- Python (numpy; rich and textual for the terminal) instead of Maple; MIT licence; 237 tests.
- One command, `signrules`, opens a full-screen app:
  - a pattern library;
  - an editor for the sign matrix, with the interaction network drawn in the article's style;
  - a choice of rules and a node limit;
  - a live dashboard of the search;
  - the verdict, with the certificate and its replay.
- Subcommands (`run`, `tree`, `replay`, `presets`, `screen`) do the same for scripts.
- Drawings of certificates, search trees and networks, an HTML search-tree viewer, and the data of every
  screen ship with the package.
