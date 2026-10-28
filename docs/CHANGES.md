# What changed since the sign-completion algorithm of our first paper

Internal note. Numbers are from [SIGNRULES.md](SIGNRULES.md), which describes how each was measured. A
one-page PDF: [changes.pdf](changes.pdf).

Terms used below:

- **Killed at the root:** the propagation reaches a contradiction before any branching.
- **Minimal impossible:** impossible, and every one-species deletion is impossible or trivial.

## 1. The starting point: the Maple algorithm of our first paper

- the 3-term Grassmann–Plücker relations;
- **R1**, weak stability: `det B > 0`, and some principal k-minor positive for every k;
- **R2n**, Cramer's rule on the full system `Bx = a`;
- branching on the open signs.

## 2. Methods added

**Rules from CIPZ's arguments:**

- **R2:** Cramer's rule on every subsystem `B_JJ x_J = w_J`, not only the full system.
- **R3:** for balanced patterns (a signature S makes S B S a Z-matrix), the principal-minor signs of an
  M-matrix.
- **R5:** Laplace expansion of every minor along a row or a column.
- **R4 (Z-blocks) and R6 (Schur blocks):** implemented, then measured to force no sign once the other rules
  run. They are idle by default.

**New rules:**

- **R7, tope rule:** `(+, sign a)` must be a tope of `[I | Bᵀ]` (the oriented-matroid Farkas lemma). It
  captures everything feasibility implies at the sign level, and contains R2.
- **R8, real eigenvalue:** suppose a signature S and a permutation σ ≠ id make S B P_σ a Z-matrix. If
  moreover S a > 0, or S a has exactly one wrong sign and σ an even number of cycles, then the Jacobian has
  a positive real eigenvalue, so the pattern is impossible.

**Search:**

- Failed-literal probing: try each value of an open sign; if one leads to a contradiction, the other is
  forced.
- Symmetry breaking with the relabellings that fix the pattern (SBDS).
- `decide` (branch on chosen coordinates only), `assume` (decide exactly whether some completion has given
  signs), `seed` (a random completion).

**Trust:**

- Every kill comes with a certificate: only the steps the contradiction depends on.
- An independent replayer re-derives each step from the identity or inequality it names. All 76 root-kill
  certificates (11 + 24 + 41) replay, and flipping any single step makes the replay fail.
- Soundness: on 75 random feasible, stable three- and four-species systems per test run, the rules never
  reject the true sign table.

**Screens:**

- Root propagation on all 248,436 non-trivial five-species patterns up to relabelling.
- Extension screens by the deletion lemma: 83 candidates at n = 5, 143 at n = 6.
- Benefactor-permutation patterns of every size n: σ ≠ id; the species fixed by σ grow alone and benefit
  from everybody; every other species i is obligate and benefits only from σ(i). All are impossible.

- **The Python API for a first run:** `Ecology.parse("+---", "+-++/-+++/+++-/++-+")` reads the sign
  strings of the command line, `print(eco)` shows the sign matrix, `result.verdict` is `impossible`,
  `open` or `undecided`, `eco.minor_signs(table)` names the minors of a completion; `complete()` documents
  its six user keywords apart from the ablation knobs; the unused `seed` keyword is gone.

## 3. Optimisations

| change | measured effect |
|---|---|
| **Trigger-driven engine:** after a sign is set, only the checks that read that coordinate run again | 5–9× faster than rerunning everything ("sweep") on the searches of Meng et al.'s eleven |
| **Flat sign tables:** checks compiled to coordinate numbers once per pattern | the 226 first-completion searches of the n = 5, 6 extension screens: 97 s → 51 s; n = 7: 10 → 5 ms per node |
| **Most-constrained branching:** the open sign read by the most checks | about 7× fewer pruned nodes than first-open, on the 32 kills at n = 4, 5, under `R1,R2n,R5,R6` |
| **Conflict-activity branching** (as VSIDS), the default | n = 6 screen (143 searches, node limit 20,000): 39,062 → 16,034 nodes; slowest candidate 290 s → 3.4 s |
| **R4 and R6 idle**, still registered so the trees do not change | 2,000-node searches on four random seven-species patterns: about 2× faster (2.4 → 1.1 s, …) |
| **Probing at every node** (optional) | slowest n = 6 candidate: 20,000 → 117 nodes, at about 2 × (open signs) propagations per node |
| **Symmetry breaking** (optional) | 24 complete searches of symmetric four-species patterns: 38,542 → 18,813 nodes |

No optimisation changes a verdict or the set of completions; only the trees and the time change.

## 4. What it decides now

Root kills of the minimal impossible patterns (n = 4: Meng et al.'s eleven; n = 5, 6: the extension
screens):

| rules | n = 4 (11) | n = 5 (24) | n = 6 (41) |
|---|---|---|---|
| R1, R2n (our first paper) | 0 | 0 | 0 |
| R1, R2 | 3 | 2 | 2 |
| R1, R7 | 8 | 13 | 21 |
| all but R8 | 11 | 21 | 34 |
| all rules | **11** | **24** | **41** |

At n = 5 there are 24 minimal impossible patterns. Ten of them are covered by no CIPZ theorem. Meng et al.
bounded the number of impossible five-species patterns by 30 (their Fig. S3); the true number is 24. At
n = 6 there are 41. Every other candidate has a tropical witness of possibility, so both lists are complete.

## 5. Achievements

| area | achievement | result |
|---|---|---|
| **n = 4** | Meng et al.'s eleven | all 11 killed at the root, in milliseconds each |
| **n = 5** | root propagation on all 248,436 non-trivial patterns, then the extension screen (83 candidates) | exactly 24 minimal impossible patterns, 10 beyond the CIPZ theorems; Meng et al.'s bound of 30 replaced by the exact count 24 |
| **n = 6** | extension screen from the 24 five-species patterns | 143 candidates: 41 minimal impossible, 102 possible |
| **all n** | benefactor-permutation patterns (σ ≠ id) | impossible for every n (rule R8) |
| **completeness** | rules plus tropical witnesses | every pattern on n ≤ 6 species decided |
| **rules** | beyond R1, R2n | R2, R3, R5 (CIPZ) and R7, R8 (new); R4, R6 force nothing beyond them |
| **trust** | certificates | all 76 root kills replay in the independent checker; any flipped step is caught |
| **soundness** | tests | the rules never reject a random feasible, stable system |
| **speed** | engine, data layout, branching | 5–9× (engine); screens 97 s → 51 s; n = 6 screen 39,062 → 16,034 nodes, slowest candidate 290 s → 3.4 s |
|  | probing, symmetry (optional) | slowest n = 6 candidate 20,000 → 117 nodes; symmetric four-species searches 38,542 → 18,813 nodes |
