# Five-species screen: results

`scripts/screen_all.py --rules R1,R2,R3,R4,R5,R6` ran chirotope propagation + rules R1–R6 at the root (no branching)
on all 248,436 non-trivial five-species ecologies up to relabeling (count matches Meng et al. 2025; run of 2026-09,
about 100 CPU-minutes on two cores). Theorem numbers refer to Çelik–Isoldi–Portakal–Zucal (arXiv:2609.09286).

**21 ecologies are proved impossible at the root**; the remaining 248,415 are not decided by root
propagation. The 59 survivors with the fewest undetermined Plücker signs (≤110 of 252) all admit a
sign-consistent completion under branching, so no further kills come from the chirotope alone at that end.

Meng et al. report an upper bound of 30 candidate impossible five-species ecologies (Fig. S3). All 21 kills
here are theorems (the rules are sound). Rule R8, added later, kills a 22nd (`parity_k2`) and, with
one wrong growth sign admitted (rule R8, docs/SIGNRULES.md §3 and §7 (d)), two more; the extension screen
(`signrules-screen 5`) shows that every other five-species ecology
is possible: exactly 24 are impossible.

## The 21 proved impossible ecologies

Notation: `i±<-jk` means species i has growth-rate sign ± and benefits from species j and k
(b_ij < 0); `.` = benefits from nobody (pure competitor). All other off-diagonal entries are competitive (+).

### Symmetric (balanced) — 8, all covered by Theorems 3.1–3.2

| # | a | benefactors | structure | kill |
|---|---|---|---|---|
| 1 | [-1, -1, -1, -1, -1] | `0-<-1234; 1-<-0234; 2-<-0134; 3-<-0124; 4-<-0123` | obligate mutualism (Thm 3.1) | R4 Z-block (M-matrix) |
| 8 | [1, -1, -1, -1, -1] | `0+<-.; 1-<-234; 2-<-134; 3-<-124; 4-<-123` | 1+4 split, quadruple obligate, singleton a>0 (Thm 3.2) | R4 Z-block (not M-matrix) |
| 9 | [1, 1, -1, -1, -1] | `0+<-1; 1+<-0; 2-<-34; 3-<-24; 4-<-23` | 2+3 split, triple obligate, pair a=(+,+) (Thm 3.2) | R4 Z-block (not M-matrix) |
| 10 | [-1, 1, -1, -1, -1] | `0-<-1; 1+<-0; 2-<-34; 3-<-24; 4-<-23` | 2+3 split, triple obligate, pair a=(-,+) (Thm 3.2) | R4 Z-block (not M-matrix) |
| 11 | [-1, -1, 1, 1, 1] | `0-<-1; 1-<-0; 2+<-34; 3+<-24; 4+<-23` | 2+3 split, pair obligate, triple a=(+,+,+) (Thm 3.2) | R3 principal minor |
| 12 | [-1, -1, -1, 1, 1] | `0-<-1; 1-<-0; 2-<-34; 3+<-24; 4+<-23` | 2+3 split, pair obligate, triple a=(-,+,+) (Thm 3.2) | R3 principal minor |
| 13 | [-1, -1, -1, -1, 1] | `0-<-1; 1-<-0; 2-<-34; 3-<-24; 4+<-23` | 2+3 split, pair obligate, triple a=(-,-,+) (Thm 3.2) | R3 principal minor |
| 14 | [-1, -1, -1, -1, -1] | `0-<-1; 1-<-0; 2-<-34; 3-<-24; 4-<-23` | 2+3 split, both groups obligate (Thm 3.2) | R3 principal minor |
### Asymmetric — 13 (three covered by Theorems 3.4–3.6, ten new at n=5)

| # | a | benefactors | structure | kill |
|---|---|---|---|---|
| 2 | [-1, -1, 1, 1, 1] | `0-<-1; 1-<-0; 2+<-0134; 3+<-0124; 4+<-0123` | obligate pair + 3 universal beneficiaries a>0 (parity, k=1) | R4 Z-block (M-matrix) |
| 3 | [1, -1, -1, 1, 1] | `0+<-.; 1-<-2; 2-<-1; 3+<-0124; 4+<-0123` | competitor + obligate pair + 2 universal beneficiaries a>0 | R5 Laplace expansion of det B |
| 4 | [-1, -1, -1, 1, 1] | `0-<-2; 1-<-0; 2-<-1; 3+<-0124; 4+<-0123` | obligate 3-cycle + 2 universal beneficiaries a>0 (parity, k=1) | R4 Z-block (M-matrix) |
| 5 | [1, -1, -1, -1, 1] | `0+<-.; 1-<-0; 2-<-3; 3-<-2; 4+<-0123` | competitor, 1 obligate on 0, obligate pair, 1 universal beneficiary a>0 | R6 Schur block, monomial minor |
| 6 | [1, -1, -1, -1, 1] | `0+<-.; 1-<-3; 2-<-1; 3-<-2; 4+<-0123` | competitor + obligate 3-cycle + universal beneficiary a>0 | R5 Laplace expansion of det B |
| 7 | [-1, -1, -1, -1, 1] | `0-<-3; 1-<-0; 2-<-1; 3-<-2; 4+<-0123` | super-predator + obligate 4-cycle (Thm 3.6) | R4 Z-block (M-matrix) |
| 15 | [-1, -1, 1, -1, -1] | `0-<-2; 1-<-0; 2+<-1; 3-<-4; 4-<-3` | 3-cycle (one species a>0) + obligate pair | R1 coefficient c_1 |
| 16 | [-1, -1, -1, -1, 1] | `0-<-2; 1-<-0; 2-<-1; 3-<-4; 4+<-3` | obligate 3-cycle + pair (one a>0) | R1 coefficient c_1 |
| 17 | [-1, -1, -1, -1, -1] | `0-<-2; 1-<-0; 2-<-1; 3-<-4; 4-<-3` | obligate 3-cycle + obligate pair (k=2) | R1 coefficient c_1 |
| 18 | [1, -1, -1, -1, -1] | `0+<-.; 1-<-0; 2-<-1; 3-<-4; 4-<-3` | competitor, obligate chain 0->1->2, obligate pair | R5 Laplace expansion of det B |
| 19 | [1, -1, -1, -1, -1] | `0+<-.; 1-<-0; 2-<-4; 3-<-2; 4-<-3` | competitor, 1 obligate on 0, obligate 3-cycle {2,3,4} | R6 Schur block, monomial minor |
| 20 | [1, -1, -1, -1, -1] | `0+<-.; 1-<-4; 2-<-1; 3-<-2; 4-<-3` | competitor + obligate 4-cycle (Thm 3.5) | R5 Laplace expansion of det B |
| 21 | [-1, -1, -1, -1, -1] | `0-<-4; 1-<-0; 2-<-1; 3-<-2; 4-<-3` | obligate 5-cycle (Thm 3.4) | R4 Z-block (M-matrix) |
## Observations

* Every killed pattern has each species benefiting from **none, exactly one, or all** of the others (or,
  in the balanced cases, from exactly its own group). No pattern with an intermediate benefactor set is
  killed at the root — either those are all possible, or they need branching / new rules.
* Patterns #2, #4, #7, #21 are the cases k = 1 of the parity family (universal beneficiaries with a > 0 plus k
  disjoint obligate cycles, impossible for k odd via det(S B P_σ) = (−1)^k det B > 0). Pattern #17 (k = 2)
  is killed instead by coefficient positivity: every principal 4-minor is forced negative, so c_1 < 0.
* The competitor (a species benefiting from nobody, a > 0) can be attached to any of the cycle/pair
  structures: #3, #5, #6, #18, #19, #20 — the mechanism of Theorem 3.5 (partial-Z rows handled by Laplace
  expansion of det B or a Schur block) applies more widely than the cycle case.
* Direct equilibrium sampling (20k samples per pattern) finds no Hurwitz instance for any of the 21,
  consistent with the proofs.

## Files

* `killed_patterns.json` — the 21 patterns with a, B, benefactor description, structure label and the
  rule that produced the contradiction.
* Full per-ecology results (`screen_results.jsonl`, 248,436 lines) are not committed; rerun
  `scripts/screen_all.py --rules R1,R2,R3,R4,R5,R6` to regenerate them (its default `--out` is this folder;
  resumable with `--resume`, time-boxed with `--seconds`).
