# Changelog

Every release is archived on Zenodo under [doi:10.5281/zenodo.24214426](https://doi.org/10.5281/zenodo.24214426).

## 0.5.3 (2026-11-24): first release on PyPI

The package is `ecological_sign_rules` on PyPI; it still imports as `signrules` and installs the `signrules`
command.

- **The classification checks from the package alone.** `signrules verify` checks four, five and six
  species from the shipped data: the four-species census (3044 classes: 715 killed at the root, 2329
  certified possible), the five- and six-species candidate lists rebuilt by the deletion lemma, every kill
  redone and its certificate replayed, and every certificate of possibility checked in exact arithmetic.
- New: `signrules.possible` (exact check of a certificate of possibility, no new dependency) and the data
  `census(4)` and `possible(5)`, `possible(6)` in `signrules.data`.
- The verdict lists the rules that run (R4 and R6, subsumed, are idle).
- Removed: the commands `signrules-live`, `signrules-tree`, `signrules-screen` (use `signrules run`, `tree`,
  `screen`); `signrules.data.killed_patterns()`; scripts that were exploration tools (`completions.py`,
  `rule_usage.py`, `search_trees_n4.py`, `search_trees_n5.py`, `sign_correlations_exact.py`) or duplicated a
  command (`extension_screen.py`).
- The records of the experiments (rule ablations, R4 subsumption, rule ladder, the exhaustive five-species
  screen) moved from the package to the `data/` folder of the repository; the wheel shrinks from 5.3 MB to
  1 MB. The source archive holds the package only.
- README rewritten for a first-time reader; docs/SIGNRULES.md restructured with a table of contents and an
  appendix on R4 and R6. `docs/CHANGES.md` is folded into the introduction of SIGNRULES.md.
- A release workflow publishes to PyPI by trusted publishing.

## 0.5.2 (2026-11-23)

- `complete(seed=k)` is back: a random branch order, so `max_completions=1` returns a random completion.
- `signrules --version` and `signrules.__version__`.
- Replay: a malformed certificate step fails the replay with `ReplayError` instead of crashing; the tamper
  test covers all eleven four-species certificates.
- CI runs the README's `uvx` install line against every pushed commit.

## 0.5.1 (2026-11-23)

- uv is the install path: `uv tool install`, `uvx`, `uv sync`; CI on uv from the lock file.
- The app (textual) is part of the core install; the README opens with it.
- Zenodo DOI in the README and CITATION.cff.

## 0.5.0 (2026-11-13): first public release

- One command, `signrules`, with subcommands and a verdict at the end of every run; the full-screen app.
- Python API for a first run: `Ecology.parse`, `SearchResult.verdict`, `minor_signs`.
- The replay tests cover every minimal impossible pattern on five and six species.
