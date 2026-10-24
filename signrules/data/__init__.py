"""
Data files shipped with the package.

* ``n5/extension_screen.jsonl``, ``n6/extension_screen.jsonl``: the candidate lists of the extension screens
  (83 and 143 patterns; fields ``code``, ``a``, ``B``, ``status`` (``killed_root`` / ``open``), ``info``,
  ``seconds``). The ``killed_root`` rows are the 24 and 41 minimal impossible patterns.
* ``n5/killed_patterns.json``, ``n5/screen_summary.md``: the exhaustive five-species root screen of 2026-09
  (``scripts/screen_all.py --rules R1,R2,R3,R4,R5,R6``): the 21 kills under R1–R6 with their structure
  (``id``, ``structure``, ``kill_rule`` added by hand), and the summary. R8 kills three more; the 24 are the
  ``killed_root`` rows of ``n5/extension_screen.jsonl`` (``root_kills(5)``).
* ``n3/sign_correlations_exact.jsonl``, ``n4/sign_correlations_exact.jsonl``: pairs of minors whose signs are
  correlated in every completion of a pattern, proved by ``scripts/sign_correlations_exact.py`` (they
  disconnect the set of realizations; used in the repository feasible_ecologies).
* ``rule_ablation/``, ``r4_subsumption/``, ``rule_ladder.jsonl``: records of ``scripts/rule_ablation.py``,
  ``scripts/r4_subsumption.py`` and ``scripts/rule_ladder.py``.
"""
import json
from pathlib import Path
from typing import List

DATA_DIR = Path(__file__).resolve().parent


def path(*parts) -> Path:
    """A file of the data directory, e.g. ``path("n5", "extension_screen.jsonl")``."""
    return DATA_DIR.joinpath(*parts)


def read_jsonl(p) -> List[dict]:
    """The records of a JSON-lines file."""
    with open(p) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def screen_records(n: int) -> List[dict]:
    """The rows of the extension screen for n species (n = 5, 6)."""
    return read_jsonl(path(f"n{n}", "extension_screen.jsonl"))


def root_kills(n: int) -> List[dict]:
    """The minimal impossible patterns on n species: the ``killed_root`` rows of the screen."""
    return [r for r in screen_records(n) if r["status"] == "killed_root"]


def killed_patterns() -> List[dict]:
    """The 21 root kills of the exhaustive five-species screen, with ``id``, ``structure``, ``a``, ``B``."""
    with open(path("n5", "killed_patterns.json")) as fh:
        return json.load(fh)


__all__ = ["DATA_DIR", "path", "read_jsonl", "screen_records", "root_kills", "killed_patterns"]
