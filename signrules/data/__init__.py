"""
Data files shipped with the package.

* ``n5/extension_screen.jsonl``, ``n6/extension_screen.jsonl``: the candidate lists of the extension screens
  (83 and 143 patterns; fields ``code``, ``a``, ``B``, ``status`` (``killed_root`` / ``open``), ``info``,
  ``seconds``). The ``killed_root`` rows are the 24 and 41 minimal impossible patterns.
* ``n5/extension_not_possible.jsonl``: the 24 five-species kills as (a, B, status).
* ``n5/killed_patterns.json``, ``n5/screen_summary.md``: the exhaustive five-species root screen
  (``scripts/screen_all.py``): the 21 kills with their structure, plus the summary.
* ``n3/``, ``n4/sign_correlations_exact.jsonl``: the exact sign-correlation census (``scripts/sign_correlations_exact.py``).
* ``rule_ablation/``, ``r4_subsumption/``, ``rule_ladder.jsonl``: records of ``scripts/rule_ablation.py``,
  ``scripts/r4_subsumption.py`` and ``scripts/rule_ladder.py``.
"""
import json
from pathlib import Path
from typing import List

DATA_DIR = Path(__file__).resolve().parent


def path(*parts) -> Path:
    return DATA_DIR.joinpath(*parts)


def read_jsonl(p) -> List[dict]:
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
