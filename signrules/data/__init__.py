"""
The classification data shipped with the package, enough to check it from the package alone
(``signrules verify``):

* ``census(4)``: every four-species pattern up to relabelling, 3044 classes, each ``killed_root`` (715) or
  ``possible`` (2329, with a tropical witness ``U``, ``T``).
* ``screen_records(n)`` (n = 5, 6): the candidate patterns of the extension screen, 83 and 143 of them (fields
  ``code``, ``a``, ``B``, ``status`` (``killed_root`` / ``open``), ``info``, ``seconds``).
* ``root_kills(n)``: the ``killed_root`` rows, the 24 and 41 minimal impossible patterns; ``signrules replay``
  re-derives the certificate of each.
* ``possible(n)``: a certificate of possibility for each ``open`` row, 59 and 102 of them (fields ``code``,
  ``a``, ``B``, and a tropical witness ``U``, ``T``); ``signrules.possible.check_witness`` checks one exactly.

Records of the experiments behind the papers (rule ablations, the R4 subsumption runs, the rule ladder, the
exhaustive five-species screen) are not part of the package: they are in the ``data/`` folder of the
repository and its Zenodo archive.
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


def census(n: int = 4) -> List[dict]:
    """Every n-species class up to relabelling, killed at the root or with a certificate of possibility (n = 4)."""
    return read_jsonl(path(f"n{n}", "census.jsonl"))


def possible(n: int) -> List[dict]:
    """The certificates of possibility of the ``open`` rows of the screen (n = 5, 6)."""
    return read_jsonl(path(f"n{n}", "possible.jsonl"))


__all__ = ["DATA_DIR", "path", "read_jsonl", "census", "screen_records", "root_kills", "possible"]
