"""Every script parses and answers ``--help`` from an installed package (guards the import rewrites)."""
import ast
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted((ROOT / "scripts").glob("*.py")) + sorted((ROOT / "docs" / "figures").glob("*.py"))
WITH_HELP = [p for p in SCRIPTS if "argparse" in p.read_text()]


@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: p.name)
def test_script_parses(path):
    ast.parse(path.read_text(), filename=str(path))


@pytest.mark.parametrize("path", WITH_HELP, ids=lambda p: p.name)
def test_script_help(path):
    r = subprocess.run([sys.executable, str(path), "--help"], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-800:]
