"""The tests import two scripts as modules: r4_subsumption and rule_ablation."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
