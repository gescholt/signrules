"""The tests import three scripts as modules: sign_correlations_exact, r4_subsumption, rule_ablation."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
