"""The tests import two scripts (sign_correlations_exact, r4_subsumption) as modules."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
