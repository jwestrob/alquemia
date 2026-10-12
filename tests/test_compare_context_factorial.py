"""Real-archive algebra and explicit corrupted-copy admission tests; no new energies."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from compare_context_factorial import collect

MACE = ROOT / "diagnostics/orca_parameter_audit_20261012/context_MACE/RESULT.json"
LOW = ROOT / "workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json"


@unittest.skipUnless(MACE.exists() and LOW.exists(), "actual completed scientific fixtures unavailable")
class RealCollection(unittest.TestCase):
    def test_different_serializations_of_same_real_geometry_admitted(self):
        result = collect(MACE, LOW)
        self.assertEqual(result["accepted_endpoints"], 40)
        self.assertLessEqual(max(r["composite_coordinate_max_error_A"] for r in result["endpoints"]), 1e-12)
        self.assertTrue(any(r["xyz"] != r["solvent_xyz"] for r in result["endpoints"]))
        self.assertIsNone(result["bands"])

    def test_missing_required_cell_has_no_successful_fallback(self):
        # Explicitly corrupted copy of completed data, not a simulated calculation.
        low = json.loads(LOW.read_text())
        bad = next(r for r in low["rows"] if r["protein"] == "q9z4j7-pqq-la_model"
                   and r["context"] == "original" and r["geometry"] == "expanded_adaptive_La" and r["metal"] == "La")
        bad.update(accepted=False, primary_solvation_hartree=None, failure_reason="explicit_test_corruption")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrupted_real_collection.json"
            path.write_text(json.dumps(low))
            result = collect(MACE, path)
        pools = {r["pool"]: r for r in result["pools"]
                 if r["protein"] == "q9z4j7-pqq-la_model" and r["context"] == "original"}
        self.assertFalse(pools["common_five"]["available"])
        self.assertIsNone(pools["common_five"]["operational"])
        self.assertTrue(pools["original_three"]["available"])

    def test_mismatched_state_rejected(self):
        low = json.loads(LOW.read_text())
        low["rows"][0]["uhf"] += 2
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrupted_real_state.json"
            path.write_text(json.dumps(low))
            with self.assertRaisesRegex(ValueError, "state"):
                collect(MACE, path)


if __name__ == "__main__":
    unittest.main()
