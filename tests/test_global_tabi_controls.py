"""Actual native controls; malformed fixtures are labelled real-output copies."""
from pathlib import Path
import copy
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from affordable_common import InvalidArtifact, read_json, verify
from affordable_tabi import parse_output, validate_native_controls


class ActualNativeControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = ROOT / "workspaces/global_electrostatic_20260916/surfaces_v1/1h4i_qm33_La/primary"
        receipt = directory / "attempt_0001/execution.json"
        if not receipt.exists():
            raise unittest.SkipTest("actual completed primary TABI output unavailable")
        execution = read_json(receipt)
        manifest = read_json(verify(execution["manifest"]))
        artifacts = execution["artifacts"]
        cls.values = parse_output(*(verify(artifacts[k]) for k in ("output", "csv", "headers")))
        cls.settings = manifest["settings"]
        cls.count = sum(line.startswith(("ATOM ", "HETATM "))
                        for line in verify(manifest["charges"]).read_text().splitlines())

    def test_actual_native_controls_and_charge_inventory_match_manifest(self):
        self.assertEqual(validate_native_controls(self.values, self.settings, self.count)["status"], "passed")
        self.assertEqual(self.count, 9141)

    def test_corrupted_real_control_echo_cannot_satisfy_the_declared_refinement(self):
        corrupt = copy.deepcopy(self.values)
        corrupt["native_csv_values"]["tree_degree"] += 2
        with self.assertRaisesRegex(InvalidArtifact, "tree_degree"):
            validate_native_controls(corrupt, self.settings, self.count)

    def test_deleted_real_control_echo_remains_invalid(self):
        corrupt = copy.deepcopy(self.values)
        del corrupt["native_csv_values"]["mesh_density"]
        with self.assertRaisesRegex(InvalidArtifact, "mesh_density"):
            validate_native_controls(corrupt, self.settings, self.count)


if __name__ == "__main__":
    unittest.main()
