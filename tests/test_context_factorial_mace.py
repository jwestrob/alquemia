"""Contract tests use real archived MACE/prepared geometry, never mock energies."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("factorial", ROOT / "scripts/context_factorial_mace.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
FIXTURE = ROOT / "workspaces/orca_parameter_audit_20261012/context_mace_v1/manifest.json"


@unittest.skipUnless(FIXTURE.exists(), "real pinned factorial fixture unavailable")
class RealReuseContract(unittest.TestCase):
    def setUp(self):
        self.manifest = json.loads(FIXTURE.read_text())
        self.task = next(t for t in self.manifest["tasks"] if t["reuse"])
        reuse = self.task["reuse"]
        self.cell = {"MACE": reuse["receipt"], "xyz": reuse["xyz"],
                     "components": {"MACE_eV": reuse["energy_eV"]}}

    def test_real_archive_and_state_match(self):
        result = M.reuse_cell(self.cell, self.task, self.manifest["model"])
        self.assertEqual(result["energy_eV"], self.task["reuse"]["energy_eV"])

    def test_changed_state_cannot_reuse(self):
        task = copy.deepcopy(self.task)
        task["charge"] += 1
        with self.assertRaisesRegex(ValueError, "state differs"):
            M.reuse_cell(self.cell, task, self.manifest["model"])

    def test_changed_method_cannot_reuse(self):
        model = copy.deepcopy(self.manifest["model"])
        model["dtype"] = "float32"
        with self.assertRaisesRegex(ValueError, "model differs"):
            M.reuse_cell(self.cell, self.task, model)

    def test_corrupted_copy_of_real_geometry_cannot_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "explicitly_corrupted_real_fixture.xyz"
            lines = Path(self.task["xyz"]["path"]).read_text().splitlines()
            atom = lines[2].split()
            atom[1] = str(float(atom[1]) + .001)
            lines[2] = " ".join(atom)
            path.write_text("\n".join(lines) + "\n")
            task = copy.deepcopy(self.task)
            task["xyz"] = M.pin(path)
            with self.assertRaisesRegex(ValueError, "coordinates differ"):
                M.reuse_cell(self.cell, task, self.manifest["model"])


if __name__ == "__main__":
    unittest.main()
