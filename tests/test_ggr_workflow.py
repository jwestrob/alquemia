"""Task/collection checks using the actual archived GGR preparation only."""
from pathlib import Path
import json
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from affordable_common import InvalidArtifact, read_json, record, verify
from affordable_workflow import dry_run
from affordable_benchmark import collect
from ggr_workflow import prepare


class GGRWorkflow(unittest.TestCase):
    def setUp(self):
        self.source = ROOT / "workspaces/affordable_challenger_20260915/verified_repairs/ggr_1glg_GGR/repair_manifest.json"
        if not self.source.exists():
            self.skipTest("real archived GGR preparation unavailable")
        self.tmp = tempfile.TemporaryDirectory(prefix="ggr_workflow_test_", dir=ROOT / "workspaces")
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.agreement = ROOT / "diagnostics/ggr_mechanism_plan_20260915/AGREEMENT.md"
        prior = read_json(ROOT / "workspaces/baseline_benchmark_20260915/run_v1/manifest.json")
        self.config = {
            "stage": "A", "approved_plan": record(ROOT / "diagnostics/ggr_mechanism_plan_20260915/PLAN.md"),
            "orca": prior["orca"], "cases": [{
                "case": "archived_GGR_software_fixture", "target_id": "GGR_1GLG",
                "source_structure_id": "1GLG", "representation": "formamide",
                "biological_group": "ggr", "evidence_stratum": "direct_same_assay_direction",
                "preparation_manifest": record(self.source), "historical_case": "ggr_1glg_GGR"}]}

    def config_path(self):
        path = self.directory / "test_configuration.json"
        path.write_text(json.dumps(self.config))
        return path

    def test_real_input_copy_and_unrun_scores_are_unavailable(self):
        out = self.directory / "prepared"
        result = prepare(self.config_path(), ROOT, out, self.agreement)
        self.assertEqual(dry_run(out / "manifest.json")["tasks"], 2)
        for task in result["tasks"]:
            self.assertEqual(verify(task["input"]).read_bytes(), verify(task["source_input"]).read_bytes())
            self.assertEqual(verify(task["xyz"]).read_bytes(), verify(task["source_xyz"]).read_bytes())
        unrun = collect(out / "manifest.json")
        self.assertEqual(unrun["completed_endpoints"], 0)
        self.assertEqual(unrun["rows"][0]["status"], "unavailable")
        self.assertIsNone(unrun["rows"][0]["score"])
        self.assertIsNone(unrun["rows"][0]["endpoints"]["Ca"]["energy_hartree"])
        with self.assertRaises(InvalidArtifact):
            prepare(self.directory / "test_configuration.json", ROOT, out, self.agreement)

    def test_corrupted_real_charge_is_rejected(self):
        actual = read_json(self.source)
        actual["outputs"]["Ca"]["charge"] = 0  # explicitly corrupted real manifest
        corrupted = self.directory / "explicitly_corrupted_real_GGR_manifest.json"
        corrupted.write_text(json.dumps(actual))
        self.config["cases"][0]["preparation_manifest"] = record(corrupted)
        with self.assertRaisesRegex(InvalidArtifact, "state mismatch"):
            prepare(self.config_path(), ROOT, self.directory / "rejected", self.agreement)


if __name__ == "__main__":
    unittest.main()
