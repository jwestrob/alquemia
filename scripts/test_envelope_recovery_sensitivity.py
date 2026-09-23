"""Read-only checks on actual primary/recovery results, with malformed copies for rejection."""
import json
import tempfile
import unittest
from pathlib import Path

from affordable_common import InvalidArtifact, record
from envelope_recovery_sensitivity import compare

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "workspaces/motion_envelope_transfer_20260923/COMPARISON_v1.json"
RECOVERY = ROOT / "workspaces/native_failed_cell_recovery_20260923/run_v1/COLLECTION.json"


class ActualRecoveryTests(unittest.TestCase):
    def test_only_one_source_and_three_declared_triples_change(self):
        before = record(PRIMARY)
        result = compare(PRIMARY, RECOVERY)
        self.assertEqual(record(PRIMARY), before)
        self.assertEqual(len(result["rows"]), 104)
        self.assertEqual(len(result["triples"]), 100)
        sources = [r for r in result["rows"] if r["primary"] != r["sensitivity"]]
        triples = [r for r in result["triples"] if r["primary"] != r["sensitivity"]]
        self.assertEqual(len(sources), 1)
        self.assertEqual({r["triple_id"] for r in triples}, {"triple_004", "triple_005", "triple_006"})
        self.assertTrue(all(r["primary"]["outcome"] == "unavailable" and
                            r["sensitivity"]["outcome"] == "correct" for r in sources + triples))
        self.assertEqual(result["counts"]["triples"]["primary"], {"correct": 91, "unavailable": 9})
        self.assertEqual(result["counts"]["triples"]["sensitivity"], {"correct": 94, "unavailable": 6})
        self.assertEqual(result["changed_source"]["pool"]["rows"]["La"]["operational_candidate"], "adaptive_Ca")

    def test_corrupted_actual_recovery_cannot_use_favorable_start(self):
        data = json.loads(RECOVERY.read_text())
        data["reported_seed"] = "adaptive_La"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrupted_seed_policy.json"
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(InvalidArtifact, "fixed-origin recovery"):
                compare(PRIMARY, path)

    def test_deleted_actual_agreement_start_is_rejected(self):
        data = json.loads(RECOVERY.read_text())
        data["rows"].pop()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrupted_missing_start.json"
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(InvalidArtifact, "two declared starts"):
                compare(PRIMARY, path)


if __name__ == "__main__":
    unittest.main()
