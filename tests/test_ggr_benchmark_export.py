"""Ledger integrity checks using the real 68-row release and executed GGR jobs."""
from copy import deepcopy
from pathlib import Path
import json
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from affordable_common import InvalidArtifact, digest, read_json, verify
from ggr_benchmark_export import export_ledger, response_record


class BenchmarkDevelopmentExport(unittest.TestCase):
    def setUp(self):
        self.parent = ROOT / "workspaces/benchmark_set_20260915/scored_release_1199508/benchmark_manifest.json"
        self.collections = [ROOT / f"workspaces/ggr_mechanism_20260915/stage_{stage}_tasks_v1/collection_{job}.json"
                            for stage, job in (("a", "1199770"), ("b", "1199802"))]
        if not all(path.is_file() for path in [self.parent, *self.collections]):
            self.skipTest("actual scored parent and Stage A/B collections unavailable")
        self.temp = tempfile.TemporaryDirectory(prefix="ggr_export_test_", dir=ROOT / "workspaces")
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def corrupted_copy(self, name, value):
        path = self.directory / ("explicitly_corrupted_" + name + ".json")
        path.write_text(json.dumps(value))
        return path

    def test_all_original_fields_groups_and_scores_preserved(self):
        original = read_json(self.parent)
        old_hash = digest(self.parent)
        result = export_ledger(self.parent, self.collections, self.directory / "export")
        self.assertEqual(digest(self.parent), old_hash)
        self.assertEqual(len(result["rows"]), 68)
        for key in original.keys() - {"rows"}:
            self.assertEqual(result[key], original[key])
        for old, new in zip(original["rows"], result["rows"]):
            self.assertEqual({key: new[key] for key in old}, old)
        rows = {row["target_id"]: row for row in result["rows"]}
        self.assertEqual(len(rows["GGR_1GLG"]["development_scores"]), 6)
        self.assertEqual(len(rows["ALACTA_BOVINE_STRONG_SITE"]["development_scores"]), 2)
        self.assertEqual(len(rows["AEQUORIN_1SL8_VECTOR"]["development_scores"]), 1)
        self.assertEqual(len(rows["AEQUORIN_1SL8_VECTOR"]["existing_scores"]), 6)
        source_rows = {row["case"]: row for path in self.collections for row in read_json(path)["rows"]}
        for target in rows.values():
            for addition in target.get("development_scores", []):
                original_score = source_rows[addition["case"]]
                self.assertEqual(addition["biological_group"], target["biological_group"])
                self.assertEqual(addition["execution_biological_group"], original_score["biological_group"])
                for key in ("score", "endpoints", "protocol_id", "tasks", "decision", "status"):
                    self.assertEqual(addition[key], original_score[key])
                self.assertFalse(addition["site_specific_affinity_label_added"])
                self.assertFalse(addition["independent_biological_observation"])
                self.assertIsNone(addition["calibrated_decision"])
                preparation = read_json(verify(addition["provenance"]["preparation_manifest"]))
                self.assertEqual(addition["preparation_record"], preparation)
        meta = result["development_export"]
        self.assertEqual(meta["development_pair_count"], 9)
        self.assertEqual(meta["completed_endpoints"], 18)
        self.assertEqual(meta["new_calibrated_decisions"], 0)
        self.assertFalse(meta["default_changed"])
        self.assertEqual(read_json(self.directory / "export/benchmark_manifest.json"), result)
        with self.assertRaisesRegex(InvalidArtifact, "new output directory"):
            export_ledger(self.parent, self.collections, self.directory / "export")

    def test_dropped_original_row_and_duplicate_stage_are_rejected(self):
        corrupted = read_json(self.parent)
        corrupted["rows"].pop()  # Intentionally corrupted real ledger; no invented label.
        path = self.corrupted_copy("parent_missing_one_record", corrupted)
        with self.assertRaisesRegex(InvalidArtifact, "68 records"):
            export_ledger(path, self.collections, self.directory / "missing_row")
        with self.assertRaisesRegex(InvalidArtifact, "duplicate stage"):
            export_ledger(self.parent, [self.collections[0]] * 2, self.directory / "duplicate")
        self.assertFalse((self.directory / "missing_row").exists())
        self.assertFalse((self.directory / "duplicate").exists())

    def test_corrupted_actual_score_and_relabeling_are_rejected(self):
        collection = read_json(self.collections[0])
        corrupted = deepcopy(collection)
        corrupted["rows"][0]["score"]["S_kcal_mol"] = 0.0  # Deliberate corruption, never scientific evidence.
        path = self.corrupted_copy("actual_score", corrupted)
        with self.assertRaisesRegex(InvalidArtifact, "score algebra"):
            export_ledger(self.parent, [path], self.directory / "bad_score")
        corrupted = deepcopy(collection)
        corrupted["rows"][0]["target_id"] = "ALACTA_BOVINE_STRONG_SITE"
        path = self.corrupted_copy("target_mapping", corrupted)
        with self.assertRaisesRegex(InvalidArtifact, "differs from its task manifest"):
            export_ledger(self.parent, [path], self.directory / "bad_target")

    def test_real_available_sensitivity_artifacts_stay_separate(self):
        manifest = ROOT / "workspaces/ggr_mechanism_20260915/stage_c_tasks_v1/manifest.json"
        if not manifest.exists():
            self.skipTest("actual approved Stage C task manifest unavailable")
        from ggr_sensitivity import collect
        # Read existing executed endpoints only. This never launches or fabricates work.
        actual = collect(manifest)
        path = self.directory / "actual_stage_c_collection_snapshot.json"
        path.write_text(json.dumps(actual))
        result = export_ledger(self.parent, self.collections, self.directory / "with_sensitivity", [path])
        ggr = next(row for row in result["rows"] if row["target_id"] == "GGR_1GLG")
        self.assertEqual(len(ggr["development_scores"]), 6)
        response = ggr["development_response_records"][0]
        self.assertEqual(response["collection_record"], actual)
        self.assertEqual(response["total_endpoints"], 20)
        self.assertIsNone(response["relaxation_correction_kcal_mol"])
        self.assertIsNone(response["entropy_correction_kcal_mol"])
        self.assertTrue(all("development_response_records" not in row for row in result["rows"] if row is not ggr))
        corrupted = deepcopy(actual)
        corrupted["relaxation_correction_kcal_mol"] = 0.0  # Explicitly corrupt unavailable-to-zero conversion.
        path = self.corrupted_copy("response_correction", corrupted)
        with self.assertRaisesRegex(InvalidArtifact, "cannot supply mechanical corrections"):
            response_record(path)


if __name__ == "__main__":
    unittest.main()
