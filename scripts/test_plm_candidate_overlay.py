"""Join checks on actual PLM exports/results; malformed copies only for rejection."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from nikasha_plm_export import ExportError, parse_tsv
from plm_candidate_overlay import FIELDS, build_overlay


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "workspaces/nikasha_plm_export_20260922/export_v1/proteins.tsv"
RESULT = ROOT / "workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json"


class RealOverlayTests(unittest.TestCase):
    def test_original_evidence_unchanged_and_unscored_missing(self):
        original = parse_tsv(TABLE.read_bytes())
        fields, rows, receipt = build_overlay(TABLE, RESULT)
        self.assertEqual(len(rows), 176)
        self.assertEqual(fields[:len(original["fields"])], original["fields"])
        for old, new in zip(original["rows"], rows):
            self.assertEqual(old, {k: new[k] for k in original["fields"]})
        unscored = [r for r in rows if r["candidate_status"] == "candidate_not_scored"]
        self.assertEqual(len(unscored), 174)
        self.assertTrue(all(not r[f] for r in unscored for f in FIELDS if f != "candidate_status"))
        self.assertFalse(receipt["biological_accuracy_evaluated"])
        self.assertEqual(receipt["new_molecular_calls"], 0)

    def test_actual_group_and_work_values_copied(self):
        _, rows, receipt = build_overlay(TABLE, RESULT)
        actual = json.loads(RESULT.read_text())
        by_id = {r["target_id"]: r for r in rows}
        self.assertEqual(receipt["candidate_source_denominator"], 6)
        for group in actual["groups"]:
            row = by_id[group["protein_id"]]
            op = group["variants"]["operational"]
            self.assertEqual(float(row["candidate_median_R_model_kcal_mol"]), op["R_model_kcal_mol"])
            self.assertEqual(json.loads(row["candidate_R_model_kcal_mol_json"]), op["member_R_model_kcal_mol"])
            self.assertEqual(row["candidate_decision"], op["decision"])
            self.assertEqual(json.loads(row["candidate_source_ids_json"]), group["members"])
            self.assertEqual(len(json.loads(row["candidate_source_details_json"])), 3)

    def test_corrupted_real_sequence_hash_rejected(self):
        original = parse_tsv(TABLE.read_bytes())
        target = json.loads(RESULT.read_text())["groups"][0]["protein_id"]
        for row in original["rows"]:
            if row["target_id"] == target:
                row["sequence_sha256"] = "corrupted_real_fixture_hash"
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory) / "corrupted_proteins.tsv"
            with table.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=original["fields"], delimiter="\t")
                writer.writeheader()
                writer.writerows(original["rows"])
            with self.assertRaisesRegex(ExportError, "full-sequence identity mismatch"):
                build_overlay(table, RESULT)

    def test_missing_real_source_rejected(self):
        damaged = json.loads(RESULT.read_text())
        damaged["rows"].pop()
        with tempfile.TemporaryDirectory() as directory:
            result = Path(directory) / "corrupted_result_missing_row.json"
            result.write_text(json.dumps(damaged))
            with self.assertRaisesRegex(ExportError, "Source denominator mismatch"):
                build_overlay(TABLE, result)


if __name__ == "__main__":
    unittest.main()
