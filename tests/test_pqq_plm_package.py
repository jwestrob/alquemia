"""Packaging checks on actual PLM structures, preparations and molecular results."""
import copy
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import pqq_plm_prepare as package
import pqq_plm_export as export
from affordable_common import verify, xyz

FRESH = ROOT / "workspaces/plm_pqq_delivery_20260923/fresh_two_v1"
OLD = ROOT / "workspaces/pqq_three_source_envelope_20260923/plm_v1"
TABLE = ROOT / "workspaces/nikasha_plm_export_20260922/export_v1/proteins.tsv"
RESULT = ROOT / "workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json"


class PLMPackageTests(unittest.TestCase):
    def test_actual_fresh_preparation_reproduces_archived_chemistry_and_coordinates(self):
        ready = package.read_json(FRESH / "prepared/READY.json")
        self.assertEqual((ready["denominator"], ready["available"]), (2, 2))
        self.assertEqual(ready["source_mode"], "fresh")
        for group in ready["groups"]:
            fresh = package.read_json(verify(group["preparation"]))
            old = package.read_json(OLD / group["protein_id"] / "PREPARATION.json")
            for a, b in zip(fresh["cases"], old["cases"]):
                self.assertEqual(a["state_key"], b["state_key"])
                self.assertEqual(a["protein_key"], b["protein_key"])
                for metal in ("Ca", "La"):
                    x, y = (r["representations"]["context"]["endpoints"][metal] for r in (a, b))
                    self.assertEqual(x["charge"], y["charge"])
                    self.assertEqual(x["multiplicity"], y["multiplicity"])
                    self.assertEqual(xyz(verify(x["xyz"])), xyz(verify(y["xyz"])))
                    self.assertEqual(package.read_json(verify(a["maps"][metal])),
                                     package.read_json(verify(b["maps"][metal])))

    def test_complete_preflight_without_new_energy_execution(self):
        d = package.dry_run(FRESH / "REQUEST.json")
        self.assertEqual(d["prospective_source_protonations"], 6)
        plan = FRESH / "prepared/scoring/plan.json"
        r = package.execution.dry_run(plan)
        self.assertEqual((r["groups"], r["sources"], r["origin_GFN2_tasks"]), (2, 6, 24))
        self.assertFalse((plan.parent / "EXECUTION_STARTED.json").exists())
        self.assertFalse(list(plan.parent.rglob("endpoint.out")))

    def test_batch_export_keeps_all_original_rows_and_actual_scores(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "export"
            result = export.export(TABLE, [RESULT], [FRESH / "prepared/READY.json"], out)
            self.assertEqual((result["proteins"], result["result_groups"], result["not_scored"]), (176, 2, 174))
            old = export.parse_tsv(TABLE.read_bytes())
            new = export.parse_tsv((out / "proteins.tsv").read_bytes())
            for a, b in zip(old["rows"], new["rows"]):
                self.assertEqual(a, {k: b[k] for k in old["fields"]})
            scored = [r for r in new["rows"] if r["candidate_status"] == "experimental_available"]
            self.assertEqual(len(scored), 2)
            self.assertTrue(all(r["candidate_decision"] == "Ca-supported" for r in scored))

    def test_duplicate_real_results_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(ValueError, "Overlapping result proteins"):
                export.export(TABLE, [RESULT, RESULT], [], Path(td) / "bad")

    def test_real_preparation_failure_remains_unavailable_in_export(self):
        # Corrupted copy of a real ledger represents a missing source; no energy is invented.
        d = package.read_json(FRESH / "prepared/READY.json")
        d["groups"] = [copy.deepcopy(d["groups"][0])]
        d["groups"][0].update(status="unavailable", reason="explicitly corrupted fixture: removed member")
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "corrupted_ledger.json"
            path.write_text(json.dumps(d))
            out = Path(td) / "export"
            r = export.export(TABLE, [], [path], out)
            self.assertEqual(r["preparation_unavailable"], 1)
            rows = export.parse_tsv((out / "proteins.tsv").read_bytes())["rows"]
            row = next(x for x in rows if x["target_id"] == d["groups"][0]["protein_id"])
            self.assertEqual(row["candidate_status"], "candidate_unavailable_preparation")
            self.assertEqual(row["candidate_median_R_model_kcal_mol"], "")


if __name__ == "__main__":
    unittest.main()
