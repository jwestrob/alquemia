"""Actual archived pool/schema/reference tests; no fabricated cross energies."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import nikasha_pool_compare as compare
from affordable_common import InvalidArtifact, read_json


class PoolComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partial = ROOT / 'workspaces/nikasha_shared_pool_20260922/pilot_v2/after_MACE_1209840.json'
        if not cls.partial.exists(): raise unittest.SkipTest('actual completed MACE-stage collection unavailable')
        cls.bundle = compare.load_results([cls.partial])
        cls.original_pin, cls.original = compare.original_reference(cls.bundle)
        cls.own_reference = read_json(ROOT / 'workspaces/accommodation_nonlinear_20260920/proposal_calibration_v1/REFERENCE.json')

    def test_actual_partial_matrix_is_unavailable_without_discarding_old_score(self):
        self.assertEqual(len(self.bundle['rows']), 4)
        for c in self.bundle['rows'].values():
            self.assertEqual(c['pool']['status'], 'unavailable')
            self.assertIsNotNone(c['old_result']['R_selected'])
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'inspect.json'; compare.inspect([self.partial], out)
            data = read_json(out)
            self.assertEqual(sum(r['expected_class'] is None for r in data['rows']), 2)
            self.assertTrue(all(r['variants']['operational']['R'] is None for r in data['rows']))

    def test_canonical_members_come_only_from_released_record(self):
        for variant in compare.VARIANTS:
            rows = compare.calibration_rows(self.bundle, self.original, variant)
            self.assertEqual(len(rows), 25)
            self.assertEqual({r['case_id'] for r in rows}, {r['case_id'] for r in self.own_reference['rows']})
            self.assertTrue(all(r['R_model_kcal_mol'] is None for r in rows))
            self.assertEqual(compare.extrema_reference(rows, variant)['status'], 'unavailable_calibration_member')

    def test_partial_pilot_cannot_create_a_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(InvalidArtifact, 'pilot4'):
                compare.calibrate([self.partial], ROOT / 'diagnostics/nikasha_shared_pool_20260922/PLAN.md', Path(tmp) / 'reference.json')

    def test_extrema_algebra_replays_real_prior_canonical_values(self):
        # Prior own-proposal values test only the unchanged extrema algebra.
        # They are never written or represented as newly computed pool energies.
        result = compare.extrema_reference(self.own_reference['rows'], 'operational')
        self.assertEqual(result['gap_model_kcal_mol'], self.own_reference['gap_model_kcal_mol'])
        self.assertEqual(result['bands']['Ca_max'], self.own_reference['bands']['Ca_supported_max_R_model_kcal_mol'])
        self.assertEqual(result['bands']['La_min'], self.own_reference['bands']['La_supported_min_R_model_kcal_mol'])
        missing = copy.deepcopy(self.own_reference['rows']); missing[0]['R_model_kcal_mol'] = None
        result = compare.extrema_reference(missing, 'mathematical')
        self.assertIsNone(result['bands']); self.assertEqual(result['available_calibration'], 24)

    def test_explicitly_corrupted_real_collection_cannot_claim_a_pool_success(self):
        damaged = read_json(self.partial)
        damaged['cases'][0]['pool']['mathematical'] = damaged['cases'][0]['old_result']['R0']
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'CORRUPTED_real_collection.json'; path.write_text(json.dumps(damaged))
            with self.assertRaisesRegex(InvalidArtifact, 'selection/score algebra'):
                compare.load_results([path])

    def test_missing_reference_preserves_actual_raw_score(self):
        row = self.own_reference['rows'][0]
        result = compare.method_record(row['R_model_kcal_mol'], None, row['expected_class'])
        self.assertEqual(result['R'], row['R_model_kcal_mol'])
        self.assertEqual(result['score_status'], 'available')
        self.assertEqual(result['reference_status'], 'unavailable')
        self.assertEqual(result['decision'], 'unavailable')

    def test_actual_completed_pool_matrix_when_present(self):
        candidates = sorted((self.partial.parent).glob('after_solvent_*.json'))
        candidates += sorted((self.partial.parent).glob('final*.json'))
        if not candidates: self.skipTest('fresh cross-solvent integration has not completed')
        bundle = compare.load_results([candidates[-1]])
        if not any(c['pool']['status'] == 'available' for c in bundle['rows'].values()):
            self.skipTest('no complete actual cross matrix yet')
        for c in bundle['rows'].values():
            if c['pool']['status'] != 'available': continue
            for z in ('Ca', 'La'):
                r = c['pool']['rows'][z]; work = r['work_from_origin_kcal_mol']
                self.assertEqual(work[r['mathematical_candidate']]['composite_kcal_mol'], min(w['composite_kcal_mol'] for w in work.values()))


if __name__ == '__main__': unittest.main()
