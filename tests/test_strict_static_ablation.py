"""Read-only algebra on actual strict32 output; corrupted copies are parser tests."""
import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL, read_json, record
from mace_hybrid import EV_TO_KCAL
import strict_static_ablation as analysis


class RealStaticAblation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT/'workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json'
        cls.qualification = read_json(cls.path)
        cls.reference = read_json(ROOT/'workspaces/strict_static_ablation_20260923/run_v1/REFERENCE.json')

    def test_actual_origin_algebra_and_units(self):
        for row in self.qualification['rows']:
            matrix = row['branches']['fresh']['matrix']
            ca, la = [matrix[z]['origin']['components'] for z in ('Ca', 'La')]
            direct = (ca['MACE_eV']-la['MACE_eV'])*EV_TO_KCAL + HA_TO_KCAL*(
                (ca['GFN2_ALPB_hartree']-ca['GFN2_vacuum_hartree'])-
                (la['GFN2_ALPB_hartree']-la['GFN2_vacuum_hartree']))
            self.assertAlmostEqual(analysis.origin_score(matrix), direct, places=6)

    def test_only_original_canonical_inputs_calibrate(self):
        actual = self.reference['static_reference']
        expected = [r for r in self.qualification['rows'] if r['role'] == 'calibration']
        self.assertEqual(len(actual['rows']), 25)
        self.assertEqual([r['case_id'] for r in actual['rows']], [r['case_id'] for r in expected])
        self.assertEqual(self.reference['qualification'], record(self.path))
        self.assertFalse(self.reference['transfer_rows_used_for_fit'])
        self.assertEqual(actual['bands']['Ca_max'], max(r['R_model_kcal_mol'] for r in actual['rows'] if r['expected_class'] == 'Ca'))
        self.assertEqual(actual['bands']['La_min'], min(r['R_model_kcal_mol'] for r in actual['rows'] if r['expected_class'] == 'La'))

    def test_corrupted_real_missing_origin_is_unavailable(self):
        matrix = copy.deepcopy(self.qualification['rows'][0]['branches']['fresh']['matrix'])
        matrix['Ca']['origin'].update(status='explicit_corrupted_copy', components=None)
        self.assertIsNone(analysis.origin_score(matrix))


if __name__ == '__main__':
    unittest.main()
