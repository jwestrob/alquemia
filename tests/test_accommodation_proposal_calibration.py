"""Real executed proposal fixtures; corrupted copies test calibration boundaries."""
import copy
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('proposal_calibration',
    ROOT / 'diagnostics/accommodation_nonlinear_20260920/calibrate_proposals.py')
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class Calibration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / 'workspaces/accommodation_nonlinear_20260920/proposals_v1/final_1204171.json'
        if not path.exists():
            raise unittest.SkipTest('actual executed proposal fixture unavailable')
        cls.result = MOD.read_json(path)
        cls.manifest = MOD.read_json(MOD.verify(cls.result['manifest']))
        cls.old = MOD.read_json(MOD.verify(cls.manifest['frozen_comparison']))

    def test_real_component_scores_and_original_membership(self):
        r = MOD.build(self.result, self.manifest, self.old)
        self.assertEqual(r['status'], 'available')
        self.assertEqual(r['calibration_denominator'], 25)
        self.assertEqual(r['calibration_class_counts'], {'La': 11, 'Ca': 14})
        self.assertAlmostEqual(r['gap_model_kcal_mol'], 5.8443969720974565)
        self.assertFalse(r['noncanonical_folds_used_for_calibration'])
        self.assertEqual(sum(x['role'] == 'consumed_crystal_transfer'
                             for x in r['initial_source_decisions']), 3)

    def test_excluding_transfer_values_leaves_calibration_identical(self):
        # No fabricated energy: remove all real noncalibration source records.
        subset = copy.deepcopy(self.result)
        ids = {r['case_id'] for r in self.old['rows']
               if r['representation'] == 'context' and r['role'] == 'calibration'}
        subset['cases'] = [c for c in subset['cases'] if c['case_id'] in ids]
        subset['endpoints'] = [c for c in subset['endpoints'] if c['case_id'] in ids]
        self.assertEqual(MOD.build(subset, self.manifest, self.old)['bands'],
                         MOD.build(self.result, self.manifest, self.old)['bands'])

    def test_incomplete_calibration_cannot_supply_bands(self):
        bad = copy.deepcopy(self.result)
        bad['cases'][0]['status'] = 'unavailable'
        r = MOD.build(bad, self.manifest, self.old)
        self.assertIsNone(r['bands'])
        self.assertEqual(r['status'], 'unavailable_calibration_member')

    def test_corrupted_label_selection_or_identity_rejected(self):
        for change in ('label', 'selection', 'duplicate'):
            with self.subTest(change=change):
                bad = copy.deepcopy(self.result)
                if change == 'label': bad['cases'][0]['expected_class'] = 'Ca'
                elif change == 'selection':
                    e = bad['endpoints'][0]['selection']
                    e['selected'] = 'proposal' if e['selected'] == 'origin' else 'origin'
                else: bad['cases'].append(copy.deepcopy(bad['cases'][0]))
                with self.assertRaises(MOD.InvalidArtifact):
                    MOD.build(bad, self.manifest, self.old)


if __name__ == '__main__':
    unittest.main()
