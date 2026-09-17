"""Exact real DFT fixtures and response-scope guards; no invented energies."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,write_new,xyz
from mace_omol_response import validate,collect,POINTS

class ActualResponse(unittest.TestCase):
    def setUp(self):
        self.manifest=ROOT/'workspaces/mace_omol_20260917/masked_response_v1/manifest.json'
        if not self.manifest.exists():self.skipTest('requires pinned real mechanics preparation and response manifest')

    def test_forty_actual_DFT_geometries_and_eight_analytic_centers(self):
        self.assertEqual(validate(self.manifest)['tasks'],40)
        m=read_json(self.manifest)
        self.assertEqual(sum(not t['energy_only'] for t in m['tasks']),8)
        for t in m['tasks']:
            self.assertIn(t['point'],POINTS)
            p=read_json(t['core_mapping']['path'])
            self.assertEqual(t['xyz'],p['grids'][t['point']]['endpoints'][t['metal']]['xyz'])
            if t['point']=='center':
                self.assertEqual(len(t['DFT_gradient']['gradient_kcal_mol_per_A']),len(xyz(t['xyz']['path'])))
                self.assertEqual(t['DFT_gradient']['quantity'],'gradient_not_force')

    def test_corrupted_real_manifest_cannot_drop_a_failure_direction(self):
        m=copy.deepcopy(read_json(self.manifest));m['tasks'].pop()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'corrupted_manifest.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'finite inventory'):validate(p)

    def test_actual_complete_response_keeps_numerical_corrections_unavailable(self):
        p=self.manifest.parent.parent/'masked_response_report_v1/result.json'
        if not p.exists():self.skipTest('scientific integration has not completed')
        actual=collect(self.manifest);self.assertEqual(actual['status'],'complete')
        self.assertEqual(len(actual['comparisons']),24)
        self.assertIsNone(actual['relaxation_correction_kcal_mol'])
        self.assertIsNone(actual['calibrated_class'])
        self.assertFalse(actual['predictive_improvement_claimed'])
        # Results may fail any scientific gate; tests must not require a win.
        for r in actual['comparisons']:
            self.assertEqual(set(r['direct_errors']),{'minus','plus'})
            self.assertEqual(set(r['gates']),{'analytic_derivative','curvature','DFT_anchored','direct_response'})

if __name__=='__main__':unittest.main()
