"""Real five-site inputs, complete comparison inventory and actual execution checks."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from mace_density_panel import validate
from mace_density_gk_hybrid import collect
from mace_multisite_density_report import report
BASE=ROOT/'workspaces/mace_omol_20260917'
W=BASE/'multisite_density_native_v1'

class MultisitePanel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('real multisite manifest unavailable')
        cls.m=validate(W/'manifest.json')

    def test_all_sites_and_exact_paired_physical_boundaries(self):
        m=self.m
        self.assertEqual((len(m['static_tasks']),len(m['response_tasks'])),(63,75))
        self.assertEqual((m['new_DFT_calls'],m['new_MACE_calls']),(0,0))
        for case in m['case_ids']:
            for var in ('primary','rigid','radius_minus','radius_plus'):
                ts=[t for t in m['static_tasks'] if t['case_id']==case and t['variant']==var]
                self.assertEqual({t['state'] for t in ts},{'Ca','La','environment'})
                xyz=[np.loadtxt(verify(t['xyz']),skiprows=1,usecols=(2,3,4)) for t in ts]
                for coords in xyz[1:]:np.testing.assert_array_equal(coords,xyz[0])
                boundaries=[read_json(verify(t['boundary'])) for t in ts]
                self.assertTrue(all(b==boundaries[0] for b in boundaries))

    def test_missing_scientific_outputs_remain_unavailable(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'real_manifest_without_outputs.json';p.write_text(json.dumps(self.m))
            r=collect(p,Path(d)/'collection.json')
            self.assertFalse(r['complete']);self.assertFalse(r['numerical_pass'])
            self.assertIsNone(r['calibrated_class']);self.assertIsNone(r['reference'])
            self.assertTrue(all(t['status']=='unavailable' for t in r['tasks'].values()))

    def test_actual_completed_science_and_ordered_report(self):
        paths=sorted(W.glob('collection_job_*.json'))
        if not paths:self.skipTest('multisite native calculations not yet executed')
        p=paths[-1];r=read_json(p)
        self.assertTrue(r['complete'])
        self.assertEqual((r['actual_energy_calls'],r['actual_field_queries'],r['actual_response_solves']),(63,60,75))
        self.assertEqual(len(r['checks']),39)
        for case in r['variants']['primary']['cases'].values():
            self.assertAlmostEqual(case['R_kcal'],case['endpoints']['Ca']['total_kcal']-case['endpoints']['La']['total_kcal'],places=7)
        with tempfile.TemporaryDirectory() as d:
            result=report(p,Path(d)/'report')
        self.assertEqual(len(result['supporting_comparisons']),6)
        self.assertEqual([v['case_id'] for v in result['ordered_vectors']['aequorin']],['AEQ_1SL8_EF1','AEQ_1SL8_EF3','AEQ_1SL8_EF4'])
        self.assertIsNone(result['aequorin_site_classifications'])
        self.assertIsNone(result['calibrated_class'])
        self.assertIn('retained',result['prior_GGR_failure'])

if __name__=='__main__':unittest.main()
