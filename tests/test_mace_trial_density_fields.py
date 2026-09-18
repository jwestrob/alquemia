"""Real prepared observations, cache separation, and actual utility results."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from mace_trial_density_fields import validate,key
from mace_qm_field import accepted

W=ROOT/'workspaces/mace_omol_20260917/trial_density_gk_fields_v1'


class TrialObservations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('real trial observation preparation absent')
        cls.m=validate(W/'manifest.json')

    def test_real_paired_centers_and_stencil_blocks(self):
        for t in self.m['tasks']:
            d=read_json(verify(t['probes']));p=np.loadtxt(verify(t['points']),skiprows=1);n=len(d['physical_ids'])
            self.assertEqual(p.shape,(49*n,3))
            np.testing.assert_allclose(p[:n],d['positions_bohr'],atol=2e-13,rtol=0)
            self.assertEqual(len(set(d['physical_ids'])),n)
            self.assertFalse(set(d['physical_ids'])&set(d['source_support_ids']))
            # All derivative pairs have the same actual physical center.
            e=p[37*n:].reshape(2,n,3,2,3)
            for step in range(2):
                for axis in range(3):
                    np.testing.assert_allclose(.5*(e[step,:,axis,0]+e[step,:,axis,1]),p[:n],atol=2e-13,rtol=0)

    def test_vacuum_receipt_cannot_satisfy_responsive_task(self):
        for t in self.m['tasks']:
            old=W.parent/'density_multipole_coupling_v1/execution'/t['task_id']/'attempt_0001'
            self.assertIsNone(accepted(old,t,W/'manifest.json'))
            damaged=copy.deepcopy(t);damaged['files']['gbw']['sha256']='explicitly-corrupted-real-wavefunction-pin'
            self.assertNotEqual(key(damaged,self.m),t['cache_key'])

    def test_actual_eight_utilities_and_numeric_checks(self):
        reports=sorted(W.glob('report_job_*/result.json'))
        if not reports:self.skipTest('eight actual responsive density utility outputs unavailable')
        r=read_json(reports[-1]);self.assertTrue(r['complete']);self.assertEqual(len(r['rows']),8)
        for t in self.m['tasks']:
            row=r['rows'][t['task_id']]
            self.assertIsNotNone(accepted(verify(row['execution_receipt']).parent,t,W/'manifest.json'))
            with np.load(verify(row['arrays'])) as arr:
                self.assertTrue(all(np.all(np.isfinite(arr[name])) for name in arr.files))
                self.assertEqual(arr['H'].shape,(len(arr['phi']),3,3))
                self.assertEqual(arr['exact'].shape,(len(arr['phi']),3))
        expected=all(x['numerical_pass'] for x in [*r['rows'].values(),*r['pairs'].values()])
        self.assertEqual(r['numerical_pass'],expected)
        self.assertIsNone(r['numerical_score'])


if __name__=='__main__':unittest.main()
