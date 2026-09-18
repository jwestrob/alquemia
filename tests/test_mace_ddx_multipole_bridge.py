"""Real multipole/state tests and actual native source-interface outputs."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import BOHR_TO_A,read_json,verify
from mace_frozen_response import source_arrays
from mace_ddx_multipole_bridge import distributions,preflight
W=ROOT/'workspaces/mace_omol_20260917'


class MultipoleBridge(unittest.TestCase):
    def test_real_two_state_distributions_reconstruct_original_dipoles(self):
        m=preflight(W/'ddx_multipole_bridge_v1/manifest.json')
        self.assertEqual(m['new_solver_calls'],0)
        for g in m['groups']:
            for e in g['endpoints'].values():
                s=read_json(verify(e['static']));r=read_json(verify(e['response']))
                _,md,mp,_,_=source_arrays(s,r);p=np.array([a['global_'] for a in s['moments']]);d=distributions(s,r)
                a,b=d['average'],d['difference']
                self.assertTrue(np.array_equal(a[:,0],p[:,0]));self.assertTrue(np.all(b[:,0]==0))
                self.assertLess(np.max(abs((a[:,1:4]+b[:,1:4])*BOHR_TO_A-p[:,1:4]-md)),1e-14)
                self.assertLess(np.max(abs((a[:,1:4]-b[:,1:4])*BOHR_TO_A-p[:,1:4]-mp)),1e-14)
                self.assertLess(np.max(abs(a[:,4:]*BOHR_TO_A**2-p[:,4:])),1e-14)
                self.assertTrue(np.all(b[:,4:]==0))

    def test_actual_dense_native_properties_and_rigid_checks(self):
        paths=list((W/'ddx_multipole_bridge_v1').glob('collection_job_*.json'))
        if not paths:self.skipTest('actual native source-property integration is not yet complete')
        r=read_json(paths[-1]);self.assertEqual(r['new_solver_calls'],0);self.assertIsNone(r['new_score'])
        self.assertEqual(len(r['groups']),4)
        for group in r['groups'].values():
            with np.load(verify(group['arrays']),allow_pickle=False) as a:
                for label,state in group['states'].items():
                    if state['status']=='computed':
                        self.assertTrue(np.isfinite(a[label+'_phi']).all())
                        self.assertTrue(np.isfinite(a[label+'_psi']).all())
                        self.assertEqual(state['checks_pass'],state['potential_error_au']<=1e-10 and state['psi_error']<=1e-12)
                    else:self.assertIn('failure_reason',state)
        self.assertEqual(r['checks_pass'],all(c['pass_'] for c in r['checks']))


if __name__=='__main__':unittest.main()
