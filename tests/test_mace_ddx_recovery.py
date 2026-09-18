"""Actual solver replay and immutable recovery inputs; no simulated energies."""
from pathlib import Path
import math
import re
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,record,verify,write_new
from mace_ddx_convergence import preflight as logged_preflight
from mace_ddx_recovery import preflight,retained_state
W=ROOT/'workspaces/mace_omol_20260917'


class DDXRecovery(unittest.TestCase):
    def test_actual_logged_replay_and_energy(self):
        p=W/'ddx_convergence_v1/collection_job_1201280.json';r=read_json(p);m=logged_preflight(verify(r['manifest']))
        self.assertEqual(r['status'],'computed');self.assertEqual(m['model']['maxiter'],1200)
        item=dict(kind='logged',receipt=record(p),state='La',arrays=r['arrays'],maxiter=1200)
        row,a=retained_state(item)
        self.assertFalse(row['new_forward_solve_started']);self.assertEqual(row['energy_hartree'],r['energy_hartree'])
        self.assertLess(abs(.5*math.fsum((a['psi']*a['x']).ravel())*HA_TO_KCAL-r['energy_kcal_mol']),1e-8)
        steps=[(int(n),float(e)) for n,e in re.findall(r'Iteration:\s*(\d+) Rel\. diff:\s*(\S+)',verify(r['solver_log']).read_text())]
        starts=[i for i,(n,e) in enumerate(steps) if n==1]+[len(steps)]
        self.assertEqual(len(starts),3)
        self.assertEqual([steps[i-1][0] for i in starts[1:]],[337,79])
        self.assertTrue(all(steps[i-1][1]<=1e-10 for i in starts[1:]))

    def test_real_recovery_reuses_success_and_preserves_state(self):
        p=W/'ddx_source_recovery_v1/manifest.json'
        if not p.exists():self.skipTest('recovery inventory awaits completed original job')
        m=preflight(p);self.assertEqual(m['requested_new_forward_solves']+m['reused_state_roles'],20)
        for g in m['groups']:
            self.assertFalse(set(g['reuse'])&set(g['new_states']))
            for label,item in g['reuse'].items():
                row,a=retained_state(item);self.assertEqual(row['status'],'computed')
                d=read_json(verify(g['source']));q=np.zeros(len(d['radii_A'])) if label=='zero' else np.array(d['endpoints'][label.split('_')[0]]['charge_e'])
                self.assertLess(np.max(abs(a['psi'][0]-np.sqrt(4*np.pi)*q)),1e-12)

    def test_corrupted_real_iteration_recovery_rejected(self):
        p=W/'ddx_source_recovery_v1/manifest.json'
        if not p.exists():self.skipTest('recovery inventory awaits completed original job')
        m=read_json(p);m['model']['solvent_epsilon']=4.
        with tempfile.TemporaryDirectory() as d:
            corrupt=Path(d)/'corrupted_real_recovery.json';write_new(corrupt,m)
            with self.assertRaises(InvalidArtifact):preflight(corrupt)


if __name__=='__main__':unittest.main()
