"""Algebra on an actual MACE grid; these are NOT scientific model-validation gates."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from mace_hybrid import accepted_attempt,EV_TO_KCAL
from mace_mechanics_assess import matrix,square_extrema,stationary_point
from mace_mechanics import H
MANIFEST=ROOT/'workspaces/mace_mechanics_20260916/core_v1/manifest.json'


@unittest.skipUnless(MANIFEST.exists(),'requires actual finite MACE grid')
class MechanicsAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        m=read_json(MANIFEST);rows={}
        for t in m['tasks']:
            if t['case_id']!='ALPHA_1F6S' or t['metal']!='La':continue
            valid=[accepted_attempt(p,t,MANIFEST) for p in (MANIFEST.parent/'execution'/t['task_id']).glob('attempt_*')]
            valid=[r for r in valid if r is not None]
            if not valid:raise unittest.SkipTest('the real alpha La grid is not yet complete')
            rows[t['point']]=valid[-1]
        if len(rows)!=17:raise unittest.SkipTest('complete real17-point grid required')
        center=rows['center']['energy_eV']
        values={k:(v['energy_eV']-center)*EV_TO_KCAL for k,v in rows.items()}
        cls.a=matrix(values,True);cls.coarse=matrix(values)
        prep=read_json(verify(m['preparation']));case=read_json(verify(prep['cases']['ALPHA_1F6S']))
        jac=np.load(verify(case['core_jacobians']));force=np.load(verify(rows['center']['forces']))
        cls.b=-H*np.einsum('nij,ij->n',jac,force)*EV_TO_KCAL

    def test_exact_square_extrema_bound_dense_grid(self):
        # A dense grid is an independent lower bound on the exact maximum.
        delta=self.a-self.coarse;exact=square_extrema(delta)
        axis=np.linspace(-1,1,201);x,y=np.meshgrid(axis,axis)
        sampled=.5*(delta[0,0]*x*x+2*delta[0,1]*x*y+delta[1,1]*y*y)
        lower=float(np.max(np.abs(sampled)))
        self.assertLessEqual(lower,exact['max_absolute_kcal_mol']+1e-12)
        self.assertLessEqual(exact['max_absolute_kcal_mol']-lower,1e-3*(1+np.linalg.norm(delta)))

    def test_coordinate_sign_convention_cannot_change_extrema_or_stationary_energy(self):
        sign=np.diag([-1.,1.]);a=sign@self.a@sign;b=sign@self.b;coarse=sign@self.coarse@sign
        self.assertAlmostEqual(square_extrema(a)['max_absolute_kcal_mol'],square_extrema(self.a)['max_absolute_kcal_mol'],places=12)
        first=stationary_point(self.a,self.b,self.coarse);second=stationary_point(a,b,coarse)
        self.assertEqual(first['status'],second['status'])
        if first['normalized_optimum'] is not None:
            np.testing.assert_allclose(sign@first['normalized_optimum'],second['normalized_optimum'],atol=1e-12)
        if first['predicted_relaxation_kcal_mol'] is not None:
            self.assertAlmostEqual(first['predicted_relaxation_kcal_mol'],second['predicted_relaxation_kcal_mol'],places=12)
            self.assertLessEqual(first['predicted_relaxation_kcal_mol'],0)
        else:self.assertIsNone(second['predicted_relaxation_kcal_mol'])


if __name__=='__main__':unittest.main()
