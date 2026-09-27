"""Real fixture parameter/geometry tests only. No Context or molecular calls."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from metal_environment_component_checks import (load_ledger,build_system,geometry_pool,
    prepare,validate,collect,execute,rotation)
LEDGER=ROOT/'workspaces/metal_environment_response_20260926/hybrid_preparation_v1/ledger_v1/LEDGER.json'
PLAN=ROOT/'diagnostics/metal_environment_response_20260926/COMPONENT_CHECK_PLAN.md'

class RealLedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.ledger,cls.artifacts=load_ledger(LEDGER)

    def test_native_force_support_and_exceptions(self):
        import openmm as mm
        with patch.object(mm,'Context',side_effect=AssertionError('no energy calls in tests')):
            system,count=build_system(self.ledger,self.artifacts,'Ca')
        self.assertEqual(count['cross_replacements'],67)
        self.assertEqual(count['particle_count'],9141)
        self.assertEqual(count['constraints'],0)
        self.assertEqual(count['bonded_counts'],self.ledger['counts']['bonded_retained'])
        forces={f.getName():f for f in system.getForces()}
        cross=forces['QM_MM_LJ_generic'];qs,ms=cross.getInteractionGroupParameters(0)
        self.assertEqual((len(qs),len(ms)),(51,9090))
        self.assertFalse(set(qs)&set(ms))
        self.assertEqual(cross.getNumExclusions(),49903)
        from openmm import unit
        for j,e in enumerate(self.artifacts['exceptions.json']):
            a,b,prod,_,_=forces['MM_Coulomb'].getExceptionParameters(j)
            self.assertEqual([a,b],e['atoms'])
            expected=e['chargeprod_e2'] if e['region']=='MM' else 0.
            self.assertAlmostEqual(prod.value_in_unit(unit.elementary_charge**2),expected,places=14)

    def test_real_geometry_directions(self):
        pool,modes=geometry_pool(self.ledger,self.artifacts)
        self.assertEqual(len(pool),16)
        x=pool['A']
        np.testing.assert_array_equal(x,pool['repeat'])
        np.testing.assert_allclose(pool['rigid'],x@rotation().T+(.173,.117,.231),atol=1e-13,rtol=0)
        for mode,desc in modes.items():
            i=desc['physical_index']
            for h in (.001,.0005):
                plus,minus=pool[f'{mode}_{h:+.4f}'],pool[f'{mode}_{-h:+.4f}']
                changed=np.flatnonzero(np.any(plus!=x,axis=1))
                np.testing.assert_array_equal(changed,[i])
                np.testing.assert_allclose((plus[i]-minus[i])/(2*h),desc['direction'],atol=2e-7,rtol=0)
                if mode=='hydroxyl':
                    pivot=desc['pivot']
                    self.assertAlmostEqual(np.linalg.norm(plus[i]-x[pivot]),np.linalg.norm(x[i]-x[pivot]),places=12)

    def test_prepare_and_missing_collection(self):
        import openmm as mm
        with tempfile.TemporaryDirectory() as td, patch.object(mm,'Context',side_effect=AssertionError('no energy calls')):
            dest=Path(td)/'prepared'
            receipt=prepare(LEDGER,PLAN,dest)
            self.assertEqual(receipt['tasks'],32)
            m=validate(dest/'manifest.json')
            self.assertEqual(m['platform'],'Reference')
            report=collect(dest/'manifest.json')
            self.assertEqual(report['completed'],0)
            self.assertEqual(report['status'],'incomplete_or_failed')
            self.assertEqual(report['checks'],[])
            self.assertFalse(report['full_hybrid_qualified'])

    def test_login_energy_guard(self):
        with patch.dict(os.environ,{},clear=True):
            with self.assertRaisesRegex(ValueError,'Slurm'):
                execute(Path('nonexistent.json'),1)

if __name__=='__main__':unittest.main()
