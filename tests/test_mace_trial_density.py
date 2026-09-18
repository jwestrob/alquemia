"""Saved real embedded states; no substitute scientific backend or output."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import BOHR_TO_A,HA_TO_KCAL,energy,read_json,verify
from density_embedding import parse_potential

W=ROOT/'workspaces/mace_omol_20260917/trial_density_gk_source_audit_v1.json'


class ActualTrialSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not W.exists():raise unittest.SkipTest('actual trial-density source audit absent')
        cls.r=read_json(W)

    def test_old_generating_field_removed_before_current_environment(self):
        for row in self.r['rows'].values():
            task=row['source_task'];qt=row['quantum_task'];groups=read_json(verify(task['probe_groups']))
            pc=np.loadtxt(verify(qt['pointcharges']),skiprows=1)
            points=np.loadtxt(verify(task['points']),skiprows=1)
            phi=parse_potential(verify(row['existing_center_potential']),points)
            ids=groups['environment_indices']
            np.testing.assert_allclose(points[ids]*BOHR_TO_A,pc[:,1:],atol=1e-11,rtol=0)
            direct=np.dot(pc[:,0],phi[ids]);embedded=energy(verify(task['source_output']))
            self.assertAlmostEqual(embedded-direct,row['intrinsic_core_hartree'],places=10)
            self.assertNotEqual(embedded,row['intrinsic_core_hartree'])
            self.assertGreaterEqual(row['intrinsic_polarization_cost_kcal'],-.05)

    def test_actual_wavefunctions_and_exact_nuclei_are_available_without_new_calls(self):
        self.assertEqual(len(self.r['rows']),8)
        for row in self.r['rows'].values():
            self.assertTrue(row['physical_state_exact']);self.assertTrue(row['core_coordinates_exact'])
            self.assertTrue(row['wavefunctions_available'])
            for key,pin in row['source_task']['files'].items():
                self.assertTrue(verify(pin).is_file())
                self.assertEqual(pin['sha256'],row['source_task']['source_wavefunctions'][key]['sha256'])
        self.assertIsNone(self.r['numerical_score'])
        self.assertTrue(all(self.r[k]==0 for k in ('new_DFT_calls','new_charge_fits','new_potential_utilities','new_native_solves')))


if __name__=='__main__':unittest.main()
