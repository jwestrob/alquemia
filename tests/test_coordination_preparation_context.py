"""Actual pinned contexts; geometry/parser checks do not imply energy validation."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import coordination_preparation_context as c
from affordable_common import read_json,verify,xyz,InvalidArtifact
from mace_site_kinematics import Kinematics


class RealContextMaps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'workspaces/coordination_preparation_20260919/prepared_v3/manifest.json'
        cls.manifest=read_json(cls.path)

    def test_all_real_source_hashes_and_endpoint_scope(self):
        actual=c.validate(self.path)
        self.assertEqual(actual['tasks'],10)
        self.assertEqual(actual['modes']['1H4I__Ca'],10)
        self.assertEqual(actual['modes']['4MAE__La'],12)

    def test_all_real_maps_preserve_bonds_and_exact_fixed_atoms(self):
        for task in self.manifest['tasks']:
            geo=read_json(verify(task['geometry']))
            q=np.array([.005*(-1)**i for i in range(len(geo['core']['modes']))])
            with self.subTest(task=task['task_id']):
                result=c.proposal_checks(geo,q,xyz(verify(task['context_xyz'])),xyz(verify(task['core']['xyz'])))
                self.assertTrue(result['pass'])
                self.assertTrue(result['water_coordinates_unchanged'])

    def test_trust_constraint_analytic_jacobian_on_all_real_maps(self):
        for task in self.manifest['tasks']:
            kin=Kinematics(read_json(verify(task['geometry']))['context']);n=len(kin.modes)
            q=np.full(n,.008);_,jac=c.constraints(kin,q)
            for i in range(n):
                d=np.zeros(n);d[i]=1e-6
                observed=(c.constraints(kin,q+d)[0]-c.constraints(kin,q-d)[0])/2e-6
                np.testing.assert_allclose(jac[:,i],observed,atol=1e-8,rtol=0)

    def test_unphysical_metal_step_rejected_on_real_fixture(self):
        task=self.manifest['tasks'][0];geo=read_json(verify(task['geometry']))
        q=np.zeros(len(geo['core']['modes']));q[0]=.21
        with self.assertRaisesRegex(InvalidArtifact,'outside physical domain'):
            c.proposal_checks(geo,q,xyz(verify(task['context_xyz'])),xyz(verify(task['core']['xyz'])))

    def test_corrupted_cap_mapping_is_rejected(self):
        task=self.manifest['tasks'][0];geo=copy.deepcopy(read_json(verify(task['geometry'])))
        cap=next(a for a in geo['core']['core_links'] if a[1]=='cap');cap[1]='independent_fake_cap'
        kin=Kinematics(geo['core'])
        with self.assertRaisesRegex(InvalidArtifact,'unsupported physical core link'):
            kin.evaluate(np.zeros(len(kin.modes)))

    def test_pqq_baseline_and_h_prepared_energy_receipts_are_real(self):
        for task in self.manifest['tasks'][:4]:
            core=task['core'];actual=c.endpoint(core['output'],core['receipt'],core['xyz'],core['input'])
            self.assertEqual(actual['energy_hartree'],core['energy_hartree'])
            self.assertIn('pqq_h_dft_v1',core['output']['path'])
        self.assertNotIn('reference',self.manifest)

    def test_actual_native_proposals_preserve_physical_map_and_water(self):
        collection=read_json(self.path.parent/'collection_1202426.json')
        self.assertEqual(len(collection['results']),10)
        for pin in collection['results']:
            result=read_json(verify(pin));task=next(t for t in self.manifest['tasks'] if t['task_id']==result['task_id'])
            self.assertEqual(result['cache_key'],task['cache_key'])
            self.assertEqual(result['status'],'computed')
            self.assertTrue(result['scipy_success'])
            self.assertLessEqual(result['context_energy_eV'],result['source_context_energy_eV']+1e-7)
            geo=read_json(verify(task['geometry']));cr=xyz(verify(task['core']['xyz']))
            checks=c.proposal_checks(geo,result['q'],xyz(verify(task['context_xyz'])),cr)
            self.assertTrue(checks['water_coordinates_unchanged'])
            expected=Kinematics(geo['core']).evaluate(result['q'])[1]
            actual=np.array([a[1:] for a in xyz(verify(result['core_xyz']))])
            np.testing.assert_allclose(actual,expected,atol=5.01e-10,rtol=0)

    def test_changed_geometry_dft_uses_original_recipe_and_charges(self):
        manifest=read_json(ROOT/'workspaces/coordination_preparation_20260919/dft_v1/manifest.json')
        self.assertIsNone(manifest['reference']);self.assertIsNone(manifest['calibrated_decision'])
        self.assertEqual(len(manifest['tasks']),10)
        for t in manifest['tasks']:
            text=verify(t['input']).read_text()
            self.assertEqual(text.splitlines()[0],'! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3')
            self.assertEqual(t['charge'],t['original_core']['charge'])
            self.assertEqual([a[0] for a in xyz(verify(t['xyz']))],
                             [a[0] for a in xyz(verify(t['original_core']['xyz']))])


if __name__=='__main__':unittest.main()
