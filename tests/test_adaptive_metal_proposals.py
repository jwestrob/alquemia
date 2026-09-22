"""Real archived coordinates/forces; geometry differences are not new energies."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import adaptive_metal_proposals as joint
from affordable_common import InvalidArtifact,read_json,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics


class JointMetalPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'workspaces/adaptive_metal_20260922/prepared_v1/manifest.json'
        cls.m=read_json(cls.path)
        cls.angular=read_json(verify(cls.m['angular_manifest']))
        cls.q=np.array([.005,-.004,.003,.002,-.003,.0015,-.0025])

    def kin(self,t):return Kinematics(read_json(verify(t['mapping']))['context'])

    def test_exact_sources_common_seven_units_and_original_selector(self):
        result=joint.validate(self.path)
        self.assertEqual((result['cases'],result['tasks'],result['molecular_calls']),(4,8,0))
        for t,old in zip(self.m['tasks'],self.angular['tasks']):
            self.assertEqual(t['active_mode_ids'],joint.METAL_IDS+old['active_mode_ids'])
            self.assertEqual(t['selector'],old['selector'])
            self.assertEqual(t['xyz'],old['xyz'])
            self.assertEqual(t['active_coordinate_units'],['angstrom']*3+['radian']*4)

    def test_constraint_jacobian_all_eight_real_maps(self):
        h=1e-6
        for t in self.m['tasks']:
            k=self.kin(t);values,jac=joint.constraints(k,t,self.q)
            numeric=np.column_stack([(joint.constraints(k,t,self.q+np.eye(7)[j]*h)[0]-
                                      joint.constraints(k,t,self.q-np.eye(7)[j]*h)[0])/(2*h) for j in range(7)])
            self.assertTrue(np.all(values>=0))
            np.testing.assert_allclose(jac,numeric,atol=1e-8,rtol=1e-6)

    def test_real_force_projection_sign_units_and_geometry_chain_rule(self):
        h=1e-6
        for t in self.m['tasks']:
            k=self.kin(t);oldpath=verify(self.m['angular_manifest']).parent/'proposals'/t['task_id']/'result.json'
            old=read_json(oldpath)['origin'];forces=np.load(verify(old['forces']),allow_pickle=False)
            origin=joint.projected_gradient_eV(k,t,np.zeros(7),forces)
            np.testing.assert_allclose(origin[:3],-forces[0],atol=1e-12,rtol=0)
            np.testing.assert_allclose(origin[3:]*EV_TO_KCAL,old['gradient_kcal_mol_rad'],atol=1e-8,rtol=0)
            # Frozen archived Cartesian covector tests only the exact mapping.
            # This is not a finite-difference MACE energy or force calculation.
            numeric=[]
            for j in range(7):
                dx=(k.evaluate(joint.full_q(t,self.q+np.eye(7)[j]*h))[1]-
                    k.evaluate(joint.full_q(t,self.q-np.eye(7)[j]*h))[1])/(2*h)
                numeric.append(-np.sum(forces*dx))
            np.testing.assert_allclose(joint.projected_gradient_eV(k,t,self.q,forces),numeric,atol=2e-8,rtol=1e-6)

    def test_cartesian_box_does_not_replace_physical_metal_sphere(self):
        q=np.array([.6,.6,0,0,0,0,0])
        for t in self.m['tasks']:
            k=self.kin(t);s=joint.physical_status(k,t,q)
            self.assertFalse(s['physical_feasible'])
            self.assertAlmostEqual(s['metal_displacement_A'],np.sqrt(.72))
            with self.assertRaisesRegex(InvalidArtifact,'physical displacement'):
                joint.final_geometry(k,t,q,[a[0] for a in xyz(verify(t['xyz']))])

    def test_fixed_atoms_and_actual_metal_motion(self):
        for t in self.m['tasks']:
            k=self.kin(t);result=joint.final_geometry(k,t,self.q,[a[0] for a in xyz(verify(t['xyz']))])
            self.assertTrue(result['all_other_physical_atoms_fixed'])
            self.assertFalse(result['water_PQQ_scaffold_modes_added'])
            p=k.evaluate(joint.full_q(t,self.q))[0]
            np.testing.assert_allclose(p[k.metal]-k.positions[k.metal],self.q[:3],atol=1e-14,rtol=0)

    def test_explicitly_corrupted_real_coordinate_units_rejected(self):
        bad=copy.deepcopy(self.m);bad['tasks'][0]['active_coordinate_units']=['radian']*7
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_manifest.json';write_new(p,bad)
            with self.assertRaisesRegex(InvalidArtifact,'units'):
                joint.validate(p)

    def test_unrun_manifest_retains_all_missing_candidates_and_scores(self):
        self.assertFalse(any(self.path.parent.glob('proposals/*/result.json')))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'collection.json';joint.collect(self.path,p);r=read_json(p)
            self.assertEqual(r['available_candidates'],0)
            self.assertEqual(len(r['endpoints']),8)
            for e in r['endpoints']:
                self.assertEqual(e['reason'],'not_run');self.assertIsNone(e['candidate'])
                self.assertIsNone(e['candidate_GFN2_vacuum_hartree'])
                self.assertIsNone(e['candidate_GFN2_ALPB_hartree'])
            for c in r['cases']:self.assertIsNone(c['score'])


if __name__=='__main__':unittest.main()
