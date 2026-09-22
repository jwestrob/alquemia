"""Prepared real contexts and kinematic checks; no optimization or model mocks."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import adaptive_angular_proposals as adaptive
from affordable_common import InvalidArtifact,read_json,verify,write_new,xyz
from mace_site_kinematics import Kinematics


class AdaptiveAngularPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'workspaces/adaptive_accommodation_20260922/proposals_v1/manifest.json'
        cls.m=read_json(cls.path)

    def test_real_four_contexts_common_modes_and_source_identity(self):
        result=adaptive.validate(self.path)
        self.assertEqual((result['cases'],result['tasks'],result['molecular_calls']),(4,8,0))
        for cid in adaptive.CASES:
            a,b=[next(t for t in self.m['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
            self.assertEqual(a['active_mode_ids'],b['active_mode_ids'])
            self.assertEqual(len(a['active_indices']),4)
            self.assertEqual(xyz(verify(a['xyz']))[1:],xyz(verify(b['xyz']))[1:])

    def test_constraint_jacobian_on_all_real_maps(self):
        q=np.array([.01,-.008,.005,-.003]);step=1e-6
        for t in self.m['tasks']:
            kin=Kinematics(read_json(verify(t['mapping']))['context'])
            values,jac=adaptive.constraints(kin,t,q)
            self.assertTrue(np.all(values>=0))
            numeric=np.column_stack([(adaptive.constraints(kin,t,q+np.eye(4)[j]*step)[0]-
                                      adaptive.constraints(kin,t,q-np.eye(4)[j]*step)[0])/(2*step) for j in range(4)])
            np.testing.assert_allclose(jac,numeric,atol=1e-8,rtol=1e-6)

    def test_angle_box_alone_does_not_qualify_physical_geometry(self):
        violations=0
        for t in self.m['tasks']:
            kin=Kinematics(read_json(verify(t['mapping']))['context']);q=np.full(4,.8)
            if not adaptive.physical_status(kin,t,q)['physical_feasible']:
                violations+=1
                with self.assertRaisesRegex(InvalidArtifact,'physical displacement'):
                    adaptive.final_geometry(kin,t,q,[a[0] for a in xyz(verify(t['xyz']))])
        self.assertGreater(violations,0)

    def test_small_feasible_motion_keeps_inactive_physical_atoms_fixed(self):
        for t in self.m['tasks']:
            kin=Kinematics(read_json(verify(t['mapping']))['context'])
            result=adaptive.final_geometry(kin,t,np.array([.002,-.001,.0015,-.002]),
                                           [a[0] for a in xyz(verify(t['xyz']))])
            self.assertTrue(result['physical_feasible'])
            self.assertTrue(result['all_other_physical_atoms_fixed'])
            self.assertFalse(result['metal_water_PQQ_scaffold_modes_added'])

    def test_explicitly_corrupted_selected_modes_rejected(self):
        bad=copy.deepcopy(self.m);bad['tasks'][0]['active_mode_ids'][0]='metal_x'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'corrupted_real_manifest.json';write_new(path,bad)
            with self.assertRaisesRegex(InvalidArtifact,'selector'):
                adaptive.validate(path)

    def test_actual_unrun_collection_has_no_fabricated_candidate_or_score(self):
        if any(self.path.parent.glob('proposals/*/result.json')):
            self.skipTest('actual optimization already started; original prelaunch receipt retained')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'collection.json';adaptive.collect(self.path,path);r=read_json(path)
            self.assertEqual(r['available_candidates'],0)
            self.assertEqual(len(r['endpoints']),8)
            for e in r['endpoints']:
                self.assertIsNone(e['candidate']);self.assertIsNone(e['candidate_GFN2_ALPB_hartree'])
                self.assertEqual(e['reason'],'not_run')
            for c in r['cases']:self.assertIsNone(c['score'])

    def test_actual_incomplete_shared_pool_snapshot_blocks_execution(self):
        path=ROOT/'workspaces/nikasha_shared_pool_20260922/pilot_v2/after_MACE_1209840.json'
        with self.assertRaisesRegex(InvalidArtifact,'complete compatible'):
            adaptive.pool_gate(self.path,path)


if __name__=='__main__':unittest.main()
