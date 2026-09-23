"""Real archived geometry/force checks; no invented scientific calculations."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import six_angle_accommodation as six
import adaptive_angular_proposals as original
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_site_kinematics import Kinematics
from mace_hybrid import EV_TO_KCAL
from affordable_common import HA_TO_KCAL

ROOT=Path(__file__).resolve().parents[1]
MP=ROOT/'workspaces/six_angle_accommodation_20260923/run_v2/manifest.json'


class RealSixAngleFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=read_json(MP)

    def test_private_dimension_and_state(self):
        self.assertEqual(original.SETTINGS['selected_angular_modes'],4)
        self.assertEqual(six.SETTINGS['selected_angular_modes'],6)
        self.assertEqual(len(self.m['tasks']),18)
        for t in self.m['tasks']:
            q=np.asarray(t['warm_cached']['full_q']);active=q[t['active_indices']]
            np.testing.assert_array_equal(six.full_q(t,active),q)
            np.testing.assert_array_equal(active[4:],0)
            with self.assertRaises(InvalidArtifact):six.full_q(t,active[:4])

    def test_coupled_constraint_jacobian_all_real_warm_starts(self):
        maximum=0
        for t in self.m['tasks']:
            k=Kinematics(six.data(t['mapping'])['context']);q=np.asarray(t['warm_cached']['full_q'])[t['active_indices']]
            values,j=six.angular.constraints(k,t,q);h=1e-6
            numerical=[]
            for i in range(6):
                d=np.zeros(6);d[i]=h
                numerical.append((six.angular.constraints(k,t,q+d)[0]-six.angular.constraints(k,t,q-d)[0])/(2*h))
            error=float(np.max(np.abs(j-np.asarray(numerical).T)));maximum=max(maximum,error)
            self.assertLess(error,2e-7)
            positions=k.positions_only(six.full_q(t,q))
            np.testing.assert_allclose(values,.8**2-np.sum((positions[k.heavy]-k.positions[k.heavy])**2,axis=1),atol=1e-12,rtol=0)
        print('maximum_constraint_J_error_A2_per_radian',maximum)

    def test_cached_gradient_reprojection_matches_real_cartesian_work(self):
        t=next(t for t in self.m['tasks'] if 'c5ax' in t['case_id'] and 'sample-3' in t['case_id'] and t['metal']=='La')
        p=six.cached_point(t,'warm_cached');q=np.asarray(p['full_q']);k=Kinematics(six.data(t['mapping'])['context'])
        f=np.load(verify(p['forces']),allow_pickle=False);h=1e-6;direction=np.array([.1,-.2,.3,-.1,.4,-.3]);dq=six.full_q(t,direction)*h
        plus=k.evaluate(q+dq)[1];minus=k.evaluate(q-dq)[1]
        cartesian=-float(np.sum(f*(plus-minus)/(2*h)))*EV_TO_KCAL
        projected=float(np.asarray(p['gradient_kcal_mol_rad'])@direction)
        self.assertAlmostEqual(cartesian,projected,places=5)
        self.assertGreater(abs(p['gradient_kcal_mol_rad'][-1]),1.)

    def test_old_pool_replay_and_missing_new_cell(self):
        for c in self.m['cases']:
            self.assertEqual(six.choose_rows(c['old_matrix'],six.CANDIDATES[:3]),c['old_pool'])
            # An intentionally incomplete copy of a real matrix must not pass
            # as a completed extension, even with an available old baseline.
            bad=copy.deepcopy(c['old_matrix'])
            self.assertEqual(six.choose_rows(bad,six.CANDIDATES)['status'],'unavailable')
        self.assertEqual(sum(len(c['old_matrix'][z])*2 for c in self.m['cases'] for z in ('Ca','La')),108)

    def test_actual_origin_and_warm_receipt_conventions(self):
        for t in self.m['tasks']:
            for name in ('origin_cached','warm_cached'):
                p=six.cached_point(t,name)
                self.assertEqual(p['MACE_eV'],t[name]['MACE_eV'])
                self.assertEqual(len(p['gradient_kcal_mol_rad']),6)
                self.assertEqual(p['coordinate'],t[name]['coordinate'])

    def test_actual_complete_five_geometry_matrix_and_sign(self):
        result=read_json(MP.parent/'scalar_v3/COLLECTION.json')
        self.assertEqual(result['complete_cases'],9)
        self.assertEqual(result['complete_cells'],72)
        for c in result['cases']:
            old=next(x for x in self.m['cases'] if x['case_id']==c['case_id'])
            for z in ('Ca','La'):
                for q in six.CANDIDATES[:3]:self.assertEqual(c['matrix'][z][q],old['old_matrix'][z][q])
                oldmin=min(v['composite_kcal_mol'] for v in old['old_pool']['rows'][z]['work_from_origin_kcal_mol'].values())
                newmin=min(v['composite_kcal_mol'] for v in c['pool']['rows'][z]['work_from_origin_kcal_mol'].values())
                self.assertLessEqual(newmin,oldmin+1e-8)
            for mode in ('mathematical','operational'):
                energies={}
                for z in ('Ca','La'):
                    q=c['pool']['rows'][z][mode+'_candidate'];v=c['matrix'][z][q]['components']
                    energies[z]=v['MACE_eV']*EV_TO_KCAL+(v['GFN2_ALPB_hartree']-v['GFN2_vacuum_hartree'])*HA_TO_KCAL
                self.assertAlmostEqual(energies['Ca']-energies['La'],c['pool'][mode]['composite_R_model_kcal_mol'],delta=1e-7)
        self.assertFalse(result['new_calibration'])
        self.assertEqual(result['reference'],self.m['reference'])

    def test_actual_new_proposals_respect_original_domain(self):
        gpu=read_json(MP.parent/'GPU_COLLECTION.json')
        self.assertEqual(len(gpu['results']),18)
        self.assertEqual(sum(c['status']=='complete' for c in gpu['cells']),36)
        for pin in gpu['results']:
            r=six.data(pin);self.assertEqual(r['status'],'proposal_available')
            self.assertLessEqual(r['final_geometry']['maximum_heavy_displacement_A'],.8+1e-7)
            self.assertTrue(r['final_geometry']['all_other_physical_atoms_fixed'])
            self.assertLessEqual(r['MACE_added_work_kcal_mol'],1e-7*EV_TO_KCAL)


if __name__=='__main__':unittest.main()
