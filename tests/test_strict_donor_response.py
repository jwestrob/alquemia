"""Checks use only the real archived donor paths and their prepared strict cells."""
import json,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import strict_donor_response as s
from affordable_common import read_json,verify,xyz,HA_TO_KCAL
from affordable_common import energy
from mace_site_kinematics import Kinematics

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'workspaces/strict_donor_response_20260923/run_v1/manifest.json'

class ActualPaths(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read_json(MANIFEST)

    def test_full_finite_scope(self):
        v=s.validate(MANIFEST)
        self.assertEqual((v['fresh_calls'],v['exact_reuses'],v['total_cells']),(32,4,36))
        self.assertEqual({(r['case_id'],r['angle_radian']) for r in self.m['reused']},{('4MAE',0.)})

    def test_all_actual_coordinate_maps(self):
        for p in self.m['points']:
            kin=Kinematics(s.data(p['mapping'])['context']);q=np.asarray(p['q']);coord=kin.evaluate(q)[1]
            actual=np.asarray([v[1:] for v in xyz(verify(p['xyz']))])
            np.testing.assert_allclose(coord,actual,atol=1e-12,rtol=0)
            self.assertTrue(p['geometry_checks']['pass'])

    def test_archived_DFT_work_sign_and_scale(self):
        points={(p['case_id'],p['angle_radian'],p['metal']):p for p in self.m['points']}
        old=read_json(verify(self.m['components']))
        for r in old['rows']:
            work={z:(points[r['case_id'],r['angle_radian'],z]['DFT']['energy_hartree']-points[r['case_id'],0.,z]['DFT']['energy_hartree'])*HA_TO_KCAL for z in ('Ca','La')}
            self.assertAlmostEqual(work['Ca']-work['La'],r['DFT_delta_R'],places=8)
            for z in work:self.assertAlmostEqual(work[z],r['endpoints'][z]['DFT'],places=8)

    def test_only_declared_numerical_recipe_change(self):
        for t in self.m['tasks']:
            fresh=verify(t['input']).read_text();old=verify(t['source']['input']).read_text()
            self.assertEqual(fresh.replace(' MaxIter 500\n','').replace(' TolE 1e-10\n',''),old)
            self.assertEqual(verify(t['xyz']).read_bytes(),verify(t['source']['xyz']).read_bytes())
            self.assertEqual(t['multiplicity'],1)
            self.assertIsNone(t['seed_source'])

    def test_actual_complete_receipts_and_native_states(self):
        c=read_json(MANIFEST.parent/'COLLECTION.json')
        self.assertEqual((len(c['rows']),c['complete_cells']),(36,36))
        for r in c['rows']:
            self.assertEqual(r['status'],'complete')
            self.assertEqual(r['initial_guess'],['SAD'])
            self.assertEqual(r['observed_TolE_hartree'],1e-10)
            self.assertTrue(r['details']['native_mixer_observed'])
            self.assertEqual(r['energy_hartree'],energy(verify(r['actual']['output'])))
            self.assertEqual(s.data(r['actual']['receipt'])['parallelism']['nprocs'],1)
        for p in self.m['points']:
            state=p['MACE']['state_check']
            self.assertEqual((state['charge'],state['spin_multiplicity']),(p['charge'],1.))

    def test_full_raw_matrix_work_reconstruction(self):
        c=read_json(MANIFEST.parent/'COLLECTION.json');cells={(r['case_id'],r['angle_radian'],r['metal'],r['medium']):energy(verify(r['actual']['output'])) for r in c['rows']}
        for r in c['endpoint_works']:
            key=(r['case_id'],r['angle_radian'],r['metal']);zero=(r['case_id'],0.,r['metal'])
            work=(cells[key+('alpb',)]-cells[key+('vacuum',)]-cells[zero+('alpb',)]+cells[zero+('vacuum',)])*HA_TO_KCAL
            self.assertAlmostEqual(work,r['strict_solvent'],places=8)
            self.assertAlmostEqual(r['MACE']+work,r['strict_composite'],places=8)
        for r in c['differentials']:
            pair={w['metal']:w for w in c['endpoint_works'] if (w['case_id'],w['angle_radian'])==(r['case_id'],r['angle_radian'])}
            self.assertAlmostEqual(pair['Ca']['strict_composite']-pair['La']['strict_composite'],r['strict_composite'],places=8)

if __name__=='__main__':unittest.main()
