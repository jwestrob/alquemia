"""Parser/algebra checks on consumed real source and saved native Dy gradient."""
import sys,unittest,copy
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import read_json,xyz,verify
import metal_environment_compact_forces as m
ROOT=Path(__file__).resolve().parents[1]
MAN=ROOT/'workspaces/metal_environment_response_20260926/lady_compact_exchange_v1/manifest.json'

class RealCompactForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.man=read_json(MAN);cls.origin=cls.man['reused_origin']['rows']['Dy_origin'];cls.om=read_json(verify(cls.man['reused_origin']['manifest']));cls.ref=cls.om['tasks'][0]['xyz'];cls.at=xyz(verify(cls.ref));cls.direction=m.geometry(cls.at)
    def test_real_nearest_oxygen_and_rigid_geometry(self):
        d=self.direction;q=np.array([a[1:] for a in self.at]);self.assertEqual(self.at[d['oxygen_index']][0],'O');self.assertAlmostEqual(np.linalg.norm(d['unit_vector']),1)
        translated=[[a[0],*(q[i]+[.173,.117,.231])] for i,a in enumerate(self.at)];dt=m.geometry(translated);self.assertEqual(d['oxygen_index'],dt['oxygen_index']);np.testing.assert_allclose(d['unit_vector'],dt['unit_vector'],atol=1e-14)
    def test_archived_dy_gradient_sign_and_projection(self):
        r=m.endpoint(self.origin,self.ref,self.direction);self.assertEqual(r['status'],'available');g=np.array(self.origin['gradient_kcal_mol_per_A']);expect=float(g[self.direction['metal_index']]@self.direction['unit_vector']);self.assertAlmostEqual(r['gradient_load_kcal_mol_A'],expect,places=10);self.assertEqual(r['force_load_kcal_mol_A'],-r['gradient_load_kcal_mol_A']);np.testing.assert_allclose(r['translation_gradient_sum_kcal_mol_A'],g.sum(axis=0),atol=1e-10)
    def test_missing_and_corrupted_real_gradient_record(self):
        self.assertEqual(m.endpoint({},self.ref,self.direction)['status'],'unavailable');bad=copy.deepcopy(self.origin);bad['gradient_kcal_mol_per_A'][0][0]+=1;self.assertEqual(m.endpoint(bad,self.ref,self.direction)['status'],'invalid')
if __name__=='__main__':unittest.main()
