"""Real 1H4I geometry/Jacobian and actual archived-gradient tests, no new science."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,record
import metal_environment_force_checks as fc
INPUTS=ROOT/'workspaces/metal_environment_response_20260926/preparation/scout_v3/INPUTS.json'
COLLECTION=ROOT/'workspaces/metal_environment_response_20260926/reference_scout_v1/FINAL_COLLECTION_1219207.json'


class ForceGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not INPUTS.exists() or not COLLECTION.exists():raise unittest.SkipTest('Actual completed reference fixture absent')
        cls.config=read_json(INPUTS);cls.modes=fc.construct_modes(cls.config)
        _,_,cls.x=fc.xyz_data(verify(cls.config['endpoints']['Ca']['xyz']))
        _,cls.q,cls.y=fc.pc_data(verify(cls.config['environments']['A']['pointcharges']))

    def test_twenty_cells_and_exact_repeat(self):
        self.assertEqual(len(fc.inventory()),20);self.assertEqual(len(set(fc.inventory())),20)
        for metal in ('Ca','La'):
            x,p=fc.rendered_geometry(self.config,metal,self.modes,'repeat',0.)
            self.assertEqual(x,verify(self.config['endpoints'][metal]['xyz']).read_text())
            self.assertEqual(p,verify(self.config['environments']['A']['pointcharges']).read_text())

    def test_hydroxyl_identity_and_analytic_tangent(self):
        m=self.modes['mm'];xp,yp=fc.displaced(self.x,self.y,self.modes,'mm',1e-5)
        xm,ym=fc.displaced(self.x,self.y,self.modes,'mm',-1e-5)
        np.testing.assert_array_equal(xp,self.x)
        self.assertEqual(np.flatnonzero(np.any(yp!=self.y,axis=1)).tolist(),[m['h_index']])
        np.testing.assert_allclose((yp[m['h_index']]-ym[m['h_index']])/2e-5,m['tangent'],atol=5e-10,rtol=0)
        self.assertAlmostEqual(np.linalg.norm(yp[m['h_index']]-yp[m['o_index']]),np.linalg.norm(self.y[m['h_index']]-self.y[m['o_index']]),places=12)

    def test_boundary_source_jacobian_and_removed_mm_charge(self):
        m=self.modes['boundary'];self.assertIsNone(m['mm_index'])
        ledger=read_json(verify(self.config['boundary_mapping']))
        glu=next(x for x in ledger['ledgers'] if x['fragment_id']=='A:GLU177')
        self.assertIn(2671,glu['removed_source_indices'])
        xp,yp=fc.displaced(self.x,self.y,self.modes,'boundary',1e-5)
        xm,ym=fc.displaced(self.x,self.y,self.modes,'boundary',-1e-5)
        np.testing.assert_array_equal(yp,self.y)
        self.assertEqual(np.flatnonzero(np.any(xp!=self.x,axis=1)).tolist(),[m['cap_index']])
        np.testing.assert_allclose((xp[m['cap_index']]-xm[m['cap_index']])/2e-5,m['cap_tangent'],atol=1e-9,rtol=0)
        x0,_=fc.displaced(self.x,self.y,self.modes,'boundary',0.)
        np.testing.assert_array_equal(x0,self.x)
        axis=np.asarray(m['omitted_A'])-m['retained_A']
        self.assertAlmostEqual(float(axis@np.asarray(m['direction'])),0.,places=12)

    def test_joint_transform_distances(self):
        x,y=fc.displaced(self.x,self.y,self.modes,'rigid')
        np.testing.assert_allclose(np.linalg.norm(x[:,None,:]-y[None,:20,:],axis=2),np.linalg.norm(self.x[:,None,:]-self.y[None,:20,:],axis=2),atol=3e-14,rtol=0)

    def test_actual_center_gradients_and_missing_map_failure(self):
        rows,_=fc.centers(COLLECTION,INPUTS)
        for row in rows.values():
            x,y=fc.gradients(row,len(self.y));self.assertEqual(x.shape,self.x.shape);self.assertEqual(y.shape,self.y.shape)
        config=copy.deepcopy(self.config)
        with tempfile.TemporaryDirectory() as directory:
            mapping=read_json(verify(config['core_mapping']));mapping=[x for x in mapping if x.get('id')!='cap/A:GLU177']
            path=Path(directory)/'corrupt_real_map.json';path.write_text(json.dumps(mapping));config['core_mapping']=record(path)
            with self.assertRaisesRegex(InvalidArtifact,'Glu177 cap'):fc.construct_modes(config)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'corrupt_real_collection.json';c=read_json(COLLECTION);c['rows']['Ca_A']['status']='failed';path.write_text(json.dumps(c))
            with self.assertRaisesRegex(InvalidArtifact,'unavailable'):fc.centers(path,INPUTS)

    def test_real_manifest_rejects_inventory_corruption(self):
        plan=ROOT/'diagnostics/metal_environment_response_20260926/FORCE_CHECK_PLAN.md'
        if not plan.exists():self.skipTest('Frozen parent agreement unavailable')
        with tempfile.TemporaryDirectory() as directory:
            result=fc.prepare(INPUTS,COLLECTION,plan,Path(directory)/'prepared',workers=20,mpi_ranks=17)
            path=verify(result['manifest']);manifest=read_json(path)
            self.assertEqual(len(manifest['tasks']),20)
            self.assertEqual(manifest['execution_resources'],{'mpi_ranks':17,'concurrent_tasks':20})
            fc.validate(path)
            manifest['tasks'].pop();path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(InvalidArtifact,'inventory changed'):fc.validate(path)


if __name__=='__main__':unittest.main()
