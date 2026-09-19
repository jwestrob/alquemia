"""Physical coordinate/force checks using real completed water-containing states."""
from pathlib import Path
import sys
import unittest
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import hydration_basin as h
from hydration_basin_coordinates import WaterCoordinates,MASSES

class RealWaterBasins(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=ROOT/'workspaces/hydration_occupancy_20260918/dft_v1/collection_1201910.json'
        if not p.exists():raise unittest.SkipTest('actual occupancy calculations absent')
        _,cls.centers,_,_=h.source_centers(p)

    def test_complete_real_state_coverage(self):
        self.assertEqual(len(self.centers),20)
        self.assertEqual(sum(c['water_count'] for c in self.centers.values()),32)
        for c in self.centers.values():
            self.assertEqual(c['charge'],-2 if c['metal']=='Ca' else -1)

    def test_rigid_physical_coordinates_preserve_every_other_atom(self):
        for c in self.centers.values():
            f=WaterCoordinates(h.xyz(h.verify(c['xyz'])),c['groups'])
            q=np.tile([.01,-.02,.03,.06,-.07,.08],len(f.groups));coords=f.positions(q)
            np.testing.assert_array_equal(coords[f.fixed],f.initial[f.fixed])
            for w in f.groups:
                a=f.initial[w['indices']];b=coords[w['indices']]
                np.testing.assert_allclose(np.linalg.norm(a[:,None]-a[None,:],axis=2),
                                           np.linalg.norm(b[:,None]-b[None,:],axis=2),atol=2e-14)

    def test_projected_actual_gradient_includes_rotation_chain_rule(self):
        for key in ('1F6S__11__La','6IP9__110__Ca'):
            c=self.centers[key];a=h.xyz(h.verify(c['xyz']));f=WaterCoordinates(a,c['groups'])
            g=h.checked_gradient(c['DFT_gradient'],c['DFT_result'],a)
            q=np.tile([.01,-.02,.03,.06,-.07,.08],len(f.groups));projected=f.gradient(q,g)
            for i in range(f.dimension):
                d=np.zeros(f.dimension);d[i]=1e-5
                derivative=np.sum(g*(f.positions(q+d)-f.positions(q-d)))/(2e-5)
                self.assertAlmostEqual(projected[i],derivative,places=6)

    def test_mass_metric_matches_actual_atom_jacobian(self):
        for c in self.centers.values():
            f=WaterCoordinates(h.xyz(h.verify(c['xyz'])),c['groups']);j=f.jacobian()
            masses=np.zeros(len(f.initial))
            for w in f.groups:
                for i in w['indices']:masses[i]=MASSES[f.symbols[i]]
            direct=np.einsum('aci,a,acj->ij',j,masses,j)
            np.testing.assert_allclose(direct,f.mass_matrix(),atol=2e-12)
            self.assertTrue(np.min(np.linalg.eigvalsh(direct))>0)

    def test_actual_constant_cartesian_anchor_has_rotational_curvature(self):
        c=self.centers['6IP9__110__La'];a=h.xyz(h.verify(c['xyz']));f=WaterCoordinates(a,c['groups'])
        g=h.checked_gradient(c['DFT_gradient'],c['DFT_result'],a)
        raw=-np.load(h.verify(c['MACE_forces']))*h.EV_TO_KCAL;correction=g-raw
        k=h.projected_hessian(lambda q:f.gradient(q,correction),f.dimension,.0005)
        self.assertGreater(np.linalg.norm(k),0.01)
        self.assertLess(np.linalg.norm(k-k.T),1e-5)
        np.testing.assert_allclose(k[:3],0.,atol=1e-9)

    def test_physical_coordinates_and_actual_gradient_transform_together(self):
        from scipy.spatial.transform import Rotation
        c=self.centers['6IP9__110__Ca'];a=h.xyz(h.verify(c['xyz']));f=WaterCoordinates(a,c['groups'])
        rotation=Rotation.from_rotvec([.4,-.2,.3]).as_matrix();shift=np.array([10.,-7.,4.])
        moved=[(z,*(rotation@np.array(row[1:])+shift)) for z,row in zip(f.symbols,a)]
        rf=WaterCoordinates(moved,c['groups']);q=np.tile([.01,-.02,.03,.06,-.07,.08],len(f.groups))
        rq=(q.reshape(-1,3)@rotation.T).reshape(-1)
        np.testing.assert_allclose(rf.positions(rq),f.positions(q)@rotation.T+shift,atol=2e-14)
        g=h.checked_gradient(c['DFT_gradient'],c['DFT_result'],a)
        expected=(f.gradient(q,g).reshape(-1,3)@rotation.T).reshape(-1)
        np.testing.assert_allclose(rf.gradient(rq,g@rotation.T),expected,atol=3e-12)


class RealBasinDirections(unittest.TestCase):
    def test_actual_coupled_directions_are_rigid_and_independent(self):
        from hydration_basin_validation import choose_directions
        root=ROOT/'workspaces/hydration_basin_20260919/mace_v1'
        files=list((root/'execution/1F6S__11__Ca').glob('attempt_*/result.json'))
        if not files:self.skipTest('actual computed coupled curvature not yet available')
        result=h.read_json(files[-1]);p=h.read_json(root/'preparation.json');c=p['centers']['1F6S__11__Ca']
        f=WaterCoordinates(h.xyz(h.verify(c['xyz'])),c['groups']);directions=choose_directions(f,result)
        self.assertEqual(set(directions),{'soft','response'})
        for v in directions.values():
            for sign in (-1,1):
                coords=f.positions(sign*.05*v)
                self.assertLessEqual(np.max(np.linalg.norm(coords-f.initial,axis=1)),.050000001)
                self.assertLessEqual(np.max(np.linalg.norm((.05*v).reshape(-1,6)[:,3:],axis=1)),.100000001)
                np.testing.assert_array_equal(coords[f.fixed],f.initial[f.fixed])
        a,b=directions.values();m=f.mass_matrix()
        cosine=abs(a@m@b)/np.sqrt((a@m@a)*(b@m@b))
        self.assertLessEqual(cosine,.95000001)

    def test_real_boundary_results_cannot_supply_minimum_entropy(self):
        root=ROOT/'workspaces/hydration_basin_20260919/mace_v1'
        if not (root/'collection_job_1202050.json').exists():self.skipTest('actual coupled pilot absent')
        c=h.collect(root/'manifest.json')
        self.assertEqual(c,h.read_json(root/'collection_job_1202050.json'))
        self.assertIsNone(c['occupancy_probabilities'])
        boundary=[r for r in c['rows'] if not r['minimum_interior']]
        self.assertEqual(len(boundary),16)
        self.assertTrue(all(not r['minimum_eligible'] for r in boundary))

    def test_geometry_audit_of_actual_proposals(self):
        from hydration_basin_report import geometry_audit
        root=ROOT/'workspaces/hydration_basin_20260919/mace_v1'
        if not (root/'collection_job_1202050.json').exists():self.skipTest('actual coupled pilot absent')
        p=h.read_json(root/'preparation.json');c=h.read_json(root/'collection_job_1202050.json')
        for row in c['rows']:
            r=h.read_json(h.verify(row['result']));g=geometry_audit(p['centers'][row['task_id']],r['minimum']['xyz'])
            self.assertTrue(g['fixed_coordinates_exact']);self.assertTrue(g['rigid_water_shapes_preserved'])
            self.assertEqual(len(g['waters']),p['centers'][row['task_id']]['water_count'])

    def test_corrupted_real_manifest_coverage_fails(self):
        path=ROOT/'workspaces/hydration_basin_20260919/mace_v1/manifest.json'
        if not path.exists():self.skipTest('real coupled manifest absent')
        m=h.read_json(path);m['tasks']=m['tasks'][:-1]
        with tempfile.TemporaryDirectory() as directory:
            corrupt=Path(directory)/'explicitly_corrupted_real_manifest.json';h.write_new(corrupt,m)
            with self.assertRaisesRegex(h.InvalidArtifact,'coverage'):h.validate(corrupt)

    def test_missing_execution_receipts_never_become_success(self):
        path=ROOT/'workspaces/hydration_basin_20260919/mace_v1/manifest.json'
        if not path.exists():self.skipTest('real coupled manifest absent')
        with tempfile.TemporaryDirectory() as directory:
            isolated=Path(directory)/'actual_manifest_without_execution.json';h.write_new(isolated,h.read_json(path))
            c=h.collect(isolated)
            self.assertEqual(c['status'],'incomplete')
            self.assertEqual(len(c['rows']),20)
            self.assertTrue(all(r['status']=='unavailable' for r in c['rows']))
            self.assertIsNone(c['occupancy_probabilities'])


class RealNativeValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=ROOT/'diagnostics/hydration_basin_20260919/VALIDATION_v1.json'
        if not p.exists():raise unittest.SkipTest('actual completed native coupled checks absent')
        cls.result=h.read_json(p)

    def test_failed_soft_curvature_remains_failed_while_response_passes(self):
        checks=self.result['direction_checks']
        self.assertEqual(len(checks),8);self.assertFalse(self.result['all_direction_checks_pass'])
        failed=[c for c in checks if not c['pass']]
        self.assertEqual([(c['center_id'],c['kind']) for c in failed],[('1F6S__11__Ca','soft')])
        self.assertAlmostEqual(failed[0]['even_error'],-.009747602821140785,places=10)
        self.assertIsNone(self.result['occupancy_probabilities'])

    def test_response_contribution_sign_uses_actual_paired_energies(self):
        rows={r['center_id']:r for r in self.result['rows'] if r['kind']=='proposal'}
        delta=rows['1F6S__11__Ca']['actual_change_kcal_mol']-rows['1F6S__11__La']['actual_change_kcal_mol']
        self.assertAlmostEqual(delta,-1.745033282177037,places=8)
        self.assertTrue(all(r['energy_pass'] for r in rows.values()))
        self.assertTrue(all(r['actual_change_kcal_mol']<0 for r in rows.values()))

    def test_recenter_uses_native_anchors_and_retains_entropy_failure(self):
        p=ROOT/'workspaces/hydration_basin_20260919/recenter_g1/preparation.json'
        if not p.exists():self.skipTest('actual recenter preparation absent')
        prep=h.read_json(p);self.assertEqual(prep['generation'],1)
        self.assertEqual(len(prep['centers']),4)
        for key,c in prep['centers'].items():
            atoms=h.xyz(h.verify(c['xyz']))
            h.checked_gradient(c['DFT_gradient'],c['DFT_result'],atoms)
            self.assertEqual(c['entropy_curvature_qualified'],key!='1F6S__11__Ca')

    def test_conditional_covariance_obeys_equipartition_with_actual_coupled_curvature(self):
        from hydration_basin_report import thermal_extent
        from scipy.constants import R
        p=h.read_json(h.verify(self.result['preparation']))
        for s in p['selection'].values():
            r=h.read_json(h.verify(s['basin_result']));k=np.array(r['minimum']['curvature']['symmetric_curvature'])
            e=thermal_extent(r['minimum']['curvature'],298.15);cov=np.array(e['covariance'])
            np.testing.assert_allclose(k@cov,(R*298.15/4184)*np.eye(len(k)),atol=1e-13)
            self.assertIsNone(e['harmonic_free_energy_kcal_mol'])
            self.assertIsNone(e['physical_fluctuation_estimate'])

if __name__=='__main__':unittest.main()
