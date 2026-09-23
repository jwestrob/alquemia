"""Real pinned geometry and archived-energy algebra; no generated science."""
from pathlib import Path
import sys
import unittest
import tempfile
import numpy as np
from scipy.integrate import simpson
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import local_basin_breadth as b
from affordable_common import InvalidArtifact, read_json, verify, xyz
from mace_site_kinematics import Kinematics
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'workspaces/local_basin_breadth_20260922'


class BreadthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.design=read_json(WORK/'grid_v1/design.json')
        cls.archive=read_json(ROOT/'workspaces/accommodation_torsion_20260920/primary_result_v1.json')
        cls.real=[p for p in cls.archive['points'] if p['case_id']=='4MAE' and p['metal']=='Ca'
                  and p['point'] in ('origin','extra_acidic_ligand_homolog_m0p4',
                      'extra_acidic_ligand_homolog_m0p2','extra_acidic_ligand_homolog_0p2',
                      'extra_acidic_ligand_homolog_0p4')]
        cls.real.sort(key=lambda p:p['angle_radian'])

    def test_all_real_coordinates_and_declared_domain(self):
        self.assertEqual(tuple(c['case_id'] for c in self.design['cases']),b.CASES)
        for c in self.design['cases']:
            self.assertEqual(len(c['points']),65)
            self.assertTrue(all(p['geometry_checks']['pass'] for p in c['points']))
            self.assertLessEqual(max(p['geometry_checks']['maximum_heavy_displacement_A'] for p in c['points']),.8+1e-10)
            self.assertEqual(c['points'][32]['angle_radian'],0.)
            self.assertEqual(c['points'][32]['xyz'],c['source_tasks'][0]['xyz'])
            for p in c['points']:
                self.assertEqual(len([v for v in p['full_q'] if v]),0 if p['index']==32 else 1)

    def test_actual_metric_and_jacobian_replay_and_fixed_metals(self):
        for c in self.design['cases']:
            t=c['source_tasks'][0];kin=Kinematics(read_json(verify(t['mapping']))['context'])
            rows=xyz(verify(t['xyz']))
            q,coords,check=b.geometry(kin,c['mode_index'],c['domain_radian'][1],[r[0] for r in rows])
            self.assertTrue(check['pass']);self.assertTrue(check['physical_metric_constant'])
            self.assertTrue(np.array_equal(coords[0],kin.core[0]))
            self.assertNotIn(kin.metal,c['mode']['moving_indices'])

    def test_nested_grids_and_measure_exact_extent(self):
        for c in self.design['cases']:
            q=np.array([p['angle_radian'] for p in c['points']])
            for x in (q,q[::2],q[8:57],q[8:57:2]):
                self.assertAlmostEqual(float(sum(b.simpson_weights(x))),float(x[-1]-x[0]),places=12)
            self.assertAlmostEqual(q[56],.75*q[-1],places=14)

    def test_actual_archived_energy_algebra_and_simpson(self):
        self.assertEqual(len(self.real),5)
        q=np.array([p['angle_radian'] for p in self.real]);e=np.array([p['composite_energy_kcal_mol'] for p in self.real]);e-=e[2]
        for p in self.real:self.assertAlmostEqual(b.component_energy(p),p['composite_energy_kcal_mol'],places=8)
        r=b.integral(q,e);direct=simpson(np.exp(-(e-e.min())/b.RT),x=q)
        self.assertAlmostEqual(r['width_radian'],float(direct),places=13)
        self.assertAlmostEqual(r['F_relative_kcal_mol'],float(e.min()-b.RT*np.log(direct/(2*np.pi))),places=12)
        # This is only algebraic gauge invariance of measured/computed values,
        # not a new scientific curve or fabricated endpoint calculation.
        shifted=b.integral(q,e+1000.)
        self.assertAlmostEqual(shifted['width_radian'],r['width_radian'],places=11)
        self.assertAlmostEqual(shifted['F_relative_kcal_mol']-r['F_relative_kcal_mol'],1000.,places=10)

    def test_actual_Ca_minus_La_breadth_sign_and_component_algebra(self):
        ids={p['point'] for p in self.real}
        la=sorted([p for p in self.archive['points'] if p['case_id']=='4MAE' and p['metal']=='La'
                   and p['point'] in ids],key=lambda p:p['angle_radian'])
        self.assertEqual(len(la),5)
        q=np.array([p['angle_radian'] for p in self.real])
        values=[]
        for points in (self.real,la):
            e=np.array([p['composite_energy_kcal_mol'] for p in points]);e-=e[2]
            values.append(b.integral(q,e))
        ca,la=values
        minima=ca['minimum_relative_kcal_mol']-la['minimum_relative_kcal_mol']
        width=-b.RT*np.log(ca['width_radian']/la['width_radian'])
        self.assertAlmostEqual(ca['F_relative_kcal_mol']-la['F_relative_kcal_mol'],minima+width,places=12)
        self.assertEqual(np.sign(width),np.sign(la['width_radian']-ca['width_radian']))

    def test_missing_corrupted_copy_never_integrates(self):
        e=np.array([p['composite_energy_kcal_mol'] for p in self.real]);q=np.array([p['angle_radian'] for p in self.real]);e[1]=np.nan
        with self.assertRaises(InvalidArtifact):b.integral(q,e)
        with self.assertRaises(InvalidArtifact):b.paired_curves(q,e,e)

    def test_actual_partial_collection_keeps_all_missing_nodes_explicit(self):
        partial=WORK/'scoring_v1/after_MACE_1210107.json'
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'result.json'
            b.analyze(WORK/'grid_v1/design.json',partial,path)
            data=read_json(path)
        self.assertEqual(len(data['cases']),4)
        for c in data['cases']:
            self.assertEqual(c['result']['status'],'conditional_correction_unavailable')
            self.assertEqual(len(c['result']['missing']),128)
            self.assertIsNone(c['result']['qualified_deltaR_kcal_mol'])
            self.assertIsNone(c['result']['qualified_breadth_deltaR_kcal_mol'])

    def test_actual_complete_collection_serializes_and_recomputes_integrals(self):
        collection=WORK/'scoring_v1/after_solvent_3_1210125.json'
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'result.json';b.analyze(WORK/'grid_v1/design.json',collection,path);data=read_json(path)
        self.assertEqual(len(data['cases']),4)
        self.assertEqual(len(data['cases'][0]['result']['missing']),2)
        for c in data['cases'][1:]:
            result=c['result'];self.assertIn('raw_deltaR_kcal_mol',result)
            for z in ('Ca','La'):
                q=np.linspace(*c['domain_radian'],65);e=np.array(result[z]['relative_energy_kcal_mol'])
                direct=simpson(np.exp(-(e-e.min())/b.RT),x=q)
                self.assertAlmostEqual(result[z]['fine']['width_radian'],float(direct),places=10)
            self.assertAlmostEqual(result['raw_deltaR_kcal_mol'],result['raw_minima_deltaR_kcal_mol']+
                                   result['raw_breadth_deltaR_kcal_mol'],places=10)

    def test_scoring_spec_declares_all_nodes_and_origin_alias(self):
        spec=read_json(WORK/'SPECIFICATION_v1.json')
        self.assertEqual(spec['maximum_candidates_per_case'],65)
        for c in spec['cases']:
            self.assertEqual(len(c['candidates']),65)
            self.assertEqual(c['candidates'][32]['id'],'basin_32')
            self.assertTrue(all(x==0 for x in c['candidates'][32]['full_q']))
            self.assertNotIn('origin',[p['id'] for p in c['candidates']])


if __name__=='__main__':unittest.main()
