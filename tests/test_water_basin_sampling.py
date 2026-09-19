"""Real archived water-coordinate fixtures; no fabricated scientific outputs."""
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from water_basin_sampling import (WaterCoordinates, coordinate_of, domain_mask, draws,
    extent, integral, left_jacobian, log_haar, read_json, verify, xyz, SETTINGS, validate, R)


class WaterBasinSamplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'workspaces/water_basins_20260919/sampling_v1/manifest.json'
        cls.manifest=read_json(cls.path)
        cls.preparation=read_json(verify(cls.manifest['preparation']))

    def test_pinned_sources_and_finite_manifest(self):
        self.assertEqual(validate(self.path)['tasks'],4)

    def test_inverse_preserves_real_water_geometries(self):
        for c in self.preparation['centers'].values():
            frame=WaterCoordinates(xyz(verify(c['common_xyz'])),c['groups'])
            actual=np.array([a[1:] for a in xyz(verify(c['xyz']))])
            q=coordinate_of(frame,actual)
            np.testing.assert_allclose(frame.positions(q),actual,atol=2e-12,rtol=0)
            np.testing.assert_array_equal(frame.positions(q)[frame.fixed],frame.initial[frame.fixed])

    def test_equivalent_H_labels_share_physical_frame(self):
        for prefix in ('1F6S__11','6IP9__110'):
            ca=self.preparation['centers'][prefix+'__Ca'];la=self.preparation['centers'][prefix+'__La']
            a=np.array([r[1:] for r in xyz(verify(ca['common_xyz']))]);b=np.array([r[1:] for r in xyz(verify(la['common_xyz']))])
            for ids in la['common_physical_frame_H_permutations']:b[ids]=b[ids[::-1]]
            np.testing.assert_allclose(a,b,atol=1e-9,rtol=0)
            self.assertLess(extent(np.array(la['anchor_q'])[None,:])[0],.2)
        self.assertEqual(len(self.preparation['centers']['1F6S__11__La']['common_physical_frame_H_permutations']),2)

    def test_haar_volume_matches_exact_rotation_jacobian(self):
        for c in self.preparation['centers'].values():
            q=np.array(c['anchor_q'])
            determinant=np.prod([np.linalg.det(left_jacobian(v[3:])) for v in q.reshape(-1,6)])
            self.assertAlmostEqual(np.exp(log_haar(q[None,:])[0]),determinant,places=13)

    def test_actual_covariance_positive_and_draws_reproducible(self):
        for c in self.preparation['centers'].values():
            covariance=np.array(c['proposal_covariance']);self.assertGreater(np.linalg.eigvalsh(covariance).min(),0)
            q,logp,seeds=draws(np.array(c['anchor_q']),covariance);p=np.load(verify(c['samples']))
            np.testing.assert_array_equal(q,p['q']);np.testing.assert_array_equal(logp,p['logp']);np.testing.assert_array_equal(seeds,p['seeds'])
            self.assertEqual(len(q),2048);self.assertTrue(np.isfinite(logp).all())

    def test_nested_physical_domains_and_rejection_denominator(self):
        for c in self.preparation['centers'].values():
            q=np.load(verify(c['samples']))['q'];masks=[domain_mask(q,*d) for d in SETTINGS['domains_A_radian']]
            self.assertTrue(np.all(~masks[0]|masks[1]));self.assertTrue(np.all(~masks[1]|masks[2]))
            self.assertLess(masks[2].sum(),len(q));self.assertGreater(masks[2].sum(),0)
            self.assertEqual(len(q),2048) # rejected draws are retained, never resampled

    def test_actual_domain_has_no_water_permutation_overlap(self):
        for c in self.preparation['centers'].values():
            frame=WaterCoordinates(xyz(verify(c['common_xyz'])),c['groups'])
            self.assertGreater(np.linalg.norm(frame.centers[0]-frame.centers[1]),1.2)
            self.assertLess(2*SETTINGS['domains_A_radian'][-1][1],np.pi)

    def test_corrupted_real_shape_fails(self):
        c=next(iter(self.preparation['centers'].values()));frame=WaterCoordinates(xyz(verify(c['common_xyz'])),c['groups'])
        corrupted=frame.initial.copy();corrupted[frame.groups[0]['indices'][1],0]+=.01
        with self.assertRaises(ValueError):coordinate_of(frame,corrupted)

    def test_rejected_draws_keep_denominator_using_real_archived_energies(self):
        # Pure quadrature algebra. These four archived DFT displacements are not
        # represented as Monte Carlo scientific evidence or new calculations.
        v=read_json(ROOT/'diagnostics/hydration_basin_20260919/VALIDATION_v1.json')
        p=read_json(verify(v['preparation']));tasks={t['task_id']:t for t in p['tasks']}
        rows=[r for r in v['rows'] if r['center_id']=='1F6S__11__Ca' and r['kind'] in ('soft','response')]
        q=np.array([tasks[r['task_id']]['q'] for r in rows]);energies=np.array([r['actual_change_kcal_mol'] for r in rows])
        seeds=np.repeat(SETTINGS['seeds'],2);logp=np.zeros(4)
        first=integral(q,logp,seeds,energies)[-1]['conditional_F_config_kcal_mol']
        c=self.preparation['centers']['1F6S__11__Ca'];samples=np.load(verify(c['samples']))
        outside=samples['q'][~domain_mask(samples['q'],.6,1.1)][:4]
        second=integral(np.r_[q,outside],np.r_[logp,np.zeros(4)],np.r_[seeds,seeds],np.r_[energies,np.full(4,np.nan)])[-1]['conditional_F_config_kcal_mol']
        self.assertAlmostEqual(second-first,(R*298.15/4184)*np.log(2),places=12)

    def test_completed_native_sampling_integral_reproduces(self):
        r=read_json(ROOT/'workspaces/water_basins_20260919/sampling_v1/execution/1F6S__11__Ca/attempt_001/result.json')
        a=np.load(verify(r['samples']));actual=integral(a['q'],a['logp'],a['seeds'],a['changes'])
        self.assertEqual(r['status'],'computed')
        for expected,observed in zip(r['integral'],actual):
            self.assertAlmostEqual(expected['conditional_F_config_kcal_mol'],observed['conditional_F_config_kcal_mol'],places=11)
            self.assertEqual(expected['denominator_all_draws'],2048)
            self.assertEqual(expected['numerical_gate_pass'],observed['numerical_gate_pass'])

    def test_native_validation_uses_exact_selected_real_geometry(self):
        m=read_json(ROOT/'workspaces/water_basins_20260919/native_1F6S_11_Ca_v2/manifest.json')
        p=read_json(verify(m['preparation']))
        for t in m['tasks']:
            c=p['selection'][t['center_id']]['center'];frame=WaterCoordinates(xyz(verify(c['common_xyz'])),c['groups'])
            actual=np.array([a[1:] for a in xyz(verify(t['xyz']))])
            np.testing.assert_allclose(frame.positions(t['q']),actual,atol=1e-12,rtol=0)
            np.testing.assert_array_equal(actual[frame.fixed],frame.initial[frame.fixed])
            self.assertIn(t['kind'],SETTINGS['representative_rules'])

    def test_all_actual_native_results_and_units(self):
        from affordable_common import HA_TO_KCAL
        from hydration_square import endpoint
        final=read_json(ROOT/'diagnostics/water_basins_20260919/export_v1/result.json')
        self.assertEqual(sum(len(c['native_checks']) for c in final['centers']),16)
        for pin in final['validation']:
            collection=read_json(verify(pin));manifest=read_json(verify(collection['manifest']))
            prep=read_json(verify(manifest['preparation']));tasks={t['task_id']:t for t in manifest['tasks']}
            for row in collection['rows']:
                t=tasks[row['task_id']];r=row['result'];center=prep['selection'][row['center_id']]['center']
                self.assertEqual(endpoint(r['output'],r['receipt'],t['xyz'],t['input']),r)
                self.assertAlmostEqual(row['actual_change_kcal_mol'],(r['energy_hartree']-center['DFT_result']['energy_hartree'])*HA_TO_KCAL,places=10)

    def test_failures_and_missing_thermochemistry_remain_explicit(self):
        final=read_json(ROOT/'diagnostics/water_basins_20260919/export_v1/result.json')
        self.assertTrue(all(not c['native_representative_gate_pass'] for c in final['centers']))
        self.assertIsNone(final['occupancy_probabilities']);self.assertIsNone(final['absolute_entropy'])
        self.assertIsNone(final['DFT_reweighted_integral']);self.assertFalse(final['baseline_changed'])


if __name__=='__main__':unittest.main()
