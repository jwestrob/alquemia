"""Physical radial-coordinate tests on real, pinned water-containing cores."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import hydration_water_motion as h


class WaterMotionFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = ROOT/'workspaces/hydration_occupancy_20260918/dft_v1/collection_1201910.json'
        if not cls.collection.exists():
            raise unittest.SkipTest('actual completed occupancy fixtures absent')
        _, cls.centers, _, _ = h.centers_from(cls.collection)

    def test_rigid_motion_preserves_shape_and_all_other_atoms(self):
        for c in self.centers.values():
            atoms = h.xyz(h.verify(c['xyz'])); jac = h.radial_jacobian(atoms, c['groups'])
            variable = {i for w in c['groups'] if w['role']=='variable' for i in w['indices']}
            for amplitude in h.GRID.values():
                moved = h.displaced(atoms, jac, amplitude)
                for i in range(len(atoms)):
                    if i not in variable:self.assertEqual(atoms[i], moved[i])
                for w in c['groups']:
                    before = np.array([atoms[i][1:] for i in w['indices']])
                    after = np.array([moved[i][1:] for i in w['indices']])
                    np.testing.assert_allclose(np.linalg.norm(before[:,None]-before[None,:],axis=2),
                                               np.linalg.norm(after[:,None]-after[None,:],axis=2),atol=2e-14)
                    if w['role']=='variable':
                        original_distance = np.linalg.norm(before[0]-np.array(atoms[0][1:]))
                        new_distance = np.linalg.norm(after[0]-np.array(atoms[0][1:]))
                        self.assertAlmostEqual(new_distance-original_distance, amplitude, places=12)

    def test_real_gradient_force_sign_and_coordinate_chain_rule(self):
        for c in self.centers.values():
            atoms = h.xyz(h.verify(c['xyz'])); jac = np.array(c['jacobian'])
            gradient = h.checked_gradient(c['DFT_gradient'], c['DFT_result'], atoms)
            delta = np.array([a[1:] for a in h.displaced(atoms, jac, .05)])-np.array([a[1:] for a in atoms])
            self.assertAlmostEqual(float(np.sum(gradient*delta)), .05*c['DFT_g0_kcal_mol_A'], places=10)
            forces = np.load(h.verify(c['MACE_forces']), allow_pickle=False)
            self.assertAlmostEqual(-float(np.sum(forces*jac))*h.EV_TO_KCAL, c['MACE_g0_kcal_mol_A'], places=12)

    def test_coordinate_rotates_with_actual_structure_and_gradient(self):
        from scipy.spatial.transform import Rotation
        c = self.centers[h.CENTERS[0]]; atoms = h.xyz(h.verify(c['xyz']))
        rotation = Rotation.from_rotvec([.2, -.1, .3]).as_matrix()
        transformed = [(a[0], *(rotation@np.array(a[1:])+[12.,-3.,7.])) for a in atoms]
        jac = h.radial_jacobian(atoms, c['groups']); moved = h.radial_jacobian(transformed, c['groups'])
        np.testing.assert_allclose(moved, jac@rotation.T, atol=1e-14)
        g = h.checked_gradient(c['DFT_gradient'], c['DFT_result'], atoms)
        self.assertAlmostEqual(float(np.sum(g*jac)), float(np.sum((g@rotation.T)*moved)), places=10)

    def test_corrupted_real_task_inventory_and_coordinate_are_rejected(self):
        mp = ROOT/'workspaces/hydration_motion_20260918/pilot_v2/mace/manifest.json'
        if not mp.exists():self.skipTest('real prepared motion manifest absent')
        h.validate_mace(mp)
        for corrupt in ('missing_task', 'wrong_displacement'):
            m = h.read_json(mp)
            if corrupt=='missing_task':m['tasks'].pop()
            else:m['tasks'][0]['amplitude_A'] *= -1
            with tempfile.TemporaryDirectory() as directory:
                p = Path(directory)/'explicitly_corrupted_real_manifest.json'; h.write_new(p, m)
                with self.assertRaises(h.InvalidArtifact):h.validate_mace(p)

    def test_real_dft_manifest_is_contained_and_matches_mace_coordinates(self):
        mp = ROOT/'workspaces/hydration_motion_20260918/pilot_v2/dft/manifest.json'
        if not mp.exists():self.skipTest('real prepared motion manifest absent')
        from affordable_workflow import dry_run
        dry_run(mp)
        m = h.read_json(mp)
        other = h.read_json(mp.parent.parent/'mace/manifest.json')
        xyz_by_id = {t['task_id']: h.verify(t['xyz']).read_bytes() for t in other['tasks']}
        self.assertEqual(len(m['tasks']), 8)
        for t in m['tasks']:
            path = h.verify(t['xyz'])
            self.assertTrue(path.is_relative_to(mp.parent))
            self.assertEqual(path.read_bytes(), xyz_by_id[t['task_id']])

    def test_actual_startup_retry_preserves_science_and_invalidates_code_cache(self):
        prior = ROOT/'workspaces/hydration_motion_20260918/pilot_v2/mace/manifest.json'
        retry = ROOT/'workspaces/hydration_motion_20260918/mace_retry_v1/manifest.json'
        if not retry.exists():self.skipTest('real technical retry manifest absent')
        a, b = h.read_json(prior), h.read_json(retry)
        h.validate_mace(retry)
        for field in ('model', 'software', 'preparation', 'motion_tolerances'):
            self.assertEqual(a[field], b[field])
        self.assertEqual(b['previous_manifest'], h.record(prior))
        for old, new in zip(a['tasks'], b['tasks']):
            self.assertNotEqual(old['cache_key'], new['cache_key'])
            self.assertEqual({k:v for k,v in old.items() if k!='cache_key'},
                             {k:v for k,v in new.items() if k!='cache_key'})

    def test_actual_response_contrast_and_missing_occupancy(self):
        rp = ROOT/'diagnostics/hydration_motion_20260918/RESULT_v2.json'
        if not rp.exists():self.skipTest('actual DFT motion integration not completed')
        from affordable_common import energy
        r = h.read_json(rp); c = h.read_json(h.verify(r['DFT_collection']))
        m = h.read_json(h.verify(c['manifest'])); tasks = {t['task_id']:t for t in m['tasks']}
        for pair in r['contrasts']:
            group = {z:next(v for v in r['rows'] if v['case']==pair['case'] and v['metal']==z)
                     for z in ('Ca', 'La')}
            for label in ('m050', 'p050'):
                energies, centers = {}, {}
                for z, row in group.items():
                    energies[z] = energy(Path(tasks[row['center_id']+'__'+label]['output_path']))
                    centers[z] = self.centers[row['center_id']]['DFT_result']['energy_hartree']
                direct = ((energies['Ca']-energies['La'])-(centers['Ca']-centers['La']))*h.HA_TO_KCAL
                self.assertAlmostEqual(direct, pair['actual_delta_R_kcal_mol'][label], places=8)
        self.assertIsNone(r['entropy_correction'])
        self.assertIsNone(r['relaxation_correction'])
        self.assertIsNone(r['occupancy_probabilities'])
        self.assertEqual(r['response_model_status'], 'response_model_not_validated')


if __name__=='__main__':unittest.main()
