"""Integrity/algebra tests on pinned real 1H4I artifacts, not model-output mocks."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import mace_hybrid as mh
from affordable_common import InvalidArtifact, cache_key, energy, read_json, verify, xyz


@unittest.skipUnless((ROOT/'workspaces/mace_hybrid_20260916/pilot_v1/manifest.json').exists(),
                     'requires pinned real prepared 1H4I pilot artifacts')
class RealPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mp = ROOT/'workspaces/mace_hybrid_20260916/pilot_v1/manifest.json'
        cls.m = read_json(cls.mp)
        cls.tasks = {t['task_id']: t for t in cls.m['tasks']}

    def test_complete_manifest_pins(self):
        self.assertEqual(mh.dry_run(self.mp)['tasks'], 12)

    def test_paired_physical_coordinates_and_parity(self):
        for variant in ('primary', 'repeat', 'translate', 'rotate'):
            tasks = [self.tasks[f'1h4i_full_{metal}_{variant}'] for metal in ('La', 'Ca')]
            rows = [xyz(verify(t['xyz'])) for t in tasks]
            self.assertEqual(len(rows[0]), 9141)
            self.assertEqual(tasks[0]['charge']-tasks[1]['charge'], 1)
            self.assertEqual([a[1:] for a in rows[0]], [a[1:] for a in rows[1]])
            self.assertEqual(sum(a[0] != b[0] for a,b in zip(*rows)), 1)
            for t, atoms in zip(tasks, rows):
                self.assertEqual(mh.check_atoms(atoms, t['charge'])['all_electron_count'] % 2, 0)

    def test_real_core_parity_corruption_is_rejected(self):
        t = self.tasks['1h4i_qm33_Ca']
        with self.assertRaises(InvalidArtifact):
            mh.check_atoms(xyz(verify(t['xyz'])), t['charge']+1)

    def test_caps_have_no_full_atom_mapping(self):
        mappings = read_json(verify(self.m['atom_mappings']))
        for rows in mappings.values():
            self.assertTrue(any(a['kind'] == 'cap' for a in rows))
            for a in rows:
                self.assertEqual(a['physical_index'] is None, a['kind'] == 'cap')

    def test_rigid_transform_and_repeat_inventory(self):
        a = np.array([r[1:] for r in xyz(verify(self.tasks['1h4i_full_La_primary']['xyz']))])
        b = np.array([r[1:] for r in xyz(verify(self.tasks['1h4i_full_La_repeat']['xyz']))])
        self.assertTrue(np.array_equal(a, b))
        c = np.array([r[1:] for r in xyz(verify(self.tasks['1h4i_full_La_translate']['xyz']))])
        np.testing.assert_allclose(c-a, np.broadcast_to([10., -7., 3.], a.shape), atol=2e-14, rtol=0)
        matrix = mh.rotation()
        np.testing.assert_allclose(matrix @ matrix.T, np.eye(3), atol=1e-15)
        self.assertAlmostEqual(np.linalg.det(matrix), 1.)
        self.assertNotEqual(self.tasks['1h4i_full_La_primary']['cache_key'], self.tasks['1h4i_full_La_repeat']['cache_key'])

    def test_archive_energy_and_sign_replay(self):
        c = read_json(verify(self.m['archived_endpoints']))
        rows = {r['task_id']: r for r in c['rows']}
        for row in rows.values():
            self.assertEqual(energy(verify(row['output'])), row['energy_hartree'])
        contrasts = {}
        for size in (33, 36):
            contrasts[size] = (rows[f'1h4i_qm{size}_Ca']['energy_hartree'] - rows[f'1h4i_qm{size}_La']['energy_hartree'])*mh.HA_TO_KCAL
        self.assertAlmostEqual(contrasts[36]-contrasts[33], 61.738603357, places=6)
        self.assertAlmostEqual(mh.EV_TO_KCAL, 1.602176634e-19 * 6.02214076e23 / 4184, places=13)
        self.assertIsNone(self.m['reference'])

    def test_missing_attempt_is_never_cached(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            self.assertIsNone(mh.accepted_attempt(p, self.tasks['1h4i_qm33_Ca'], self.mp))
            # Explicit malformed receipt; no invented successful scientific output.
            (p/'receipt.json').write_text(json.dumps({'result': None}))
            self.assertIsNone(mh.accepted_attempt(p, self.tasks['1h4i_qm33_Ca'], self.mp))

    def test_model_settings_change_cache_identity(self):
        base = self.m['model']; original = cache_key(base)
        for field, value in [('dtype', 'float32'), ('pbc_handling', 'pbc'),
                             ('assembly', 'different_assembly'), ('microstate', 'different_microstate'),
                             ('explicit_waters', ['corrupted_inventory'])]:
            changed = copy.deepcopy(base); changed[field] = value
            self.assertNotEqual(cache_key(changed), original)


if __name__ == '__main__':
    unittest.main()
