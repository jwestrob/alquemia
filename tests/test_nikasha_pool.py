"""Real pinned source/matrix checks; no invented molecular energies."""
import copy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import nikasha_pool as pool


class SharedPool(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / 'workspaces/nikasha_shared_pool_20260922/pilot_v2/manifest.json'
        if not path.exists(): raise unittest.SkipTest('actual prepared pilot fixture unavailable')
        cls.manifest = pool.read_json(path)

    def test_actual_common_geometry_and_target_electronic_states(self):
        m = self.manifest
        self.assertEqual({c['case_id'] for c in m['cases']}, set(pool.PILOT))
        self.assertEqual(len(m['tasks']), 8)
        for c in m['cases']:
            self.assertEqual(c['status'], 'prepared')
            self.assertEqual(set(c['matrix']['Ca']), set(c['matrix']['La']))
            for z in ('Ca', 'La'):
                own = c['matrix'][z]['origin']
                for q, cell in c['matrix'][z].items():
                    if cell['reused']: continue
                    t = next(t for t in m['tasks'] if t['task_id'] == cell['task_id'])
                    candidate = next(k for k in c['candidates'] if k['id'] == q)
                    a, b = [pool.xyz(pool.verify(p)) for p in (t['xyz'], candidate['xyz'])]
                    self.assertEqual(a[0][0], z)
                    self.assertEqual(a[0][1:], b[0][1:])
                    self.assertEqual(a[1:], b[1:])
                    self.assertEqual(t['multiplicity'], 1)
                self.assertEqual(own['status'], 'complete')

    def test_identity_pool_replays_actual_original_contrast(self):
        for c in self.manifest['cases']:
            matrix = {z: {'origin': c['matrix'][z]['origin']} for z in ('Ca', 'La')}
            result = pool.choose_rows(matrix, ['origin'])
            self.assertEqual(result['status'], 'available')
            self.assertEqual(result['mathematical'], c['old_result']['R0'])
            self.assertEqual(result['operational'], result['mathematical'])

    def test_missing_actual_cross_cells_are_not_baseline_success(self):
        for c in self.manifest['cases']:
            result = pool.choose_rows(c['matrix'], [q['id'] for q in c['candidates']])
            self.assertEqual(result['status'], 'unavailable')
            self.assertIsNone(result['mathematical'])
            self.assertIsNotNone(c['old_result']['R_selected'])

    def test_coordinate_identity_ignores_only_the_intended_metal(self):
        c = self.manifest['cases'][0]
        a = pool.xyz(pool.verify(c['matrix']['Ca']['origin']['xyz']))
        b = pool.xyz(pool.verify(c['matrix']['La']['origin']['xyz']))
        self.assertTrue(pool.same_geometry(a, b))
        bad = copy.deepcopy(b)
        bad[1] = ('X', *bad[1][1:])  # Explicitly corrupted copy of the real fixture.
        self.assertFalse(pool.same_geometry(a, bad))
        self.assertFalse(pool.same_geometry(a, b[:-1]))

    def test_real_own_cells_obey_existing_selection_policy(self):
        for c in self.manifest['cases']:
            for z in ('Ca', 'La'):
                origin = c['matrix'][z]['origin']['components']
                own = c['matrix'][z]['proposal_' + z]['components']
                delta = pool.relative_components(own, origin)['composite_kcal_mol']
                chosen = 'proposal' if delta < -.1 else 'origin'
                self.assertEqual(c['old_result']['endpoint_selection'][z]['selected'], chosen)


if __name__ == '__main__': unittest.main()
