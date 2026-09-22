"""Actual adaptive candidates and common matrix; no synthetic scientific results."""
import copy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import nikasha_pool as pool


class AdaptivePool(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = ROOT / 'workspaces/adaptive_accommodation_20260922/common_pool_v1'
        path = cls.root / 'manifest.json'
        if not path.exists(): raise unittest.SkipTest('real adaptive pool not prepared locally')
        cls.m = pool.read_json(path)

    def test_identical_candidate_membership_retains_previous_real_pool(self):
        self.assertEqual(self.m['protocol_id'], pool.ADAPTIVE_PROTOCOL)
        self.assertEqual(len(self.m['cases']), 4)
        for c in self.m['cases']:
            self.assertEqual(c['status'], 'prepared')
            ids = [q['id'] for q in c['candidates']]
            self.assertEqual(ids, pool.ADAPTIVE_SETTINGS['candidate_order'])
            self.assertEqual(set(c['matrix']['Ca']), set(ids))
            self.assertEqual(set(c['matrix']['La']), set(ids))
            old_ids = pool.SETTINGS['candidate_order']
            old = pool.choose_rows(c['matrix'], old_ids)
            self.assertEqual(old, c['prior_pool'])

    def test_actual_own_native_receipts_reused_at_exact_coordinates(self):
        reused = [t for t in self.m['tasks'] if t.get('native_reuse')]
        self.assertEqual(len(reused), 8)
        self.assertEqual(self.m['new_MACE_cells'], 8)
        self.assertEqual(self.m['maximum_new_GFN2_calls'], 32)
        for t in reused:
            self.assertEqual(pool.native_reuse(t, self.m)['status'], 'complete')
            damaged = copy.deepcopy(t)  # Explicitly corrupted state of this real task.
            damaged['charge'] += 1
            with self.assertRaises(pool.InvalidArtifact): pool.native_reuse(damaged, self.m)

    def test_missing_expanded_matrix_does_not_inherit_old_success(self):
        for c in self.m['cases']:
            result = pool.choose_rows(c['matrix'], [q['id'] for q in c['candidates']])
            self.assertEqual(result['status'], 'unavailable')
            self.assertIsNone(result['operational'])
            self.assertEqual(c['prior_pool']['status'], 'available')

    def test_actual_completed_expansion_cannot_raise_mathematical_row_minimum(self):
        paths = sorted(self.root.glob('after_solvent_*.json'))
        if not paths: self.skipTest('actual expanded solvent evaluations not complete')
        data = pool.read_json(paths[-1])
        self.assertEqual(data['available'], 3)
        for c in data['cases']:
            if c['pool']['status'] == 'unavailable':
                self.assertEqual(c['case_id'], '4MAE')
                self.assertIsNone(c['pool']['mathematical'])
                self.assertEqual(c['prior_pool']['status'], 'available')
                self.assertEqual(c['matrix']['La']['adaptive_Ca']['low']['vacuum']['status'], 'unavailable')
                continue
            for z in ('Ca', 'La'):
                new = c['pool']['rows'][z]; old = c['prior_pool']['rows'][z]
                a = new['work_from_origin_kcal_mol'][new['mathematical_candidate']]['composite_kcal_mol']
                b = old['work_from_origin_kcal_mol'][old['mathematical_candidate']]['composite_kcal_mol']
                self.assertLessEqual(a, b)


if __name__ == '__main__': unittest.main()
