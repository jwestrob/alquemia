"""Actual joint proposals; distinguish a failed search from an available pool."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import nikasha_pool as pool


class JointPool(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = ROOT / 'workspaces/adaptive_metal_20260922/common_pool_v1'
        path = cls.root / 'manifest.json'
        if not path.exists():
            raise unittest.SkipTest('actual joint proposal pool is not prepared locally')
        cls.m = pool.read_json(path)
        cls.base = pool.read_json(pool.verify(cls.m['base_pool']))

    def test_all_previous_real_energies_and_candidates_retained(self):
        self.assertEqual(self.m['protocol_id'], pool.JOINT_PROTOCOL)
        self.assertEqual(len(self.m['cases']), 4)
        for c, old in zip(self.m['cases'], self.base['cases']):
            self.assertEqual(c['case_id'], old['case_id'])
            for metal in ('Ca', 'La'):
                for candidate in old['candidates']:
                    name = candidate['id']
                    self.assertEqual(c['matrix'][metal][name]['components'],
                                     old['matrix'][metal][name]['components'])
            if c['status'] == 'prepared':
                ids = [q['id'] for q in c['candidates']]
                self.assertEqual(ids, pool.JOINT_SETTINGS['candidate_order'])
                self.assertEqual(set(c['matrix']['Ca']), set(c['matrix']['La']))
                self.assertEqual(set(ids), set(c['matrix']['Ca']))

    def test_real_unsafe_search_remains_unavailable_without_substitution(self):
        bad = [c for c in self.m['cases'] if c['status'] != 'prepared']
        self.assertEqual([c['case_id'] for c in bad], ['PQQSEQ_83440678cbbd658047c9'])
        self.assertEqual(bad[0]['reason'], 'adaptive candidate unavailable: La')
        self.assertEqual(bad[0]['prior_pool']['status'], 'available')
        self.assertFalse(any(t['case_id'] == bad[0]['case_id'] for t in self.m['tasks']))
        actual = pool.read_json(pool.verify(self.m['adaptive_proposals']))
        failed = next(e for e in actual['endpoints'] if e['case_id'] == bad[0]['case_id'] and e['metal'] == 'La')
        self.assertIsNone(failed['candidate'])
        self.assertEqual(failed['reason'], 'unsupported physical proposal geometry')

    def test_charge_and_coordinates_of_real_reused_native_endpoints(self):
        tasks = [t for t in self.m['tasks'] if t.get('native_reuse')]
        self.assertEqual(len(tasks), 6)
        for t in tasks:
            self.assertEqual(pool.native_reuse(t, self.m)['status'], 'complete')
            self.assertEqual(pool.xyz(pool.verify(t['xyz']))[0][0], t['metal'])
        self.assertEqual(self.m['new_MACE_cells'], 6)
        self.assertEqual(self.m['maximum_new_GFN2_calls'], 24)
        self.assertEqual(self.m['GFN2_maxiter'], 500)

    def test_executed_joint_matrix_and_row_minima(self):
        paths = sorted(self.root.glob('after_solvent_*.json'))
        if not paths:
            self.skipTest('actual joint solvent calculations have not completed')
        d = pool.read_json(paths[-1])
        self.assertEqual(d['available'], 3)
        for c in d['cases']:
            if c['status'] != 'prepared':
                self.assertEqual(c['pool']['status'], 'unavailable')
                self.assertIsNone(c['pool']['operational'])
                continue
            self.assertEqual(c['pool'], pool.choose_rows(c['matrix'], [q['id'] for q in c['candidates']]))
            for metal in ('Ca', 'La'):
                new = c['pool']['rows'][metal]
                old = c['prior_pool']['rows'][metal]
                self.assertLessEqual(new['work_from_origin_kcal_mol'][new['mathematical_candidate']]['composite_kcal_mol'],
                                     old['work_from_origin_kcal_mol'][old['mathematical_candidate']]['composite_kcal_mol'])


if __name__ == '__main__':
    unittest.main()
