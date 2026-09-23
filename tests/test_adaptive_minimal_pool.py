"""Real archived adaptive matrices and frozen-reference replay, no model calls."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import adaptive_minimal_pool as m


class MinimalAdaptive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.canonical = m.read_json(ROOT / 'workspaces/adaptive_completion_20260922/original30_pool_v1/final_collection.json')
        cls.reference_path = ROOT / 'workspaces/adaptive_minimal_pool_20260923/REFERENCE_v1.json'
        cls.reference = m.read_json(cls.reference_path)
        cls.result = m.read_json(ROOT / 'workspaces/adaptive_minimal_pool_20260923/COMPARISON_v1.json')

    def test_actual_canonical_subset_and_minimum_algebra(self):
        for case in self.canonical['cases']:
            reduced = m.reduced_case(case)
            self.assertEqual(reduced['pool']['status'], 'available')
            for metal in ('Ca', 'La'):
                self.assertEqual(set(reduced['pool']['rows'][metal]['work_from_origin_kcal_mol']), set(m.CANDIDATES))
                self.assertIn(reduced['pool']['rows'][metal]['operational_candidate'], m.CANDIDATES)

    def test_corrupted_real_required_cell_stays_unavailable(self):
        case = copy.deepcopy(self.canonical['cases'][0])
        del case['matrix']['La']['adaptive_Ca']
        result = m.reduced_case(case)['pool']
        self.assertEqual(result['status'], 'unavailable')
        self.assertIsNone(result['operational'])

    def test_labels_never_enter_geometry_selection(self):
        case = copy.deepcopy(self.canonical['cases'][0])
        original = m.reduced_case(case)['pool']
        case['old_result']['expected_class'] = 'deliberately_corrupted_label'
        self.assertEqual(m.reduced_case(case)['pool'], original)

    def test_reference_replay_and_reject_corrupted_canonical_label(self):
        signature = {'signature': self.reference['signature']}
        self.assertEqual(m.checked_reference(self.reference_path, signature), self.reference)
        altered = copy.deepcopy(self.reference)
        altered['variants']['operational']['rows'][0]['expected_class'] = 'corrupted'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'corrupted_reference.json'
            m.write_new(path, altered)
            with self.assertRaises(m.InvalidArtifact):
                m.checked_reference(path, signature)

    def test_complete_transfer_denominators_and_strict_missing_groups(self):
        d = self.result
        self.assertEqual((len(d['rows']), len(d['pools']), len(d['triples'])), (225, 75, 100))
        index = {r['case_id']: r for r in d['rows']}
        for group in d['triples']:
            values = [index[c]['methods']['minimal_operational']['R'] for c in group['members']]
            if any(v is None for v in values):
                self.assertIsNone(group['methods']['minimal_operational']['R'])
        for case in d['rows']:
            pool = case['minimal_pool']
            if pool['status'] != 'available':
                self.assertIsNone(case['methods']['minimal_operational']['R'])

    def test_saved_actual_counts_and_matching_decisions(self):
        d = self.result
        counts = m.summarize_rows(d['rows'], d['bands'])
        self.assertEqual(counts, d['counts']['all225'])
        for case in d['rows']:
            # This tests the observed completed replay, not an assumed scientific fixture.
            self.assertEqual(case['methods']['minimal_operational']['decision'],
                             case['methods']['full_adaptive']['decision'])
        self.assertEqual(d['new_molecular_calls'], 0)
        self.assertFalse(d['new_numerical_qualification'])


if __name__ == '__main__':
    unittest.main()
