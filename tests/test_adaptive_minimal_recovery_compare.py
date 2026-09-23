"""Completed real recovery ledger, without fabricated molecular output."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import read_json, verify
from accommodation_fold_proposals import summarize_rows


class RecoveredLedger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = read_json(ROOT / 'workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_v1.json')
        cls.base = read_json(verify(cls.result['base_comparison']))

    def test_exact_population_and_reference(self):
        d = self.result
        self.assertEqual((len(d['rows']), len(d['pools']), len(d['triples'])), (225, 75, 100))
        self.assertEqual(len({r['case_id'] for r in d['rows']}), 225)
        self.assertEqual(d['reference'], self.base['reference'])
        self.assertFalse(d['thresholds_refitted'])

    def test_every_historical_method_field_is_unchanged(self):
        for kind in ('rows', 'pools', 'triples'):
            for before, after in zip(self.base[kind], self.result[kind]):
                for method, value in before['methods'].items():
                    self.assertEqual(after['methods'][method], value)

    def test_only_two_declared_unavailable_rows_gain_actual_pools(self):
        changed = []
        for row in self.result['rows']:
            old = row['methods']['minimal_operational']
            new = row['methods']['minimal_recovered']
            if row['recovery_result']:
                changed.append(row['case_id'])
                self.assertIsNone(old['R'])
                self.assertEqual(new['outcome'], 'correct')
                self.assertEqual(row['recovery_result']['case']['pool']['status'], 'available')
            else:
                self.assertEqual(new, old)
        self.assertEqual(set(changed), set(self.result['changed_source_ids']))
        self.assertEqual(len(changed), 2)

    def test_failed_source_and_strict_groups_remain_unavailable(self):
        index = {r['case_id']: r for r in self.result['rows']}
        failed = index['q60ar6-pqq-la_model__conditioned_La__seed-1_sample-0']
        self.assertIsNone(failed['methods']['minimal_recovered']['R'])
        for group in self.result['triples']:
            if any(index[c]['methods']['minimal_recovered']['R'] is None for c in group['members']):
                self.assertIsNone(group['methods']['minimal_recovered']['R'])

    def test_counts_come_from_all_actual_rows(self):
        d = self.result
        self.assertEqual(summarize_rows(d['rows'], d['bands']), d['counts']['all225'])
        self.assertEqual(d['counts']['all225']['minimal_recovered'],
                         {'denominator':225, 'correct':204, 'wrong':1, 'inconclusive':1, 'unavailable':19})
        self.assertEqual(d['new_molecular_calls_in_comparison'], 0)


if __name__ == '__main__':
    unittest.main()
