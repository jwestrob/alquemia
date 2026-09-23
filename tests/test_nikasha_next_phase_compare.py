"""Join real prior reports; corrupted-copy checks are not scientific evidence."""
import copy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json
from nikasha_next_phase_compare import join_rows, key


class ComparisonJoin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        a = ROOT / 'workspaces/nikasha_recovery_20260922/proposal_comparison.json'
        b = ROOT / 'workspaces/nikasha_shared_pool_20260922/primary225_v1/final_comparison.json'
        if not a.exists() or not b.exists(): raise unittest.SkipTest('actual completed comparison fixtures absent')
        cls.base, cls.pool = read_json(a), read_json(b)

    def test_real_rows_and_strict_groups_retain_membership_and_missing_scores(self):
        for kind, count in (('rows', 225), ('pools', 75), ('triples', 100)):
            result = join_rows(self.base[kind], self.pool[kind], kind,
                               'pool_operational_new', 'actual_prior_common_pool')
            self.assertEqual(len(result), count)
            saved = {key(r, kind): r for r in self.pool[kind]}
            for original, joined in zip(self.base[kind], result):
                self.assertEqual(joined['methods']['context_composite'], original['methods']['context_composite'])
                self.assertEqual(joined['methods']['actual_prior_common_pool'],
                                 saved[key(original, kind)]['methods']['pool_operational_new'])

    def test_corrupted_real_labels_or_denominators_cannot_join(self):
        damaged = copy.deepcopy(self.pool['rows'])
        damaged[0]['expected_class'] = 'CORRUPTED_LABEL'
        with self.assertRaisesRegex(InvalidArtifact, 'label changed'):
            join_rows(self.base['rows'], damaged, 'rows', 'pool_operational_new', 'invalid')
        with self.assertRaisesRegex(InvalidArtifact, 'populations differ'):
            join_rows(self.base['rows'], self.pool['rows'][:-1], 'rows', 'pool_operational_new', 'invalid')


if __name__ == '__main__': unittest.main()
