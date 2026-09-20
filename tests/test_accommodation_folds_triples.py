"""Three-fold descriptor checks using the actual completed reference outputs."""
from itertools import combinations
from pathlib import Path
import statistics
import sys
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify


class ThreeFolds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=read_json(ROOT/'workspaces/accommodation_goal_20260920/folds_v1/triples_v1/result.json')
        cls.original=read_json(verify(cls.result['comparison']))

    def test_every_declared_triple_is_present_without_selection(self):
        self.assertEqual(len(self.result['rows']),200)
        for p in self.result['proteins']:
            expected={frozenset(x) for x in combinations(p['declared_La4_members'],3)}
            self.assertEqual({frozenset(t['members']) for t in p['triples']},expected)
            self.assertEqual(len(p['triples']),4)
            for t in p['triples']:self.assertNotIn(t['omitted_member'],t['members'])

    def test_actual_unrounded_medians_and_missing_values(self):
        original={(r['case_id'],r['representation']):r for r in self.original['rows']}
        for row in self.result['rows']:
            for method,field in [('native','native_R_model_kcal_mol'),('composite','composite_R_model_kcal_mol')]:
                values=[original[c,row['representation']][field] for c in row['members']]
                actual=row['methods'][method]
                if any(v is None for v in values):
                    self.assertIsNone(actual['median_R_model_kcal_mol'])
                    self.assertEqual(actual['outcome'],'unavailable')
                else:self.assertEqual(actual['median_R_model_kcal_mol'],statistics.median(values))

    def test_actual_context_coverage_and_worst_outcomes(self):
        context=self.result['counts']['context']
        self.assertEqual(context['composite']['triples'],{'correct':94,'wrong':0,'inconclusive':0,'unavailable':6})
        self.assertEqual(context['composite']['proteins_all_four_correct'],23)
        self.assertEqual(context['native']['triples']['wrong'],2)
        for root in ('mmol_1770-pqq-la_model','q88jh5-pqq-la_model'):
            p=next(p for p in self.result['proteins'] if p['representation']=='context' and p['root_case_id']==root)
            outcome=p['methods']['composite']
            self.assertEqual(outcome['worst_outcome'],'unavailable')
            self.assertEqual(outcome['worst_scored_outcome'],'correct')
            self.assertEqual(outcome['counts']['unavailable'],3)
            self.assertFalse(outcome['all_four_correct'])

    def test_old_bands_original_results_and_simple_comparator_retained(self):
        self.assertEqual(self.result['frozen_bands'],self.original['frozen_bands'])
        self.assertEqual(self.result['counts']['core']['native']['triples'],self.result['counts']['context']['composite']['triples'])
        self.assertEqual(self.result['new_molecular_calls'],0)
        self.assertFalse(self.result['production_changed'])
        self.assertFalse(self.result['new_PLM_validation'])
        self.assertFalse(self.result['threshold_refitted'])


if __name__=='__main__':unittest.main()
