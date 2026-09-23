"""Run only after the actual full100 final comparison exists."""
from pathlib import Path
import sys
import statistics
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome
class EnvelopeFinalTransfer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.d=read_json(ROOT/'workspaces/motion_envelope_transfer_20260923/COMPARISON_v1.json')
    def test_full_ledger_and_original_correlated_groups(self):
        d=self.d;self.assertEqual(len(d['rows']),104);self.assertEqual(len(d['triples']),100)
        self.assertEqual(len(d['source_context_occurrences']),300);self.assertEqual(len(d['groups']),25)
        self.assertEqual(len({r['case_id']for r in d['rows']}),98)
        self.assertEqual(sum(t['status']!='prepared'for t in d['triples']),6)
        self.assertTrue(all(len(g['triple_ids'])==4 for g in d['groups']))
        self.assertFalse(d['new_thresholds_fitted']);self.assertEqual(sum(r['pool_reused']for r in d['rows']),3)
    def test_all_strict_medians_and_frozen_decisions(self):
        d=self.d
        for t in d['triples']:
            members=[r for r in d['source_context_occurrences']if r['triple_id']==t['triple_id']]
            self.assertEqual([r['case_id']for r in members],t['members'])
            for method,cell in t['methods'].items():
                values=[r['methods'][method]['R']for r in members]
                expected=None if any(v is None for v in values)else statistics.median(values)
                self.assertEqual(cell['R'],expected)
                self.assertEqual(cell['decision'],decision(expected,d['bands'][method]))
                self.assertEqual(cell['outcome'],outcome(cell['decision'],t['expected_class']))
    def test_actual_old_context_join_per_occurrence(self):
        d=self.d;old=read_json(verify(d['old_threefold_baseline']));triples={t['triple_id']:t for t in old['triples']}
        pairs={(r['selection_id'],r['case_id']):r for r in old['rows']}
        for occurrence in d['source_context_occurrences']:
            prior=triples[occurrence['triple_id']];key=(prior['selection_id'],occurrence['case_id'])
            self.assertEqual(occurrence['old_threefold_selection_id'],prior['selection_id'])
            if key in pairs:self.assertEqual(occurrence['methods']['threefold3p5'],pairs[key]['methods']['triple_union'])
            else:self.assertIsNone(occurrence['methods']['threefold3p5']['R'])
    def test_real_strict_cells_and_accommodation_sign(self):
        for r in self.d['rows']:
            if r['status']!='available':continue
            delta=r['methods']['envelope_operational']['R']-r['methods']['envelope_origin']['R']
            self.assertAlmostEqual(delta,r['selected_work']['Ca']['composite_kcal_mol']-r['selected_work']['La']['composite_kcal_mol'],places=6)
            for cells in r['matrix'].values():
                for cell in cells.values():
                    self.assertEqual(cell['status'],'complete')
                    for low in cell['low'].values():self.assertEqual(low['observed_TolE_hartree'],1e-10)
if __name__=='__main__':unittest.main()
