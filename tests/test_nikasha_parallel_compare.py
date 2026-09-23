"""Actual archived matrix comparisons, with explicitly corrupted copies only."""
import copy
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json
from nikasha_finite_candidates import sources
from nikasha_parallel_compare import result_row,pinned,choose_rows


class PilotComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=ROOT/'diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json'
        if not p.exists():raise unittest.SkipTest('actual pilot pins unavailable')
        cls.inputs=read_json(p);cls.reference=pinned(cls.inputs['adaptive_reference'])

    def test_actual_identity_pool_preserves_both_band_policies_and_zero_work(self):
        for src in self.inputs['cases']:
            _,base,_,_=sources(self.inputs,src['case_id'])
            row=result_row({**base,'prior_pool':base['pool']},base,src['known_class'],self.reference)
            self.assertEqual(row['delta_R_from_adaptive'],0.)
            for z in ('Ca','La'):self.assertEqual(row['endpoint_work_from_adaptive'][z]['composite_kcal_mol'],0.)
            for b in ('released','adaptive'):self.assertEqual(row['decisions'][b]['candidate'],row['decisions'][b]['adaptive'])

    def test_missing_actual_required_cell_is_unavailable(self):
        src=self.inputs['cases'][0];_,base,_,_=sources(self.inputs,src['case_id'])
        case=copy.deepcopy(base);case['prior_pool']=base['pool']
        case['matrix']['Ca']['origin']['status']='unavailable'
        case['pool']=choose_rows(case['matrix'],[q['id'] for q in case['candidates']])
        row=result_row(case,base,src['known_class'],self.reference)
        self.assertIsNone(row['candidate_R']);self.assertIsNotNone(row['adaptive_R'])
        self.assertEqual(row['decisions']['adaptive']['candidate']['outcome'],'unavailable')

    def test_corrupted_archived_score_rejected(self):
        src=self.inputs['cases'][0];_,base,_,_=sources(self.inputs,src['case_id'])
        case=copy.deepcopy(base);case['prior_pool']=base['pool']
        case['pool']['operational']['composite_R_model_kcal_mol']+=1.
        with self.assertRaises(InvalidArtifact):result_row(case,base,src['known_class'],self.reference)


if __name__=='__main__':unittest.main()
