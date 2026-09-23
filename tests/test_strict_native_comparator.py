import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import strict_native_comparator as s
from affordable_common import read_json,verify,xyz
D=ROOT/'workspaces/strict_native_comparator_20260923/run_v1'
class ComparatorTests(unittest.TestCase):
    def test_real_scope_recipe_and_paired_coordinates(self):
        self.assertEqual(s.validate(D/'manifest.json')['fresh_calls'],48);m=read_json(D/'manifest.json')
        for cid in s.CASES:
            for q in s.CANDIDATES:
                rr={t['metal']:t for t in m['tasks'] if t['case_id']==cid and t['candidate']==q and t['medium']=='vacuum'}
                ca,la=rr['Ca'],rr['La'];a=xyz(verify(ca['xyz']));b=xyz(verify(la['xyz']))
                self.assertEqual(la['charge']-ca['charge'],1);self.assertEqual([x[0] for x in a[1:]],[x[0] for x in b[1:]])
                self.assertLessEqual(np.max(np.abs(np.asarray([x[1:] for x in a])-np.asarray([x[1:] for x in b]))),1e-12)
                for t in rr.values():self.assertIn('NoAutostart',verify(t['input']).read_text());self.assertIsNone(t['seed_source'])
    def test_actual_collection_no_new_reference_or_fallback(self):
        if not (D/'COLLECTION.json').exists():self.skipTest('48 molecular calls unrun')
        c=read_json(D/'COLLECTION.json');self.assertEqual((len(c['rows']),len(c['cases'])),(48,4));self.assertFalse(c['new_calibration'])
        from nikasha_pool import choose_rows
        from accommodation_folds_compare import decision
        ref=read_json(verify(c['reference']))['branches']['fresh']
        for case in c['cases']:
            self.assertEqual(choose_rows(case['matrix'],list(s.CANDIDATES)),case['pool'])
            for v in ('mathematical','operational'):
                val=case['pool'][v]['composite_R_model_kcal_mol'] if case['status']=='available' else None
                self.assertEqual(case['calls'][v]['decision'],decision(val,ref['variants'][v]['bands']))
        for r in c['rows']:
            if r['status']=='complete':self.assertEqual(r['observed_TolE_hartree'],1e-10);self.assertEqual(r['initial_guess'],['SAD'])
            else:self.assertIsNone(r['energy_hartree'])
if __name__=='__main__':unittest.main()
