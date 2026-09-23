"""Actual final34 matrices, canonical-only references and matched component algebra."""
from pathlib import Path
from collections import Counter
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from nikasha_pool import choose_rows
from nikasha_pool_compare import extrema_reference
from accommodation_nonlinear import contrast_components,relative_components
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome
BASE=ROOT/'workspaces/motion_envelope_20260923'
class EnvelopeFinal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection=read_json(BASE/'pilot34_pool_v1/COLLECTION_FINAL.json')
        cls.reference=read_json(BASE/'REFERENCE_PINNED_v1.json');cls.comparison=read_json(BASE/'COMPARISON_v1.json')
    def test_actual34_complete_matrices_and408_scalar_cells(self):
        self.assertEqual(self.collection['denominator'],34);self.assertEqual(self.collection['available'],34)
        self.assertEqual(self.collection['GFN2_complete'],272)
        self.assertEqual(read_json(BASE/'pilot34_scalar_origins_v1/COLLECTION.json')['complete'],136)
        for c in self.collection['cases']:
            self.assertEqual(c['pool'],choose_rows(c['matrix'],[x['id']for x in c['candidates']]))
            self.assertEqual(set(c['matrix']['Ca']),{'origin','adaptive_Ca','adaptive_La'})
            for cells in c['matrix'].values():
                for cell in cells.values():
                    self.assertEqual(cell['status'],'complete')
                    for low in cell['low'].values():self.assertEqual(low['observed_TolE_hartree'],1e-10)
    def test_canonical_only_bands_and_every_declared_role(self):
        rows=self.reference['rows'];roles=Counter(r['role']for r in rows)
        self.assertEqual(roles,{'canonical_calibration':25,'crystal_transfer':3,'separate_probe':6})
        canon=[r for r in rows if r['role']=='canonical_calibration']
        self.assertEqual(len({r['root_case_id']for r in canon}),25)
        for v in ('mathematical','operational'):
            expected=extrema_reference([{**r,'R_model_kcal_mol':r['scores'][v]}for r in canon],v)
            self.assertEqual(expected['bands'],self.reference['variants'][v]['bands'])
            for r in rows:
                bands=expected['bands'];d=decision(r['scores'][v],bands)if bands else'unavailable'
                self.assertEqual(r['own_reference'][v],{'decision':d,'outcome':outcome(d,r['expected_class'])})
        static=extrema_reference([{**r,'R_model_kcal_mol':r['origin_R_model_kcal_mol']}for r in canon],'origin')
        self.assertEqual(static['bands'],self.reference['origin_only_reference']['bands'])
    def test_matched_physical_identity_and_origin_work_decomposition(self):
        self.assertEqual(self.comparison['denominator'],34)
        for r in self.comparison['rows']:
            self.assertEqual(r['same_original_core_max_A'],{'Ca':0.,'La':0.})
            for variant,v in r['variants'].items():
                for k,total in v['pool_shift'].items():
                    self.assertAlmostEqual(total,v['origin_shift'][k]+v['differential_accommodation_shift'][k],places=8)
                for label,mx,pool in [('envelope',r['matrix'],r['pool']),('strict_tenfold',r['strict_tenfold_matrix'],r['strict_tenfold_pool'])]:
                    origin=contrast_components(mx['Ca']['origin']['components'],mx['La']['origin']['components'])
                    work=v['selected_endpoint_work'][label]
                    self.assertAlmostEqual(pool[variant]['composite_R_model_kcal_mol']-origin['composite_R_model_kcal_mol'],work['Ca']['composite_kcal_mol']-work['La']['composite_kcal_mol'],places=7)
if __name__=='__main__':unittest.main()
