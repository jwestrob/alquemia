"""Actual archived matrices and staged inputs; no fabricated scientific outputs."""
import copy,json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import strict_native_transfer as run
from affordable_common import read_json,verify,InvalidArtifact
from precision_pool_continuation import data
from nikasha_pool import choose_rows

BASE=ROOT/'workspaces/strict_native_transfer_20260923/run_v1'
REF=ROOT/'workspaces/strict_native_pool_20260923/run_v1/REFERENCES.json'
QUAL=ROOT/'workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json'

class RealFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref=run.checked_reference(REF,QUAL)
        cls.fresh=data(cls.ref['collections']['fresh'])
        cls.low={run.tid(r):r for r in cls.fresh['rows']}
        cls.inv=data(data(cls.fresh['manifest'])['inventory'])
    def test_frozen_reference_and_eight_disjoint_reuses(self):
        self.assertEqual(len(run.REUSED),8)
        self.assertEqual(len(set(run.REUSED)),8)
        self.assertEqual(self.ref['reference_id'],'strict_native_fresh_canonical25_v1')
        self.assertFalse(self.ref['crystals_or_noncanonical_used_for_fit'])
    def test_real_strict_matrix_replays_qualified_four_folds(self):
        qualified=read_json(QUAL)
        for cid in run.FOLDS:
            group=next(x for x in self.inv['cases'] if x['case_id']==cid)
            old=next(x for x in data(group['collection'])['cases'] if x['case_id']==cid)
            mat=run.matrix_from_rows(old,self.low)
            actual=next(x for x in qualified['rows'] if x['case_id']==cid)['branches']['fresh']['pool']
            self.assertEqual(choose_rows(mat,list(run.CANDIDATES)),actual)
            for z in ('Ca','La'):
                for q in run.CANDIDATES:self.assertEqual(mat[z][q]['MACE'],old['matrix'][z][q]['MACE'])
    def test_corrupted_real_missing_cell_has_no_old_energy_fallback(self):
        cid=run.FOLDS[0];group=next(x for x in self.inv['cases'] if x['case_id']==cid)
        old=next(x for x in data(group['collection'])['cases'] if x['case_id']==cid)
        corrupted=copy.deepcopy(self.low);key='__'.join((cid,'origin','Ca','vacuum'))
        corrupted[key]['status']='explicit_corrupted_copy_for_missing_parser_test'
        corrupted[key]['energy_hartree']=None
        mat=run.matrix_from_rows(old,corrupted)
        self.assertIsNone(mat['Ca']['origin']['components'])
        self.assertEqual(choose_rows(mat,list(run.CANDIDATES))['status'],'unavailable')
    def test_actual_finite_disjoint_tasks_and_fresh_input(self):
        ip=BASE/'INVENTORY.json'
        if not ip.exists():self.skipTest('real finite preparation still running; not scientific success')
        inv=read_json(ip);seen=set()
        self.assertEqual(len(inv['rows']),2496)
        self.assertEqual(len(inv['cases']),225)
        for i in range(4):
            mp=BASE/f'shard_{i}'/'manifest.json'
            if not mp.exists():self.skipTest('real staged shards not complete')
            m=read_json(mp);ids={t['task_id'] for t in m['tasks']}
            self.assertEqual(len(ids),600);self.assertFalse(seen&ids);seen|=ids
            self.assertTrue(set(m['case_ids']).isdisjoint(run.REUSED))
            for t in m['tasks']:
                self.assertEqual(verify(t['input']).read_text(),run.recipe(t['charge'],t['medium'],'fresh'))
                self.assertEqual(t['xyz']['sha256'],t['source']['xyz']['sha256'])
                self.assertIsNone(t['seed_source']);self.assertFalse(t['gradient_requested'])
        self.assertEqual(len(seen),2400)
        self.assertEqual(sum(c['status']=='preparation_unavailable' for c in inv['cases']),17)

class ActualTransfer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        dest=BASE/'COMPARISON.json'
        if not dest.exists():raise unittest.SkipTest('actual full225 terminal comparison unavailable; no fabricated results')
        cls.final=read_json(dest);cls.prior=data(cls.final['prior'])
    def test_old_methods_memberships_and_failures_unchanged(self):
        self.assertEqual(len(self.final['rows']),225)
        self.assertEqual([r['case_id'] for r in self.final['rows']],[r['case_id'] for r in self.prior['rows']])
        for old,new in zip(self.prior['rows'],self.final['rows']):
            for method,value in old['methods'].items():self.assertEqual(new['methods'][method],value)
            if old['precision_status']=='preparation_unavailable':
                self.assertEqual(new['strict_status'],'prior_preparation_unavailable')
                self.assertIsNone(new['methods']['union_precision_strict']['R'])
        for section,n in [('pools',75),('triples',100)]:
            self.assertEqual(len(self.final[section]),n)
            for old,new in zip(self.prior[section],self.final[section]):
                self.assertEqual(new['members'],old['members'])
                for method,value in old['methods'].items():self.assertEqual(new['methods'][method],value)
    def test_actual_pool_algebra_and_frozen_decisions(self):
        ref=data(self.final['reference'])['branches']['fresh']
        for case in self.final['actual_pools']:
            self.assertEqual(run.choose_rows(case['matrix'],list(run.CANDIDATES)),case['pool'])
            row=next(r for r in self.final['rows'] if r['case_id']==case['case_id'])
            for method,v in zip(run.METHODS,run.VARIANTS):
                val=case['pool'][v]['composite_R_model_kcal_mol'] if case['pool']['status']=='available' else None
                self.assertEqual(row['methods'][method]['R'],val)
                self.assertEqual(row['methods'][method]['decision'],run.decision(val,ref['variants'][v]['bands']))
        self.assertFalse(self.final['new_thresholds_fitted'])
        self.assertEqual(self.final['cell_denominator'],2496)

if __name__=='__main__':unittest.main()
