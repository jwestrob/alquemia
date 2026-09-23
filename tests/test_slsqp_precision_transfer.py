"""Pinned source/parallelism checks; no molecular calls or fabricated outputs."""
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import slsqp_precision_transfer as transfer
import adaptive_completion
from affordable_common import read_json,verify,record,write_new
from compact_solvation import input_text
from affordable_workflow import dry_run

class TransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=ROOT/'workspaces/slsqp_precision_transfer_20260923/run_v1'
        cls.top=read_json(cls.root/'SELECTION.json')
        cls.ready=read_json(cls.root/'READY.json')
        cls.manifests=[verify(p) for p in cls.ready['shards']]

    def test_all_four_real_source_preflights(self):
        totals=[read_json(p.parent/'PREFLIGHT.json') for p in self.manifests]
        for path,receipt in zip(self.manifests,totals):
            self.assertEqual(receipt['status'],'validated');self.assertEqual(receipt['manifest'],record(path))
        self.assertEqual([r['sources'] for r in totals],[51,51,50,50])
        self.assertEqual(sum(r['new_searches'] for r in totals),404)
        self.assertEqual(sum(r['GFN2_maximum'] for r in totals),1616)

    def test_full_denominator_disjointness_and_actual_six_reuses(self):
        rows=self.top['all_cases'];self.assertEqual(len(rows),225)
        fresh=[c for p in self.manifests for c in read_json(p)['declared_case_ids']]
        reused=self.top['reuse_case_ids']
        self.assertEqual((len(fresh),len(set(fresh)),len(reused)),(202,202,6))
        self.assertFalse(set(fresh)&set(reused))
        self.assertEqual(set(fresh)|set(reused),{r['case_id'] for r in rows if r['union_origin_status']=='complete'})
        self.assertEqual(sum(r['union_origin_status']!='complete' for r in rows),17)
        repaired=next(r for r in rows if r['case_id']=='mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-4')
        self.assertEqual(repaired['union_origin_status'],'preparation_unavailable')
        self.assertIsNone(repaired['precision34_reuse'])
        for row in rows:
            if not row['precision34_reuse']:continue
            coll=read_json(verify(row['precision34_reuse']));c=next(c for c in coll['cases'] if c['case_id']==row['case_id'])
            self.assertEqual(c['pool']['status'],'available')
            for v,R in row['reused_score'].items():self.assertEqual(c['pool'][v]['composite_R_model_kcal_mol'],R)

    def test_frozen_reference_and_actual_rank_qualification(self):
        self.assertEqual(self.top['reference']['sha256'],transfer.REFERENCE_HASH)
        self.assertEqual(transfer.qualification(verify(self.top['rank_qualification'])),self.top['rank_qualification'])
        ref=read_json(verify(self.top['reference']));self.assertEqual(ref['calibration_denominator'],25)
        self.assertFalse(ref['crystals_or_noncanonical_used_for_fit'])
        inv=read_json(verify(self.top['inventory']));old=read_json(verify(inv['precision34_comparison']))
        self.assertEqual((old['numerical_passes'],old['numerical_comparison_denominator']),(27,32))

    def test_private_numerical_policy_does_not_mutate_shared(self):
        old=dict(adaptive_completion.SETTINGS);u,p=transfer.engine()
        self.assertEqual({k for k in old if old[k]!=u.SETTINGS[k]},{'optimizer_ftol'})
        self.assertEqual(u.SETTINGS['optimizer_ftol'],1e-8)
        self.assertEqual(adaptive_completion.SETTINGS,old)
        self.assertIs(p.low_prepare,transfer.low_prepare)

    def test_actual_input_recipe_rank1_and_shard_containment(self):
        """Stage one real coordinate fixture; no computed energy is supplied."""
        source=read_json(self.manifests[0]);t=source['tasks'][0]
        cell={k:t[k] for k in ('task_id','case_id','metal','xyz','charge','multiplicity')}
        m={'protocol_id':transfer.pilot.POOL_PROTOCOL,'source_manifest':record(self.manifests[0]),
           'agreement':self.top['agreement'],'orca':source['orca'],'implementation':source['implementation'],
           'shard_count':1,'tasks':[cell]}
        with tempfile.TemporaryDirectory() as directory:
            mp=Path(directory)/'manifest.json';write_new(mp,m);transfer.low_prepare(mp)
            lp=Path(directory)/'solvent/shard_0/manifest.json';low=read_json(lp)
            self.assertEqual(low['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
            self.assertEqual(low['cached_origin_ranks'],8)
            for row in low['tasks']:
                expected=input_text(row['charge'],row['multiplicity'],row['medium'],'native').replace('%scf\n','%scf\n MaxIter 500\n')
                self.assertEqual(verify(row['input']).read_text(),expected)
                for key in ('input','xyz'):self.assertTrue(verify(row[key]).is_relative_to(lp.parent))
                self.assertTrue(Path(row['output_path']).is_relative_to(lp.parent))
            self.assertEqual(dry_run(lp)['tasks'],2)

    def test_actual_reused_numerical_comparison_replay(self):
        from slsqp_precision_transfer_report import numerical
        ref=read_json(ROOT/'workspaces/slsqp_precision_20260923/COMPARISON_v1.json')
        coll=read_json(verify(ref['collection']));m=read_json(verify(coll['manifest']))
        sm=read_json(verify(m['source_manifest']))
        expected=next(r for r in ref['rows'] if r['case_id']=='1H4I')
        src=next(r for r in sm['sources'] if r['case_id']=='1H4I')
        old=next(c for c in read_json(verify(src['old_pool']))['cases'] if c['case_id']=='1H4I')
        new=next(c for c in coll['cases'] if c['case_id']=='1H4I')
        actual=numerical(old,new,src['parent'],m['source_manifest'])
        self.assertEqual(actual['numerical_pass'],expected['numerical_pass'])
        self.assertEqual(actual['delta_R']['operational'],expected['delta_R'])
        self.assertEqual(actual['components'],expected['components'])

if __name__=='__main__':unittest.main()
