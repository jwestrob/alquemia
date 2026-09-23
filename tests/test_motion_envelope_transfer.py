"""Actual full100 ledger and strict finite-stage preparation checks."""
import copy
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
import motion_envelope_transfer as run
import motion_envelope_transfer_inventory as inv
from motion_envelope_transfer_compare import aggregate
BASE=ROOT/'workspaces/motion_envelope_transfer_20260923/run_v1'
class EnvelopeTransfer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.a=read_json(BASE/'INVENTORY.json')
    def test_all100_and_exact_three_reuses(self):
        a=self.a;self.assertEqual(len(a['triples']),100);self.assertEqual(sum(t['status']=='prepared'for t in a['triples']),94)
        self.assertEqual(len(a['rows']),104);reuse=[r for r in a['rows']if r['reuse']]
        self.assertEqual(len(reuse),3);self.assertEqual(sum(r['reuse']['scalar_cells_verified']for r in reuse),36)
        self.assertEqual({r['case_id']for r in reuse},{'a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4','a8r3s4-pqq-la_model__conditioned_La__seed-1_sample-1','a8r3s4-pqq-la_model__conditioned_La__seed-1_sample-3'})
        self.assertTrue(all(t['missing_members']for t in a['triples']if t['status']!='prepared'))
    def test_disjoint_finite_shards(self):
        shards=[read_json(verify(p))for p in self.a['shards']];self.assertEqual([len(s['cases'])for s in shards],[51,50])
        keys=[inv.key(r)for s in shards for r in s['cases']];self.assertEqual(len(keys),len(set(keys)));self.assertEqual(len(keys),101)
        self.assertEqual(sum(len(s['origin_tasks'])for s in shards),202);self.assertEqual(self.a['counts']['maximum_new_GFN2_total'],1212)
        for i in (0,1):run.scope(BASE/'INVENTORY.json',i)
    def test_changed_prepared_state_cannot_reuse(self):
        a=self.a;d=read_json(verify(a['pilot_collection']));pm=read_json(verify(d['manifest']));sm=read_json(verify(pm['source_manifest']));ref=read_json(verify(a['reference']))
        r=next(r for r in a['rows']if r['reuse']);c=next(c for c in d['cases']if c['case_id']==r['reuse']['case_id'])
        broken=copy.deepcopy(r['prepared']);broken['state_key']='corrupted-real-artifact'
        with self.assertRaises(InvalidArtifact):inv.audit_reuse(broken,c,verify(a['pilot_collection']),pm,sm,ref)
    def test_actual_triple_median_and_missing_are_strict(self):
        old=read_json(ROOT/'workspaces/union_triple_transfer_20260923/COMPARISON_v1.json')
        t=next(t for t in old['triples']if t['methods']['triple_union']['R']is not None)
        cells=[{'R':v}for v in t['methods']['triple_union']['member_R']]
        got=aggregate(cells,old['bands']['triple_union'],t['expected_class'])
        self.assertEqual(got['R'],t['methods']['triple_union']['R'])
        cells[1]['R']=None
        self.assertIsNone(aggregate(cells,old['bands']['triple_union'],t['expected_class'])['R'])

    def test_actual_rejected_runner_path_is_exact_byte_copy(self):
        mp=BASE/'shard_0/origins/scalar/manifest.json';m=read_json(mp);t=m['tasks'][0]
        receipt=read_json(t['output_path']+'.execution.json')
        self.assertNotEqual(receipt['task_runner']['path'],m['execution_policy']['task_runner']['path'])
        self.assertEqual(receipt['task_runner']['sha256'],m['execution_policy']['task_runner']['sha256'])
        self.assertEqual(receipt['runtime_renderer']['sha256'],m['execution_policy']['runtime_renderer']['sha256'])
        parsed=run.scalar_engine().endpoint(mp,t)
        self.assertEqual(parsed['status'],'complete');self.assertEqual(parsed['observed_TolE_hartree'],1e-10)
        self.assertEqual(parsed['receipt']['path'],t['output_path']+'.execution.json')

if __name__=='__main__':unittest.main()
