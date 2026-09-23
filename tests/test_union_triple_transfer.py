"""Real declared triple sources and actual archived pools; no synthetic energies."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new
import union_triple_transfer as transfer
import nikasha_pool as pool

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'workspaces/union_triple_transfer_20260923/run_v2'


class TripleTransfer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.a=read_json(BASE/'AUDIT.json')

    def test_declared100_and125_preserve_missing_and_exclude_stress(self):
        a=self.a;self.assertEqual(len(a['triples']),100);self.assertEqual(len(a['rows']),125)
        self.assertEqual(sum(t['status']=='prepared' for t in a['triples']),94)
        self.assertEqual(sum(t['status']!='prepared' for t in a['triples']),6)
        self.assertEqual(len({r['pair_id'] for r in a['rows']}),125)
        self.assertTrue(all(not r['source']['is_separate_stress_probe'] for r in a['rows']))
        self.assertTrue(all(r['source']['source']['source_conditioning_metal']=='La' for r in a['rows']))
        self.assertTrue(all(not r['source']['source']['canonical_coordinate_match'] for r in a['rows']))
        keys={(r['selection_id'],r['case_id']) for r in a['rows']}
        for t in a['triples']:
            if t['status']=='prepared':self.assertTrue(all((t['selection_id'],cid) in keys for cid in t['members']))

    def test_actual70_pools_and_new55_scope(self):
        self.assertEqual(self.a['counts'],{'triple_denominator':100,'complete_prepared_triples':94,'old_unavailable_triples':6,
            'source_membership_pairs':125,'full_pool_reuses':70,'new_source_pools':55,'new_origin_MACE':38,
            'new_origin_GFN2':76,'new_searches':110,'maximum_cross_MACE':110,'maximum_candidate_GFN2':440})
        for r in self.a['rows']:
            p=r['pool_reuse']
            if p:
                self.assertEqual(pool.choose_rows(p['pool']['matrix'],transfer.CANDIDATES),p['pool']['pool'])
                self.assertTrue(all(x<=1e-12 for x in p['maximum_map_coordinate_copy_difference_A'].values()))
                m=read_json(verify(p['source_manifest']));self.assertEqual(m['settings']['optimizer_ftol'],1e-8)

    def test_actual_coordinate_copy_difference_and_corrupted_graph(self):
        r=next(r for r in self.a['rows'] if r['pool_reuse'] and any(r['pool_reuse']['maximum_map_coordinate_copy_difference_A'].values()))
        original=read_json(verify(r['pool_reuse']['mapping']['Ca']));mapped=transfer.mapped(r['source'],self.a['config'])['Ca']
        self.assertLessEqual(transfer.mapping_equivalence(mapped,original),1e-12)
        corrupt=copy.deepcopy(original);corrupt['context']['core_positions_A'][0][0]+=.001
        with self.assertRaisesRegex(InvalidArtifact,'coordinates differ'):transfer.mapping_equivalence(corrupt,original)
        corrupt=copy.deepcopy(original);corrupt['context']['physical_ids'][0]='explicitly_corrupted_source_identity'
        with self.assertRaisesRegex(InvalidArtifact,'graph, cap rule or mode'):transfer.mapping_equivalence(corrupt,original)

    def test_actual76_origin_tasks_validate_and_charge_corruption_rejected(self):
        path=BASE/'origins/solvent/manifest.json';transfer.validate_origins(path);m=read_json(path)
        self.assertEqual(len(m['tasks']),76)
        self.assertEqual(m['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
        corrupt=copy.deepcopy(m);corrupt['tasks'][0]['charge']+=1;corrupt['all_tasks']=corrupt['tasks']
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_actual_charge.json';write_new(p,corrupt)
            with self.assertRaisesRegex(InvalidArtifact,'source/membership/state'):transfer.validate_origins(p)

    def test_fixed_finite_search_scope_has_no_guessed_forces(self):
        s=read_json(BASE/'MOLECULAR_SCOPE.json')
        self.assertEqual((len(s['search_tasks']),len(s['cross_MACE_requests']),len(s['candidate_GFN2_requests'])),(110,110,440))
        self.assertEqual((len(s['fully_preflighted_force_pairs']),len(s['force_pairs_waiting_declared_origins'])),(36,19))
        for t in s['search_tasks']:
            if t['origin_status']=='await_finite_origin_tasks':self.assertNotIn('active_indices',t)
            else:self.assertEqual(len(t['active_indices']),4)
        self.assertEqual(s['new_molecular_calls'],0);self.assertFalse(s['production_changed'])

    def test_actual_complete_origins_have38_native_and76_scalar_calls(self):
        d=read_json(BASE/'origins/COLLECTION.json')
        self.assertEqual((d['denominator'],d['complete']),(55,55))
        self.assertEqual(sum(not e['reused'] for r in d['rows'] for e in r['native_endpoints'].values()),38)
        self.assertEqual(sum(not e['reused'] for r in d['rows'] for z in r['solvent_endpoints'].values() for e in z.values()),76)
        self.assertTrue(all(e['component_audit']['charge_sanity_status']=='pass' for r in d['rows'] for z in r['solvent_endpoints'].values() for e in z.values()))

    def test_actual_two_shards_exactly_cover_frozen55_and110(self):
        ready=read_json(BASE/'searches/READY.json');self.assertEqual(ready['sources_per_shard'],[28,27]);ids=[];tasks=[]
        for pin in ready['shards']:
            m=read_json(verify(pin));verify(m['execution_adapter'])
            self.assertEqual(m['reference'],self.a['reference']);self.assertEqual(m['settings'],self.a['settings'])
            self.assertEqual(m['unavailable'],[]);ids.extend(m['declared_case_ids']);tasks.extend(m['tasks'])
            self.assertEqual([c['case_id'] for c in m['cases']],m['declared_source_pairs'])
        self.assertEqual(set(ids),{r['pair_id'] for r in self.a['rows'] if not r['pool_reuse']})
        self.assertEqual(len(ids),55);self.assertEqual(len({t['task_id'] for t in tasks}),110)


if __name__=='__main__':unittest.main()
