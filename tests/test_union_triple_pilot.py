"""Finite pilot uses real reconstructed inputs and actual saved precision pools."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
import union_triple_origins as origin
import nikasha_pool as pool
import union_triple_adaptive as adaptive

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'workspaces/union_triple_pilot_20260923'


class TriplePilot(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reuse=read_json(BASE/'POOL_REUSE_v1.json');cls.i=read_json(BASE/'origins_v1/INPUTS.json')
        cls.low=read_json(BASE/'origins_v1/solvent/manifest.json');cls.mace=read_json(BASE/'origins_v1/mace/manifest.json')

    def test_exact19_new_precision_pools_and_9_changed(self):
        self.assertEqual((self.reuse['denominator'],self.reuse['complete_pool_reuses'],self.reuse['new_pools']),(28,19,9))
        for r in self.reuse['rows']:
            if r['status']!='exact_new_ftol_pool_reuse':continue
            c=r['pool'];self.assertEqual(pool.choose_rows(c['matrix'],[q['id'] for q in c['candidates']]),c['pool'])
            m=read_json(verify(r['source_manifest']));self.assertEqual(m['settings']['optimizer_ftol'],1e-8)
            for z in ('Ca','La'):verify(r['mapping'][z])

    def test_finite_actual_origin_scope_and_rank1(self):
        result=origin.validate(BASE/'origins_v1/solvent/manifest.json')
        self.assertEqual(len(self.low['tasks']),24);self.assertEqual(len(self.mace['tasks']),12)
        self.assertEqual(len(self.i['cases']),29);self.assertEqual(len(self.i['new_origin_sources']),6)
        self.assertEqual(self.low['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
        self.assertEqual({t['case_id'] for t in self.low['tasks']},set(self.i['new_origin_sources']))
        self.assertTrue(result)

    def test_actual_paired_origins_and_explicit_medium_switch(self):
        for cid in self.i['new_origin_sources']:
            ts=[t for t in self.low['tasks'] if t['case_id']==cid]
            self.assertEqual({(t['metal'],t['medium']) for t in ts},{(z,s) for z in ('Ca','La') for s in ('vacuum','alpb')})
            coords=[xyz(verify(t['xyz'])) for t in ts]
            self.assertTrue(all(a[1:]==coords[0][1:] for a in coords))
            for t in ts:
                text=verify(t['input']).read_text();self.assertIn('MaxIter 500',text);self.assertIn('UseXTBMixer true',text)
                self.assertEqual('ALPB(Water)' in text,t['medium']=='alpb');self.assertNotIn('EnGrad',text)

    def test_corrupted_actual_origin_source_rejected(self):
        d=copy.deepcopy(self.low);d['tasks'][0]['charge']+=1;d['all_tasks']=d['tasks']
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted_charge.json';write_new(path,d)
            with self.assertRaisesRegex(InvalidArtifact,'origin source/state'):origin.validate(path)

    def test_actual29_origin_pairs_complete(self):
        d=read_json(BASE/'origins_v1/COLLECTION.json')
        self.assertEqual((d['denominator'],d['complete']),(29,29))
        self.assertEqual(sum(not e['reused'] for r in d['rows'] for e in r['native_endpoints'].values()),12)
        self.assertEqual(sum(not e['reused'] for r in d['rows'] for z in r['solvent_endpoints'].values() for e in z.values()),24)

    def test_exact20_prepared_searches_and_private_profile(self):
        path=BASE/'proposals_v1/manifest.json';m=read_json(path);v=adaptive.validate(path)
        self.assertEqual((v['prepared_sources'],v['new_searches'],v['maximum_GFN2']),(10,20,80))
        self.assertEqual(m['settings']['optimizer_ftol'],1e-8);self.assertIsNone(m['reference'])
        u,p=adaptive.engine();self.assertEqual(u.PROTOCOL,adaptive.PROTOCOL)
        self.assertEqual(u.POOL_PROTOCOL,adaptive.POOL_PROTOCOL);self.assertEqual(u.SETTINGS,m['settings'])
        for t in m['tasks']:
            self.assertEqual(len(t['active_indices']),4)
            self.assertEqual(t['active_mode_ids'],[r['id'] for r in t['selector']['selected']])

    def test_corrupted_actual_proposal_bounds_profile_rejected(self):
        d=read_json(BASE/'proposals_v1/manifest.json');d['settings']['optimizer_ftol']=1e-3
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted_optimizer_tolerance.json';write_new(path,d)
            with self.assertRaisesRegex(InvalidArtifact,'scope/settings'):adaptive.validate(path)

    def test_final_reference_uses_exact25_not_crystals_or_stress(self):
        d=read_json(BASE/'REFERENCE_v1.json');rows=d['rows']
        self.assertEqual(len(rows),29)
        self.assertEqual(sum(r['role']=='canonical_calibration' for r in rows),25)
        self.assertEqual(sum(r['role']=='consumed_crystal_transfer' for r in rows),3)
        self.assertEqual(sum(r['role']=='separate_stress' for r in rows),1)
        self.assertFalse(d['crystals_or_stress_used_for_fit']);self.assertFalse(d['production_changed'])
        for v in ('operational','mathematical'):
            canonical=[r for r in rows if r['role']=='canonical_calibration']
            extrema={z:([r['pool'][v]['composite_R_model_kcal_mol'] for r in canonical if r['expected_class']==z]) for z in ('Ca','La')}
            self.assertEqual((len(extrema['Ca']),len(extrema['La'])),(14,11))
            self.assertEqual(d['variants'][v]['bands'],{'Ca_max':max(extrema['Ca']),'La_min':min(extrema['La'])})
            self.assertGreater(d['variants'][v]['gap_model_kcal_mol'],d['minimum_gap_model_kcal_mol'])
            for r in rows:
                self.assertEqual(r['status'],'available')
                self.assertEqual(r['own_reference_decisions'][v],r['expected_class']+'-supported')
        self.assertEqual(sum(r['pool_reused'] for r in rows),19)

    def test_final_actual_pool_algebra_and_complete_native_cells(self):
        d=read_json(BASE/'REFERENCE_v1.json')
        for r in d['rows']:
            self.assertEqual(pool.choose_rows(r['matrix'],adaptive.CANDIDATES),r['pool'])
        c=read_json(BASE/'pool_v1/collection_final.json')
        self.assertEqual((c['denominator'],c['available'],c['GFN2_complete'],c['MACE_calls_complete']),(10,10,80,20))
        self.assertEqual(c['native_candidate_reuses'],20)
        lm=read_json(BASE/'pool_v1/solvent/shard_0/manifest.json')
        self.assertEqual(lm['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
        self.assertEqual(len(lm['tasks']),80)

    def test_final_search_receipts_keep_boundary_and_actual_call_counts(self):
        rows=[read_json(p) for p in (BASE/'proposals_v1/proposals').glob('*/result.json')]
        self.assertEqual(len(rows),20)
        self.assertEqual({r['status'] for r in rows},{'proposal_available'})
        self.assertEqual(sum(r['boundary_flag'] for r in rows),11)
        self.assertTrue(all(not r['unconstrained_minimum_claimed'] for r in rows))
        actual=[x for r in rows for x in r['requests'] if not x['reused']]
        self.assertEqual(len(actual),389)
        self.assertEqual(len({x['result']['sha256'] for x in actual}),389)
        for x in actual:verify(x['result'])


if __name__=='__main__':unittest.main()
