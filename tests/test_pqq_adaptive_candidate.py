"""Actual prepared/recovery artifacts; tests make no molecular calls."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import pqq_adaptive_candidate as a
REC=ROOT/'workspaces/adaptive_origin_recovery_20260923/run_v2'
SOURCE=ROOT/'diagnostics/pqq_fast_release_20260920/examples/1H4I_source.json'
PREP=ROOT/'workspaces/pqq_fast_release_20260920/standard_1H4I_v1/source_preparation/preparation.json'

class Candidate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(dir=ROOT/'workspaces/pqq_adaptive_candidate_20260923')
        cls.root=Path(cls.temp.name)
        req={'schema_version':a.REPLAY,'collection':a.record(REC/'pool/final_collection.json'),'case_ids':list(a.recovery.CASES)}
        a.put(cls.root/'replay.json',req)
        a.prepare(cls.root/'replay.json',cls.root/'replay')
        a.prepare(SOURCE,cls.root/'fresh')
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def test_exact_archive_replay_and_no_new_execution(self):
        plan=self.root/'replay/plan.json'
        self.assertEqual(a.execute(plan)['new_molecular_calls'],0)
        a.collect(plan,self.root/'replay_result.json')
        new=a.read_json(self.root/'replay_result.json');old=a.read_json(REC/'decisions.json')
        self.assertEqual(new['available'],2)
        for row,previous in zip(new['rows'],old['rows']):
            self.assertEqual(row['candidate']['operational']['R_model_kcal_mol'],previous['decisions']['operational']['R'])
            self.assertEqual(row['candidate']['operational']['decision'],'Ca-supported')
        a.report(self.root/'replay_result.json',self.root/'replay_report.md')

    def test_fresh_explicit_source_ready_without_energy(self):
        r=a.dry_run(self.root/'fresh/plan.json')
        self.assertEqual(r['mode'],'fresh_source_request');self.assertEqual(r['new_molecular_calls'],0)
        self.assertEqual(r['declared_maximum_calls']['GFN2'],12)
        self.assertEqual(r['fresh_preparation_and_selector_status'],'not_executed')
        a.collect(self.root/'fresh/plan.json',self.root/'not_executed.json')
        r=a.read_json(self.root/'not_executed.json')['rows'][0]
        self.assertEqual(r['status'],'unavailable');self.assertIsNone(r['candidate']['operational']['R_model_kcal_mol'])

    def test_actual_prepared_source_mapping_and_native_inputs(self):
        mp=a.origin_manifest(self.root/'fresh/plan.json',PREP,self.root/'prepared_origins')
        m=a.validate_stage(mp,self.root/'fresh/plan.json');self.assertEqual(len(m['tasks']),2);self.assertEqual(m['cases'][0]['status'],'prepared')
        low=a.read_json(mp.parent/'solvent/shard_0/manifest.json')
        self.assertEqual(len(low['tasks']),4);self.assertEqual(low['execution_resources'],{'mpi_ranks':8,'concurrent_tasks':4})
        for t in low['tasks']:
            body=a.verify(t['input']).read_text()
            self.assertEqual(body,a.pool.input_text(t['charge'],1,t['medium'],'native').replace('%scf\n','%scf\n MaxIter 500\n'))
            self.assertFalse(Path(t['output_path']).exists())
        a.paired(a.verify(m['tasks'][1]['xyz']),a.verify(m['tasks'][0]['xyz']),m['tasks'][1]['charge'],m['tasks'][0]['charge'])

    def test_actual_origin_force_selector_replays_both_cases(self):
        m=a.read_json(REC/'manifest.json')
        for case in m['cases']:
            pair=[t for t in m['tasks'] if t['case_id']==case['case_id']]
            selected=a.selected_tasks(pair,case['origins'],m)
            for new,old in zip(selected,pair):
                self.assertEqual(new['selector'],old['selector']);self.assertEqual(new['active_indices'],old['active_indices'])
                self.assertEqual(new['origin_reuse']['point']['gradient_kcal_mol_rad'],old['origin_reuse']['point']['gradient_kcal_mol_rad'])

    def test_corrupted_real_state_and_unsupported_normalization_rejected(self):
        m=a.read_json(REC/'manifest.json');case=m['cases'][0];pair=copy.deepcopy(m['tasks'][:2]);pair[0]['charge']+=2
        with self.assertRaises(a.InvalidArtifact):a.selected_tasks(pair,case['origins'],m)
        req=a.read_json(SOURCE);req['cases'][0]['normalization']='guess_metal_from_filename'
        with self.assertRaises(a.InvalidArtifact):a.check_request(req,a.standard.release(a.standard.DEFAULT_RELEASE),a.reference_data(a.DEFAULT_REFERENCE))

    def test_failed_actual_cell_cannot_become_origin_fallback(self):
        c=copy.deepcopy(a.read_json(REC/'pool/final_collection.json')['cases'][0])
        c['matrix']['La']['adaptive_Ca']['status']='unavailable'  # explicitly corrupted actual record
        selected=a.pool.choose_rows(c['matrix'],a.CANDIDATES)
        self.assertEqual(selected['status'],'unavailable');self.assertIsNone(selected['operational'])

    def test_actual_fresh_two_source_integration(self):
        run=ROOT/'workspaces/pqq_adaptive_candidate_20260923/two_crystals_v1'
        if not (run/'VERIFICATION.json').exists():self.skipTest('separately authorized fresh integration has not completed')
        v=a.read_json(run/'VERIFICATION.json');r=a.read_json(run/'result.json')
        self.assertEqual((v['available'],v['MACE_calls'],v['GFN2_attempts'],v['GFN2_complete']),(2,72,24,24))
        self.assertEqual(v['MACE_failed_requests'],0)
        self.assertLess(v['max_abs_origin_composite_difference_kcal_mol'],.01)
        self.assertLess(v['max_abs_pooled_R_difference_kcal_mol'],.1)
        self.assertEqual([x['candidate']['operational']['decision'] for x in r['rows']],['Ca-supported','La-supported'])
        for case in v['rows']:
            for ep in case['endpoints'].values():self.assertTrue(ep['source_atom_order_and_coordinates_exact'])
        for stage in ('origins','proposals','pool'):a.validate_stage(run/stage/'manifest.json',run/'plan.json')

    def test_corrupted_actual_stage_invalidates_method_cache(self):
        run=ROOT/'workspaces/pqq_adaptive_candidate_20260923/two_crystals_v1'
        m=a.read_json(run/'origins/manifest.json');m['model']['dtype']='float32'
        path=self.root/'corrupted_actual_method.json';a.put(path,m)
        with self.assertRaises(a.InvalidArtifact):a.validate_stage(path,run/'plan.json')

if __name__=='__main__':unittest.main()
