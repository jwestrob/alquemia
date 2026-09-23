"""Actual prepared requests/archived energies only; no calculator is called."""
import copy
import json
from pathlib import Path
import statistics
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import pqq_three_source_execution as e
import adaptive_completion
import pqq_adaptive_candidate

BASE=ROOT/'workspaces/pqq_three_source_execution_20260923'
PLAN=BASE/'preflight_v1/plan.json'
UNSUPPORTED=BASE/'unsupported_fixture_v1/plan.json'


class ExplicitThreeExecution(unittest.TestCase):
    def test_exact_plm_preflight_has_no_canonical_or_label_gate(self):
        p,req,_,ref=e.checked(PLAN);result=e.dry_run(PLAN)
        self.assertEqual((result['sources'],result['groups'],result['prepared_origin_sources']),(6,2,6))
        self.assertEqual((result['origin_MACE_tasks'],result['origin_GFN2_tasks']),(12,24))
        self.assertEqual(result['declared_calls']['maximum_GFN2'],72)
        self.assertEqual(result['new_molecular_calls'],0)
        self.assertEqual(ref['reference_id'],e.preparation.REFERENCE_ID)
        self.assertTrue(all(c['source_conditioning_metal']=='La' and 'canonical_coordinate_match' not in c for c in req['cases']))
        self.assertTrue(all(v['expected_class'] is None for g in req['groups'] for v in g['source_evidence'].values()))

    def test_native_origin_state_input_recipe_and_no_archived_energy(self):
        p,_,_,_=e.checked(PLAN);engine=e.engine(PLAN)
        manifest=PLAN.parent/'origins/manifest.json';m=engine.validate_stage(manifest,PLAN)
        low=e.read_json(manifest.parent/'solvent/shard_0/manifest.json')
        self.assertEqual(low['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
        self.assertTrue(all(not t.get('native_reuse') for t in m['tasks']))
        self.assertEqual({(t['cell_id'],t['medium']) for t in low['tasks']},
                         {(t['task_id'],s) for t in m['tasks'] for s in ('vacuum','alpb')})
        for t in low['tasks']:
            self.assertIn('MaxIter 500',e.verify(t['input']).read_text())
            self.assertFalse(Path(t['output_path']).exists())
        self.assertEqual(len(p['case_ids']),6)

    def test_private_optimizer_keeps_existing_default_unchanged(self):
        old=copy.deepcopy(adaptive_completion.SETTINGS);protocol=pqq_adaptive_candidate.PROTOCOL
        engine=e.engine(PLAN)
        self.assertEqual(engine.completion.SETTINGS,e.preparation.precision.SETTINGS)
        self.assertEqual(engine.completion.SETTINGS['optimizer_ftol'],1e-8)
        self.assertEqual(adaptive_completion.SETTINGS,old)
        self.assertEqual(pqq_adaptive_candidate.PROTOCOL,protocol)

    def test_real_missing_member_creates_no_origin_calls(self):
        result=e.dry_run(UNSUPPORTED)
        self.assertEqual((result['sources'],result['prepared_origin_sources'],result['origin_GFN2_tasks']),(3,0,0))
        env=e.source_envelope(UNSUPPORTED)
        self.assertEqual(len(env['cases']),3)
        self.assertTrue(all(c['status']=='group_preparation_unavailable' for c in env['cases']))
        self.assertTrue(all('sample-4' in c['reason'] for c in env['cases']))

    def test_actual_incomplete_collection_keeps_all_sources_and_strict_groups(self):
        with tempfile.TemporaryDirectory(dir=BASE,prefix='actual_incomplete_report_') as d:
            output=Path(d)/'RESULT.json';e.collect(PLAN,output);data=e.read_json(output)
            self.assertEqual((data['source_denominator'],data['source_available'],data['group_denominator']),(6,0,2))
            self.assertTrue(all(r['union_origin_R_model_kcal_mol'] is None for r in data['rows']))
            self.assertTrue(all(r['source_evidence']['expected_class'] is None for r in data['rows']))
            self.assertTrue(all(g['variants']['operational']['R_model_kcal_mol'] is None for g in data['groups']))
            self.assertTrue(all(g['variants']['operational']['decision']=='unavailable' for g in data['groups']))
            self.assertFalse((PLAN.parent/'origins/collection.json').exists())

    def test_real_archived_pool_algebra_and_explicitly_corrupted_cell(self):
        # These are measured ten-member-union energies used only to test arithmetic,
        # not substituted into this proposed three-source scientific execution.
        archive=e.read_json(ROOT/'workspaces/pqq_union_execution_20260923/union-candidate_v4/pool/collection.json')
        ref=e.read_json(e.verify(e.read_json(PLAN)['reference']));rows=[]
        cases=[c for c in archive['cases'] if 'conditioned_La' in c['case_id']][:3]
        for case in cases:
            row=e.source_result({'case_id':case['case_id']},case,None,ref);rows.append(row)
            for mode in ('mathematical','operational'):
                r=case['pool'][mode]['composite_R_model_kcal_mol']
                self.assertEqual(row['variants'][mode]['R_model_kcal_mol'],r)
                self.assertEqual(row['variants'][mode]['delta_R_from_union_origin'],r-row['union_origin_R_model_kcal_mol'])
        group={'protein_id':'a0a3f2yly8-pqq-la_model','members':[c['case_id'] for c in cases],
               'physical_state_anchor':cases[0]['case_id'],'source_evidence':{}}
        result=e.group_results([group],rows,ref)[0]
        self.assertEqual(result['variants']['operational']['R_model_kcal_mol'],statistics.median(r['variants']['operational']['R_model_kcal_mol'] for r in rows))
        broken=copy.deepcopy(cases[0]);del broken['matrix']['Ca']['adaptive_La']
        with self.assertRaises(e.InvalidArtifact):e.source_result({'case_id':broken['case_id']},broken,None,ref)

    def test_corrupted_plan_rejects_scope_or_profile_change(self):
        original=e.read_json(PLAN)
        for key,value in [('profile','unqualified_profile'),('case_ids',original['case_ids'][:-1]),('cache_policy','archive_energy_ok')]:
            broken=copy.deepcopy(original);broken[key]=value
            with tempfile.TemporaryDirectory(dir=BASE,prefix='corrupted_actual_plan_') as d:
                path=Path(d)/'plan.json';path.write_text(json.dumps(broken))
                with self.assertRaises(e.InvalidArtifact):e.checked(path)

if __name__=='__main__':unittest.main()
