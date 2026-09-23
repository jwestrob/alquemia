"""Finite fresh-execution preflight using real prepared A0A3 sources only."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import pqq_union_execution as e
from mace_site_kinematics import Kinematics
import adaptive_completion

WORK=ROOT/'workspaces/pqq_union_execution_20260923'


class UnionExecution(unittest.TestCase):
    def test_fixed_fresh_counts_and_profiles(self):
        calls=[]
        for arm in e.ARMS:
            p=WORK/(arm+'_v4')/'plan.json';v=e.dry_run(p);calls.append(v['declared_calls'])
            self.assertEqual(v['sources'],10);self.assertEqual(v['new_molecular_calls'],0)
            self.assertEqual(v['fresh_source_preparation_status'],'prepared' if (p.parent/'source_preparation/preparation.json').exists() else 'not_executed')
        self.assertEqual(sum(x['maximum_GFN2'] for x in calls),160)
        self.assertEqual(sum(x['native_MACE_scalar']+x['native_MACE_q0_gradient']+x['maximum_cross_MACE'] for x in calls),60)
        self.assertEqual(sum(x['bounded_native_MACE_searches'] for x in calls),20)

    def test_actual_archived_mapping_preflight(self):
        plan=WORK/'union-candidate_v4/plan.json';engine=e.candidate_engine(plan)
        mp=WORK/'mapping_preflight_v4/origins/manifest.json';m=engine.validate_stage(mp,plan)
        self.assertEqual(len(m['tasks']),20)
        for t in m['tasks']:
            kin=Kinematics(e.read_json(e.verify(t['mapping']))['context'])
            coords=kin.evaluate(np.zeros(len(kin.modes)))[1]
            self.assertTrue(np.allclose(coords,[a[1:] for a in e.xyz(e.verify(t['xyz']))],atol=1e-12,rtol=0))
        low=e.read_json(mp.parent/'solvent/shard_0/manifest.json')
        self.assertEqual(len(low['tasks']),40)
        self.assertEqual(low['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
        self.assertFalse(any(Path(t['output_path']).exists() for t in low['tasks']))

    def test_private_profile_preserves_original_optimizer(self):
        old=copy.deepcopy(adaptive_completion.SETTINGS)
        engine=e.candidate_engine(WORK/'union-candidate_v4/plan.json')
        self.assertEqual(engine.completion.SETTINGS['optimizer_ftol'],1e-8)
        self.assertEqual(adaptive_completion.SETTINGS,old)
        self.assertNotEqual(engine.completion.SETTINGS['optimizer_ftol'],old['optimizer_ftol'])

    def test_corrupted_actual_policy_is_rejected(self):
        base=e.read_json(WORK/'union-candidate_v4/plan.json')
        for kind in ('precision','calls','reuse'):
            m=copy.deepcopy(base)
            if kind=='precision':m['settings']['optimizer_ftol']=1e-6
            elif kind=='calls':m['declared_calls']['maximum_GFN2']=121
            else:m['cache_policy']='use archived scores'
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/'explicitly_corrupted_actual_plan.json';e.write_new(p,m)
                with self.assertRaises(e.InvalidArtifact):e.checked(p)

    def test_real_unexecuted_paths_report_missing_without_fallback(self):
        for arm in e.ARMS:
            plan=WORK/(arm+'_v3')/'plan.json'
            with tempfile.TemporaryDirectory() as d:
                out=Path(d)/'missing_collection.json';v=e.collect(plan,out);r=e.read_json(out)
                self.assertEqual(v['available'],0);self.assertEqual(r['denominator'],10)
                self.assertTrue(all(x['R'] is None for x in r['rows']))
                self.assertFalse(r['archive_energy_reuse'])

    def test_unexecuted_matched_comparison_keeps_strict_members(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);paths=[]
            for arm in e.ARMS:
                out=root/(arm+'.json');e.collect(WORK/(arm+'_v3')/'plan.json',out);paths.append(out)
            result=root/'comparison.json';e.compare(paths[0],paths[1],result);r=e.read_json(result)
            self.assertEqual(r['source_denominator'],10);self.assertEqual(r['common_coverage'],0)
            self.assertEqual(r['independent_protein_groups'],1)
            self.assertEqual(r['all_counts']['union_candidate'],{'unavailable':10})
            self.assertEqual(len(r['aggregates']),7)
            self.assertTrue(all(x['methods']['union_candidate']['R'] is None for x in r['aggregates']))

    def test_archived_preparation_cannot_claim_same_run_restart(self):
        plan=WORK/'union-candidate_v4/plan.json'
        source=WORK/'mapping_preflight_v4/ARCHIVED_PREPARATION_ONLY.json'
        with self.assertRaises(e.InvalidArtifact):e.same_run_preparation(plan,source)

    def test_released_scanner_restart_is_collect_only_after_actual_start(self):
        release=e.standard.release(e.standard.DEFAULT_RELEASE)
        actual=e.read_json(e.verify(release['artifacts']['release_result']))
        collection=e.read_json(e.verify(actual['score_collection']))
        manifest=e.verify(collection['manifest'])
        state=e.static_execution_state(manifest)
        self.assertEqual(state['status'],'already_started_collect_only')
        self.assertFalse(state['automatic_molecular_retry'])
        self.assertTrue(any(x.endswith('execution_complete.json') for x in state['evidence_paths']))
        not_started=WORK/'mapping_preflight_v4/static/scoring/manifest.json'
        self.assertEqual(e.static_execution_state(not_started)['status'],'not_started')

    def test_actual_static_component_files_are_not_top_level_results(self):
        plan=WORK/'released-static_v4/plan.json'
        if not (plan.parent/'RESULT.json').exists():self.skipTest('actual static output unavailable')
        files=list((plan.parent/'prepared_score/scoring').glob('result_*.json'))
        self.assertGreater(len(files),1)
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'actual_recollection.json';value=e.collect(plan,out)
            self.assertEqual(value['available'],10)
            self.assertEqual(e.read_json(out)['rows'],e.read_json(plan.parent/'RESULT.json')['rows'])

    def test_actual_candidate_algebra_and_all_cells(self):
        path=WORK/'union-candidate_v4/RESULT.json'
        if not path.exists():self.skipTest('actual candidate output unavailable')
        data=e.read_json(path)
        for row in data['rows']:
            self.assertEqual(row['status'],'available')
            ca,la=[row['pool']['rows'][z]['operational_candidate'] for z in ('Ca','La')]
            expected=e.shared_pool.score(row['matrix']['Ca'][ca]['components'],row['matrix']['La'][la]['components'])
            self.assertEqual(expected['composite_R_model_kcal_mol'],row['R'])
            self.assertEqual(sum(cell['status']=='complete' for v in row['matrix'].values() for cell in v.values()),6)
        costs=e.read_json(WORK/'COSTS_v1.json')
        self.assertEqual(costs['new_GFN2_attempts'],160)
        self.assertEqual(costs['new_DFT'],0)
        self.assertEqual(sum(a['GFN2_normal'] for a in costs['arms']),160)

    def test_scientific_integration_requires_actual_execution(self):
        paths=[WORK/(arm+'_v4')/'RESULT.json' for arm in e.ARMS]
        if not all(p.exists() for p in paths):self.skipTest('fresh matched molecular jobs have not produced both final collections')
        for p in paths:
            r=e.read_json(p);self.assertEqual(r['denominator'],10);self.assertEqual(r['available'],10)
            self.assertFalse(r['archive_energy_reuse']);self.assertTrue(r['execution_receipts'])
            e.same_run_preparation(e.verify(r['plan']),p.parent/'source_preparation/preparation.json')
        comp=e.read_json(WORK/'COMPARISON_v1.json')
        self.assertEqual(comp['common_coverage'],10)
        self.assertFalse(comp['calibration_refitted'])


if __name__=='__main__':unittest.main()
