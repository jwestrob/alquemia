"""Actual preparation and archived molecular receipts; no fabricated scientific output."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import statistics
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import pqq_three_source_envelope_execution as p

PLAN=ROOT/'workspaces/pqq_three_source_envelope_execution_20260923/run_v1/plan.json'


class EnvelopeExecution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan,cls.request,_,cls.reference=p.checked(PLAN)
        archive=p.read_json(p.verify(cls.reference['collection']))
        cls.real=next(c for c in archive['cases'] if c['source']['origin_row']['source_case_id']=='1H4I')

    def test_actual_six_source_preflight_scope_state_recipe_and_unknown_labels(self):
        r=p.dry_run(PLAN);self.assertEqual(r['sources'],6);self.assertEqual(r['groups'],2)
        self.assertEqual((r['origin_MACE_tasks'],r['origin_GFN2_tasks']),(12,24))
        self.assertEqual(r['declared_calls']['maximum_GFN2'],72)
        for g in self.request['groups']:
            self.assertTrue(all(v['expected_class'] is None for v in g['source_evidence'].values()))
        low=p.read_json(PLAN.parent/'origins/solvent/shard_0/manifest.json')
        self.assertEqual(low['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
        for t in low['tasks']:
            self.assertEqual(p.verify(t['input']).read_text(),p.scalar.recipe(t['charge'],t['medium'],'fresh'))

    def test_unchanged_selector_replays_actual_envelope_forces(self):
        # Existing source-to-score warm receipts use the format emitted by this
        # executor. The earlier envelope study's legacy q0 format is different.
        root=ROOT/'workspaces/pqq_union_execution_20260923/union-candidate_v4'
        archived=p.read_json(root/'proposals/manifest.json');c=p.read_json(root/'origins/collection.json')['cases'][0]
        cid=c['case_id'];origins=p.read_json(root/'origins/manifest.json')
        pair=[next(t for t in origins['tasks'] if (t['case_id'],t['metal'])==(cid,z)) for z in ('Ca','La')]
        engine=p.engine(PLAN)
        tasks=engine.selected_tasks(pair,{z:c['matrix'][z]['origin'] for z in ('Ca','La')},{'model':archived['model']})
        expected=next(t['active_mode_ids'] for t in archived['tasks'] if t['case_id']==cid)
        self.assertEqual([t['active_mode_ids'] for t in tasks],[expected,expected])
        self.assertEqual(engine.completion.SETTINGS,p.preparation.precision.SETTINGS)

    def test_actual_strict_scalar_receipt_passes_active_tolerance_rank_and_state(self):
        low=self.real['matrix']['Ca']['origin']['low']['alpb']
        result=p.engine(PLAN).pool.completed(p.verify(low['manifest']),low['task_id'])
        self.assertEqual(result['status'],'complete');self.assertEqual(result['observed_TolE_hartree'],1e-10)
        self.assertEqual(result['energy_hartree'],low['energy_hartree'])

    def test_unavailable_collection_keeps_all_sources_and_groups_without_fallback(self):
        ex=p.interface();ref=self.reference
        rows=[ex.source_result(s,None,None,ref) for s in self.request['cases']]
        groups=ex.group_results(self.request['groups'],rows,ref)
        self.assertEqual(len(rows),6);self.assertEqual(len(groups),2)
        self.assertTrue(all(r['status']=='unavailable' and r['variants']['operational']['R_model_kcal_mol'] is None for r in rows))
        self.assertTrue(all(g['available']==0 and g['variants']['operational']['R_model_kcal_mol'] is None for g in groups))

    def test_corrupted_real_plan_policy_and_scope_rejected(self):
        for name in ('profile','count','reference'):
            value=copy.deepcopy(self.plan)
            if name=='profile':value['profile']='explicitly_corrupted_loose_native'
            elif name=='count':value['declared_calls']['maximum_GFN2']=73
            else:value['reference']=value['rank_qualification']
            with tempfile.TemporaryDirectory() as td:
                path=Path(td)/'explicitly_corrupted_real_plan.json';p.candidate.put(path,value)
                with self.assertRaises(p.InvalidArtifact):p.checked(path)


@unittest.skipUnless((PLAN.parent/'RESULT_1211626.json').exists(),'actual six-source molecular collection not yet available')
class CompletedIntegration(unittest.TestCase):
    def test_actual_complete_matrices_direct_score_sign_and_unknown_labels(self):
        from mace_hybrid import EV_TO_KCAL
        from affordable_common import HA_TO_KCAL
        result=p.read_json(PLAN.parent/'RESULT_1211626.json')
        self.assertEqual((result['source_denominator'],result['source_available'],result['group_denominator']),(6,6,2))
        self.assertFalse(result['biological_accuracy_evaluated']);self.assertFalse(result['original_DFT_and_production_results_changed'])
        self.assertEqual(result['candidate_status'],'experimental_not_promoted')
        for row in result['rows']:
            self.assertIsNone(row['source_evidence']['expected_class'])
            for variant in ('mathematical','operational'):
                components=[]
                for metal in ('Ca','La'):
                    selected=row['pool']['rows'][metal][variant+'_candidate']
                    components.append(row['matrix'][metal][selected]['components'])
                ca,la=components
                direct=(ca['MACE_eV']-la['MACE_eV'])*EV_TO_KCAL+((ca['GFN2_ALPB_hartree']-ca['GFN2_vacuum_hartree'])-(la['GFN2_ALPB_hartree']-la['GFN2_vacuum_hartree']))*HA_TO_KCAL
                self.assertAlmostEqual(direct,row['variants'][variant]['R_model_kcal_mol'],places=7)
            self.assertEqual({q['id'] for q in row['candidates']},{'origin','adaptive_Ca','adaptive_La'})
            for metal in ('Ca','La'):
                for cell in row['matrix'][metal].values():
                    self.assertEqual(cell['status'],'complete')
                    for low in cell['low'].values():self.assertEqual(low['observed_TolE_hartree'],1e-10)

    def test_actual_strict_three_source_summaries_and_physical_admission(self):
        result=p.read_json(PLAN.parent/'RESULT_1211626.json');rows={r['case_id']:r for r in result['rows']}
        for group in result['groups']:
            self.assertEqual(group['available'],3);self.assertFalse(group['biological_accuracy_evaluated'])
            original=[rows[c]['union_origin_R_model_kcal_mol'] for c in group['members']]
            self.assertEqual(group['union_origin']['median_R_model_kcal_mol'],statistics.median(original))
            for variant in ('mathematical','operational'):
                values=[rows[c]['variants'][variant]['R_model_kcal_mol'] for c in group['members']]
                self.assertEqual(group['variants'][variant]['R_model_kcal_mol'],statistics.median(values))
                self.assertEqual(group['variants'][variant]['range_kcal_mol'],max(values)-min(values))
        paths=list((PLAN.parent/'proposals/proposals').glob('*/result.json'));self.assertEqual(len(paths),12)
        for path in paths:
            d=p.read_json(path);self.assertEqual(d['status'],'proposal_available')
            self.assertTrue(d['optimizer']['success']);self.assertTrue(d['final_geometry']['physical_feasible'])
            self.assertLessEqual(d['final_geometry']['maximum_heavy_displacement_A'],.8000001)


if __name__=='__main__':unittest.main()
