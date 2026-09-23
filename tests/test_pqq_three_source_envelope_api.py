"""Interface-only qualification against real prepared triples and v1 outputs."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import pqq_three_source_envelope_execution as p

BASE=ROOT/'workspaces/pqq_three_source_envelope_api_20260923/preflight_v1'
OLD=ROOT/'workspaces/pqq_three_source_envelope_execution_20260923/run_v1'
PREP=ROOT/'workspaces/pqq_three_source_envelope_20260923/plm_v1'


class ReusableEnvelopeAPI(unittest.TestCase):
    def test_actual_triples_preflight_separately_and_together(self):
        for name,n in [('07ab',3),('8344',3),('both',6)]:
            plan=BASE/name/'plan.json';r=p.dry_run(plan);data=p.read_json(plan)
            self.assertEqual(data['protocol_id'],p.PROTOCOL)
            self.assertEqual(r['sources'],n);self.assertEqual(r['groups'],n//3)
            self.assertEqual(r['origin_MACE_tasks'],2*n);self.assertEqual(r['origin_GFN2_tasks'],4*n)
            self.assertEqual(r['declared_calls']['maximum_GFN2'],12*n)
            self.assertEqual(r['submission_status'],'not_submitted')
            self.assertFalse(list(plan.parent.glob('execution_*.json')))
            self.assertFalse(list(plan.parent.rglob('endpoint.out')))
            self.assertFalse(list(plan.parent.rglob('mace_result.json')))

    def test_v2_preserves_every_old_coordinate_state_and_scalar_recipe(self):
        old=p.read_json(OLD/'origins/manifest.json');oldlow=p.read_json(OLD/'origins/solvent/shard_0/manifest.json')
        oldtasks={r['task_id']:r for r in old['tasks']};oldscalar={r['task_id']:r for r in oldlow['tasks']}
        for name in ('07ab','8344','both'):
            m=p.read_json(BASE/name/'origins/manifest.json');low=p.read_json(BASE/name/'origins/solvent/shard_0/manifest.json')
            for key in ('model','software','orca','reference','settings','proposal_settings','numerical_policy_id'):
                self.assertEqual(m[key],old[key])
            for t in m['tasks']:
                for key in ('xyz','charge','multiplicity','mapping','source_preparation'):
                    self.assertEqual(t[key],oldtasks[t['task_id']][key])
            for t in low['tasks']:
                prev=oldscalar[t['task_id']]
                self.assertEqual(p.verify(t['input']).read_bytes(),p.verify(prev['input']).read_bytes())
                self.assertEqual(p.verify(t['xyz']).read_bytes(),p.verify(prev['xyz']).read_bytes())

    def test_legacy_plan_and_result_remain_v1_and_replayable(self):
        plan,_,_,_=p.checked(OLD/'plan.json');self.assertEqual(plan['protocol_id'],p.LEGACY_PROTOCOL)
        self.assertEqual(p.dry_run(OLD/'plan.json')['sources'],6)
        data=p.read_json(OLD/'RESULT_1211626.json')
        self.assertEqual(data['protocol_id'],p.LEGACY_PROTOCOL)
        self.assertEqual(data['source_available'],6)
        engine=p.engine(OLD/'plan.json')
        for row in data['rows']:
            self.assertEqual(engine.pool.choose_rows(row['matrix'],['origin','adaptive_Ca','adaptive_La']),row['pool'])

    def test_actual_metadata_alias_has_no_canonical_or_plm_identifier_gate(self):
        # An explicitly identified metadata-only alias of the same real protein,
        # not a new structure/label or fabricated molecular result.
        original=p.read_json(PREP/(p.PROTEINS[0]+'_REQUEST.json'))
        value=copy.deepcopy(original);value['protein_id']='soil.site_07-alias'
        self.assertEqual(p.preparation.check_request(value)['protein_id'],'soil.site_07-alias')
        with tempfile.TemporaryDirectory() as td:
            td=Path(td);rp=td/'metadata_alias_of_actual_request.json';p.candidate.put(rp,value)
            pp=p.read_json(PREP/p.PROTEINS[0]/'PREPARATION.json');pp['request']=p.record(rp)
            pp['metadata_alias_only_of']=p.record(PREP/p.PROTEINS[0]/'PREPARATION.json')
            prep=td/'metadata_alias_of_actual_preparation.json';p.candidate.put(prep,pp)
            declared=p.interface().request_from_preparations([prep]);request=td/'execution_request.json';p.candidate.put(request,declared)
            plan=p.read_json(BASE/'07ab/plan.json');plan['request']=p.record(request)
            path=td/'metadata_alias_plan.json';p.candidate.put(path,plan)
            _,accepted,_,_=p.checked(path)
            self.assertEqual(accepted['groups'][0]['protein_id'],'soil.site_07-alias')
            self.assertEqual(accepted['cases'],original['cases'])
        for bad in ('../bad','bad/name','two words'):
            value['protein_id']=bad
            with self.assertRaises(p.InvalidArtifact):p.preparation.check_request(value)

    def test_real_other_sequence_and_changed_state_are_rejected(self):
        left=p.read_json(PREP/(p.PROTEINS[0]+'_REQUEST.json'));right=p.read_json(PREP/(p.PROTEINS[1]+'_REQUEST.json'))
        self.assertNotEqual(p.preparation.original.source_identity(left['cases'][0],left['config']),
                            p.preparation.original.source_identity(right['cases'][0],right['config']))
        bad=copy.deepcopy(left);bad['cases'][-1]=copy.deepcopy(right['cases'][0]);bad['cases'].sort(key=lambda x:x['case_id'])
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'explicitly_corrupted_mixed_real_proteins.json';p.candidate.put(path,{'sources':bad['cases']})
            bad['source_descriptors']=p.record(path)
            with self.assertRaisesRegex(p.InvalidArtifact,'source protein sequence/numbering/site roles/assembly differ'):
                p.preparation.check_request(bad)
            prep=p.read_json(PREP/p.PROTEINS[0]/'PREPARATION.json');prep['cases'][1]['state_key']='explicitly_corrupted'
            path=Path(td)/'explicitly_corrupted_real_state.json';p.candidate.put(path,prep)
            with self.assertRaises(p.InvalidArtifact):p.preparation.dry_run(path)

    def test_duplicate_groups_and_partial_medians_do_not_fall_back(self):
        path=PREP/p.PROTEINS[0]/'PREPARATION.json'
        with self.assertRaises(p.InvalidArtifact):p.interface().request_from_preparations([path,path])
        data=p.read_json(OLD/'RESULT_1211626.json');rows=copy.deepcopy(data['rows']);req=p.read_json(OLD/'REQUEST.json')
        rows[0]['status']='unavailable'
        for variant in ('mathematical','operational'):rows[0]['variants'][variant]['R_model_kcal_mol']=None
        groups=p.interface().group_results(req['groups'],rows,p.read_json(p.verify(data['reference'])))
        self.assertEqual(groups[0]['available'],2)
        self.assertIsNone(groups[0]['variants']['operational']['R_model_kcal_mol'])
        self.assertEqual(groups[0]['variants']['operational']['decision'],'unavailable')
        self.assertEqual(groups[1]['available'],3)


if __name__=='__main__':unittest.main()
