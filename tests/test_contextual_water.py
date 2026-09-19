"""Reusable preparation tested on pinned real alpha/GGR/PQQ artifacts."""
from __future__ import annotations
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from contextual_water_context import discover,policy,verify_proposal_geometry
from contextual_water_prepare import validate as validate_prepared,collect as collect_prepared
from contextual_water_score import validate as validate_score,collect as collect_score


class ContextualWaterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=ROOT/'workspaces/contextual_water_20260919/replay_v3'
        cls.pm=cls.root/'manifest.json';cls.pc=read_json(cls.root/'collection.json')
        cls.sm=ROOT/'workspaces/contextual_water_20260919/score_replay_v1/manifest.json'
        cls.sc=read_json(cls.sm.parent/'collection.json')
        cls.req=read_json(ROOT/'diagnostics/contextual_water_20260919/REPLAY_REQUEST.json')

    def test_real_replay_uses_no_new_evaluations(self):
        p=validate_prepared(self.pm);s=validate_score(self.sm)
        self.assertEqual((p['tasks'],p['reused_proposals'],p['zero_water_identities']),(0,4,30))
        self.assertEqual((s['tasks'],s['reused_endpoints'],s['zero_water_identities']),(0,8,30))
        self.assertEqual(self.pc,collect_prepared(self.pm));self.assertEqual(self.sc,collect_score(self.sm))

    def test_dry_pqq_and_ggr_are_exact_source_identity(self):
        dry=[c for c in self.pc['cases'] if c['operation']=='exact_zero_water_identity']
        self.assertEqual(len(dry),30)
        for c in dry:
            original=read_json(verify(c['source_preparation']))
            self.assertFalse(original['explicit_waters']);self.assertEqual(c['endpoints'],original['endpoints'])
            self.assertTrue(c['common_paired_geometry'])
            self.assertFalse(any(a['kind']=='retained_site_water' for a in original['physical_atoms']))

    def test_water_transfer_matches_actual_scanner_xyz(self):
        old=read_json(ROOT/'workspaces/hydration_scanner_20260919/mace_v2/manifest.json')
        for c in self.pc['cases']:
            if c['operation']=='exact_zero_water_identity':continue
            prep=read_json(verify(c['source_preparation']))
            for metal,e in c['endpoints'].items():
                oldtask=next(t for t in old['tasks'] if t['case_id']==c['case_id'] and t['metal']==metal and t['preparation_variant']=='contextual' and t['position']=='bound')
                self.assertEqual(xyz(verify(e['xyz'])),xyz(verify(oldtask['source_xyz'])))
                initial=xyz(verify(prep['endpoints'][metal]['xyz']));actual=xyz(verify(e['xyz']))
                mobile={a['whole_index'] for a in e['water_transfer']['whole_H_mapping']}
                self.assertEqual(len(mobile),2*c['explicit_water_count'])
                self.assertTrue(all(actual[i]==initial[i] for i in range(len(actual)) if i not in mobile))
                self.assertEqual(e['charge'],prep['endpoints'][metal]['charge'])

    def test_actual_four_state_descriptor_reproduces(self):
        old=read_json(ROOT/'workspaces/hydration_scanner_20260919/mace_v2/collection_job_1202084.json')
        for c in self.sc['cases']:
            if c['operation']=='exact_zero_water_identity':
                self.assertIsNone(c['interaction_R_model_kcal']);continue
            historical=next(r for r in old['cases'] if c['case_id']=='ALPHA_'+r['case'])['MACE_contextual']
            for key in ('interaction_R_model_kcal','bound_total_R_model_kcal','detached_environment_difference_model_kcal'):
                self.assertAlmostEqual(c[key],historical[key],places=10)
            self.assertFalse(c['common_paired_geometry'])
            self.assertNotEqual(c['detached_environment_difference_model_kcal'],0)
            self.assertLess(abs(c['component_closure_error_model_kcal']),1e-10)

    def test_generic_identifiers_rebuild_same_real_contexts(self):
        wet=[copy.deepcopy(c) for c in self.req['cases'] if c.get('reuse_proposals')]
        renamed={}
        for i,c in enumerate(wet):renamed[c['case_id']]='arbitrary_sample_'+str(i);c['case_id']=renamed[c['case_id']]
        generated=discover(wet,policy(verify(self.req['policy_source'])))
        m=read_json(self.pm)
        self.assertEqual(set(generated),set(renamed.values()))
        for c in m['cases']:
            if c['case_id'] not in renamed:continue
            old=read_json(verify(c['context']));new=generated[renamed[c['case_id']]]
            np.testing.assert_allclose([a[1:] for a in old['atoms']],[a[1:] for a in new['atoms']],rtol=0,atol=0)
            self.assertEqual(old['charges'],new['charges'])

    def test_unconverged_actual_proposal_rejected(self):
        m=read_json(self.pm);t=m['proposal_tasks'][0];r=read_json(verify(m['reused'][t['task_id']]))
        damaged=copy.deepcopy(r);damaged['optimization_eligible']=False
        with self.assertRaises(InvalidArtifact):verify_proposal_geometry(t,damaged)

    def test_corrupted_real_oxygen_is_rejected(self):
        from mace_hybrid import write_xyz
        m=read_json(self.pm);t=m['proposal_tasks'][0];r=read_json(verify(m['reused'][t['task_id']]))
        rows=xyz(verify(r['proposed_xyz']));i=t['water_groups'][0]['oxygen_index'];a=rows[i]
        corrupted=list(rows);corrupted[i]=(a[0],a[1]+.01,a[2],a[3])
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'explicitly_corrupted_real_oxygen.xyz';write_xyz(path,corrupted)
            damaged=copy.deepcopy(r);damaged['proposed_xyz']=record(path)
            with self.assertRaises(InvalidArtifact):verify_proposal_geometry(t,damaged)

    def test_missing_actual_detached_endpoint_does_not_fallback(self):
        # Explicitly withheld real receipt: recovery behavior, not new evidence.
        m=read_json(self.sm);t=next(t for t in m['scorer_tasks'] if t['position']=='detached')
        damaged=copy.deepcopy(m);damaged['reused'].pop(t['task_id']);damaged['tasks']=[t]
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'withheld_real_receipt_manifest.json';write_new(path,damaged)
            result=collect_score(path)
            row=next(r for r in result['cases'] if r['case_id']==t['case_id'])
            self.assertEqual(result['status'],'incomplete');self.assertEqual(row['status'],'unavailable')
            self.assertIsNone(row['interaction_R_model_kcal'])

    def test_changed_method_invalidates_real_manifest(self):
        m=read_json(self.pm);m['optimization']['gtol_eV_radian']=.002
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'explicitly_corrupted_method.json';write_new(path,m)
            with self.assertRaises(InvalidArtifact):validate_prepared(path)


if __name__=='__main__':unittest.main()
