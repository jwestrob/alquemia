"""No new folds, protonation or quantum calculations: immutable fixture tests."""
import copy
import os
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import prepare_batch as prep
from geometry_checks import check,unchanged_non_H

OLD=prep.A/'workspaces/plm_xoxf_fixed_core_20260916'
EQ=Path('/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue')

class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.w=prep.wrapper();cls.review=prep.read(OLD/'selection_review.json')
        cls.inventory=prep.by_id(prep.read(EQ/'xoxf_all/inventory/manifest.json')['targets'])
        cls.pp=prep.read(OLD/'prepared/prepared_pairs.json')
        cls.carves=[prep.read(c['carve_manifest']['path']) for c in cls.pp['cases']]

    def test_raw_six_models_match_existing_metrics(self):
        for previous in self.review['models']:
            entry={'sample':previous['sample'],'seed':101,'source_cif':previous['source_cif'],'summary_confidences':previous['summary_confidence']}
            current=prep.raw_model(self.w,self.inventory[previous['target_id']],entry)
            for key in ('CN','protein_CN','protein_La_iptm','direct_N','source_cif'):
                self.assertEqual(current[key],previous[key])

    def test_map_actual_roles_and_partners(self):
        expected={'PQQSEQ_242fa05e3ffc20087d42':[190,274,316,318,343],
                  'PQQSEQ_faa97386262eec4316fc':[186,252,302,304,329]}
        for selected in self.review['selected']:
            uid=selected['target_id'];roles,partners=prep.map_roles(self.w,self.inventory[uid],selected['selection_metrics'])
            self.assertEqual([roles[r]['resnum'] for r in self.w.fixed.ROLE_ORDER],expected[uid])
            self.assertEqual(roles['catalytic_asp_cationic_partner']['resname'],'ARG')
            self.assertLessEqual(partners[0]['distance_A'],3.5)

    def test_bad_identity_or_unverified_role_rejected(self):
        selected=self.review['selected'][0];target=copy.deepcopy(self.inventory[selected['target_id']])
        target['role_mappings'][0]['amino_acid']='D'
        with self.assertRaisesRegex(ValueError,'differs'):prep.map_roles(self.w,target,selected['selection_metrics'])
        target=copy.deepcopy(self.inventory[selected['target_id']]);target['role_mappings'][0]['sequence_identity_verified']=False
        with self.assertRaisesRegex(ValueError,'unverified'):prep.map_roles(self.w,target,selected['selection_metrics'])

    def test_choose_before_admission_no_replacement(self):
        rows=[{'CN':7,'protein_La_iptm':.99,'protein_CN':4,'direct_N':1,'sample':0},
              {'CN':8,'protein_La_iptm':.8,'protein_CN':5,'direct_N':1,'sample':1}]
        self.assertEqual(prep.select(rows)['sample'],1)
        self.assertFalse(prep.admission(prep.select(rows)))

    def test_sample_tie_break(self):
        rows=[{'CN':7,'protein_La_iptm':.99,'protein_CN':4,'sample':i} for i in (2,0,1)]
        self.assertEqual(prep.select(rows)['sample'],0)

    def test_original_malformed_H_rejected(self):
        path=prep.A/'workspaces/plm_adh9_af3_20260916/preparation/pairs/PQQSEQ_13d74836d4b7a3e02140_AF3_sample1/PQQSEQ_13d74836d4b7a3e02140_AF3_sample1_carve_manifest.json'
        self.assertEqual(check(prep.read(path)['qm_fragments'])['status'],'FAIL')

    def test_actual_clean_geometry_passes(self):
        for cm in self.carves:self.assertEqual(check(cm['qm_fragments'])['status'],'PASS')

    def test_frozen_atom_change_rejected(self):
        old=self.carves[0]['qm_fragments'];new=copy.deepcopy(old);new[0]['atom_records'][0]['xyz_A'][0]+=.001
        with self.assertRaisesRegex(ValueError,'Frozen'):unchanged_non_H(old,new)

    def test_hash_mismatch_and_duplicate_ids_rejected(self):
        rec=prep.record(OLD/'selection_review.json');rec['sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'changed'):prep.verify(rec)
        with self.assertRaisesRegex(ValueError,'Duplicate'):prep.by_id([{'target_id':'x'},{'target_id':'x'}])

    def test_worker_limit_single_thread_login(self):
        with patch.dict(os.environ,{'SLURM_CPUS_ON_NODE':'1'}):
            self.assertEqual(prep.worker_count('auto',174)[0],1)

    def test_missing_model_outcome_retained_without_preparation(self):
        uid=next(iter(prep.REUSE_IDS));target=self.inventory[uid]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'adapter_inputs').mkdir()
            candidate,outcome=prep.choose_target(self.w,target,None,root,{})
            self.assertIsNone(candidate);self.assertEqual(outcome['status'],'preparation_failed')
            self.assertEqual(outcome['failure_stage'],'folding_output')
            self.assertTrue(Path(outcome['selection_review']['path']).is_file())

    def test_zero_ready_batch_accounts_all176_without_protonation(self):
        inventory_path=EQ/'xoxf_all/inventory/manifest.json';approval_path=EQ/'xoxf_all/authorization.json'
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);fold_path=root/'fold_input.json';results_path=root/'fold_results.json';out=root/'output'
            fold_targets=[dict(t,action='reuse_complete' if uid in prep.REUSE_IDS else 'fold') for uid,t in self.inventory.items()]
            prep.write(fold_path,{'targets':fold_targets})
            prep.write(results_path,{'schema_version':'plm.xoxf_all.fold_results.v1','inventory':prep.record(inventory_path),
                'manifest':prep.record(fold_path),'targets':[]})
            args=['--inventory',str(inventory_path),'--inventory-sha256',prep.record(inventory_path)['sha256'],
                  '--fold-input-manifest',str(fold_path),'--fold-input-sha256',prep.record(fold_path)['sha256'],
                  '--fold-results',str(results_path),'--approval',str(approval_path),
                  '--approval-sha256',prep.record(approval_path)['sha256'],'--output',str(out)]
            with patch.dict(os.environ,{'SLURM_JOB_ID':'synthetic_no_compute_fixture'}):prep.main(args)
            result=prep.read(out/'prepared_pairs.json');outcomes=prep.read(out/'target_outcomes.json')
            self.assertEqual(result['cases'],[]);self.assertEqual(result['reused_count'],2)
            self.assertEqual(result['unsupported_count'],0);self.assertEqual(result['failed_count'],174)
            self.assertEqual(len(outcomes['per_target']),176)
            self.assertFalse(any((out/'original_preparation').iterdir()))

    def test_actual_upstream_missing_keeps_reason_and_is_not_chemistry(self):
        folds=prep.by_id(prep.read(EQ/'xoxf_all/folding/fold_results.json')['targets'])
        uid=next(uid for uid in folds if uid not in prep.REUSE_IDS)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'adapter_inputs').mkdir()
            candidate,outcome=prep.choose_target(self.w,self.inventory[uid],folds[uid],root,{})
            self.assertIsNone(candidate);self.assertEqual(outcome['status'],'preparation_failed')
            self.assertEqual(outcome['upstream_fold_state'],'missing')
            self.assertIn('GPU stage produced no target outcome',outcome['reason'])

    def test_selected_model_gate_failure_remains_unsupported(self):
        old=self.review['models'][0];target=self.inventory[old['target_id']]
        models=[{'sample':i,'seed':101,'status':'complete'} for i in range(3)]
        fold={'sequence_sha256':target['sequence_sha256'],'state':'complete','models':models}
        def low_cn(w,t,e):return dict(old,sample=e['sample'],CN=6)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'adapter_inputs').mkdir()
            with patch('prepare_batch.raw_model',side_effect=low_cn):
                candidate,outcome=prep.choose_target(self.w,target,fold,root,{})
            self.assertIsNone(candidate);self.assertEqual(outcome['status'],'unsupported')
            self.assertEqual(outcome['failure_stage'],'selected_model_admission')
            self.assertIn('no substitution',outcome['reason'])

    def test_native_worker_crash_retained_without_retry(self):
        target=prep.read(OLD/'candidate_manifest.json')['targets'][0]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'worker_dispatch').mkdir();(root/'worker_results').mkdir()
            with patch('prepare_batch.subprocess.run',return_value=subprocess.CompletedProcess([], -11)) as run:
                result=prep.dispatch_worker((target,'unused_fixture_pins',str(root),'unused_fixture_approval'))
            self.assertEqual(run.call_count,1);self.assertIsNone(result['case'])
            self.assertEqual(result['outcome']['status'],'preparation_failed')
            self.assertIn('-11',result['outcome']['reason'])
            receipt=prep.read(result['outcome']['dispatch_receipt']['path'])
            self.assertEqual(receipt['retry_count'],0);self.assertEqual(receipt['returncode'],-11)

    def test_recorded_worker_failure_retained_without_retry(self):
        target=prep.read(OLD/'candidate_manifest.json')['targets'][0]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'worker_dispatch').mkdir();(root/'worker_results').mkdir()
            expected={'case':None,'outcome':{'target_id':target['target_id'],'status':'preparation_failed','reason':'synthetic retained preparation error'},'geometry':{'status':'UNSUPPORTED'}}
            prep.write(root/'worker_results'/f'{target["target_id"]}.json',expected)
            with patch('prepare_batch.subprocess.run',return_value=subprocess.CompletedProcess([],0)) as run:
                result=prep.dispatch_worker((target,'unused_fixture_pins',str(root),'unused_fixture_approval'))
            self.assertEqual(run.call_count,1);self.assertEqual(result['outcome']['reason'],expected['outcome']['reason'])

if __name__=='__main__':unittest.main(verbosity=2)
