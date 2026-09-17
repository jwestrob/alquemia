"""Prepared-input interface tests: real states, real receipts, corrupted copies."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_omol_prepared import audit_preparation,tasks,reuse,collect,validate
from mace_omol_intact import geometry

GLOBAL=ROOT/'workspaces/mace_global_benchmark_20260916/prepared_v1'
WORK=ROOT/'workspaces/mace_omol_20260917'
SOURCE=WORK/'charge_ablation_development_v2/collection_job_1200828.json'


@unittest.skipUnless(SOURCE.exists(),'requires actual completed charge-ablated molecular calculations')
class PreparedInterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:import openmm
        except ImportError:raise unittest.SkipTest('exact preparation replay requires the existing OpenMM preparation environment')
        cls.source=read_json(SOURCE);cls.parent=read_json(verify(cls.source['manifest']))

    @cached_file_checks
    def test_real_pqq_and_ggr_replay_and_four_state_reuse(self):
        for name in ('PQQ_1H4I','GGR_1GLG'):
            p=GLOBAL/name/'preparation.json';audit=audit_preparation(p)
            self.assertEqual(audit['status'],'pass');self.assertEqual(audit['new_model_forwards'],0)
            wanted=tasks(p);refs,values=reuse(SOURCE,wanted,self.parent)
            self.assertEqual(len(wanted),4);self.assertEqual(len(refs),4);self.assertEqual(len(values),4)
            for task in wanted:
                old=next(t for t in self.parent['tasks'] if t['task_id']==refs[task['task_id']]['source_task_id'])
                self.assertEqual(geometry(task),xyz(verify(old['xyz'])))
                self.assertEqual(task['charge'],old['charge']);self.assertEqual(task['state'],old['state'])

    @cached_file_checks
    def test_corrupted_mapping_cannot_pass_source_replay(self):
        original=read_json(GLOBAL/'GGR_1GLG/preparation.json')
        corrupted=copy.deepcopy(original);corrupted['physical_atoms'][0]['xyz_A'][0]+=.1
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_real_preparation.json';write_new(p,corrupted)
            with self.assertRaisesRegex(InvalidArtifact,'exact source replay'):audit_preparation(p)

    @cached_file_checks
    def test_actual_1kb0_false_peptide_connections_rejected(self):
        p=WORK/'intact_panel_prepared_v1/1KB0/preparation.json'
        if not p.exists():self.skipTest('requires actual archived invalid1KB0 preparation')
        with self.assertRaisesRegex(InvalidArtifact,'unsupported peptide connection'):audit_preparation(p)

    @cached_file_checks
    def test_real_reused_interface_score_has_no_automatic_class(self):
        p=WORK/'prepared_interface_pqq_v1/manifest.json'
        if not p.exists():self.skipTest('requires prepared exact-reuse interface manifest')
        self.assertEqual(validate(p)['tasks'],0)
        result=collect(p);original=self.source['scores']['PQQ_1H4I']
        self.assertEqual(result['R_mask_model_kcal'],original['R_mask_model_kcal'])
        self.assertTrue(result['numerical_gate_pass']);self.assertIsNone(result['calibrated_class'])
        self.assertEqual(result['decision_status'],'compatible_calibration_not_supplied')


if __name__=='__main__':unittest.main()
