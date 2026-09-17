"""Pinned real OMOL preparations and receipts; no mocked successful inference."""
import copy
from pathlib import Path
import sys
import tempfile
import json
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_omol import validate,COMPONENT,accepted_state
from mace_omol_report import paired_energy
W=ROOT/'workspaces/mace_omol_20260917'
MANIFEST=W/'qualification_v1/manifest.json'


@unittest.skipUnless(MANIFEST.exists(),'requires actual OMOL qualification preparation')
class OMOLTests(unittest.TestCase):
    def test_exact_qualification_inventory_and_declared_model(self):
        self.assertEqual(validate(MANIFEST)['status'],'pass')
        m=read_json(MANIFEST)
        self.assertEqual(len(m['tasks']),10)
        self.assertTrue(all(t['energy_component']==COMPONENT for t in m['tasks']))
        self.assertEqual(m['model']['head'],'omol')
        self.assertIsNone(m['model']['solvent'])
        self.assertEqual({t['variant'] for t in m['tasks']},{'primary','repeat','rotate','translate'})

    def test_real_task_charge_and_checkpoint_corruption_rejected(self):
        original=read_json(MANIFEST)
        with tempfile.TemporaryDirectory() as directory:
            for kind in ('charge','checkpoint'):
                m=copy.deepcopy(original)
                if kind=='charge':m['tasks'][0]['charge']+=1
                else:m['model']['checkpoint']['sha256']='0'*64
                path=Path(directory)/(kind+'.json');path.write_text(json.dumps(m))
                with self.assertRaises(InvalidArtifact):validate(path)

    def test_actual_OMOL_receipt_cannot_claim_predicted_density_or_wrong_state(self):
        paths=sorted((W/'qualification_v1/execution').glob('*/attempt_*/result.json'))
        if not paths:self.skipTest('actual scientific OMOL output not yet available')
        r=next((read_json(p) for p in paths if read_json(p)['status']=='computed'),None)
        if r is None:self.skipTest('no successful scientific OMOL output exists')
        m=read_json(verify(r['manifest']));task=next(t for t in m['tasks'] if t['task_id']==r['task_id'])
        self.assertTrue(accepted_state(r,task));self.assertIsNone(r['density_coefficients'])
        bad=copy.deepcopy(r);bad['input_state_check']['spin_multiplicity']=0
        self.assertFalse(accepted_state(bad,task))
        bad=copy.deepcopy(r);bad['density_coefficients']={'corrupted':'actual output copy'}
        self.assertFalse(accepted_state(bad,task))

    @unittest.skipUnless((W/'benchmark_v1/manifest.json').exists(),'requires prepared benchmark')
    def test_benchmark_preserves_all_cases_and_only_qualified_cache(self):
        m=read_json(W/'benchmark_v1/manifest.json')
        self.assertEqual(validate(W/'benchmark_v1/manifest.json')['status'],'pass')
        self.assertEqual(len(m['tasks']),60)
        self.assertEqual(set(m['reused']),{'1H4I_Ca','1H4I_La','4MAE_Ca','4MAE_La'})
        self.assertEqual(sum(t['evaluation_role']=='retrospective_nonPQQ_direction' for t in m['tasks']),8)
        self.assertEqual({t['case_id'] for t in m['tasks'] if t['evaluation_role']=='retrospective_nonPQQ_direction'},
                         {'GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9'})

    @unittest.skipUnless((W/'qualification_v1/collection_job_1200797.json').exists(),'requires executed endpoints')
    def test_actual_energy_contrast_and_missing_endpoint_status(self):
        c=read_json(W/'qualification_v1/collection_job_1200797.json')
        ca=c['rows']['1H4I_Ca_primary'];la=c['rows']['1H4I_La_primary']
        self.assertAlmostEqual(paired_energy(ca,la),(ca['energy_eV']-la['energy_eV'])*23.06054783061903,places=9)
        unavailable=copy.deepcopy(la);unavailable.update(status='unavailable',energy_eV=None)
        self.assertIsNone(paired_energy(ca,unavailable))

    @unittest.skipUnless((W/'report_v1/result.json').exists(),'requires completed real benchmark report')
    def test_actual_PQQ_success_does_not_promote_failed_primary_affinity(self):
        r=read_json(W/'report_v1/result.json')
        self.assertEqual(r['calibration']['valid_calibration_count'],25)
        self.assertTrue(r['canonical_operational_gate_pass'])
        self.assertEqual(r['supported_transfer_count'],3)
        primary=[x for x in r['nonPQQ_contrasts'] if x['role']=='primary']
        self.assertEqual(len(primary),2)
        self.assertTrue(all(x['negative_case']=='GGR_extended' and not x['pass'] for x in primary))
        self.assertFalse(r['nonPQQ_primary_gate_pass'])
        self.assertFalse(r['nonPQQ_robustness_gate_pass'])
        self.assertFalse(r['broad_affinity_validated'])
        self.assertIsNone(r['solution_score'])


if __name__=='__main__':unittest.main()
