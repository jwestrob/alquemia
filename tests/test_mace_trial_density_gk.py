"""Real complete trial-density native model and matched component comparison."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,verify
from mace_density_gk_hybrid import validate,parse
from mace_density_gk_compare import compare

BASE=ROOT/'workspaces/mace_omol_20260917'
REPORT=BASE/'trial_density_gk_hybrid_v1/collection_job_1201137.json'
PARENT=BASE/'density_gk_hybrid_v1/collection_job_1201074.json'


class ActualTrialHybrid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not REPORT.exists():raise unittest.SkipTest('real responsive-trial GK calculation incomplete')
        cls.r=read_json(REPORT);cls.m=validate(verify(cls.r['manifest']))

    def test_new_calls_and_reuses_have_separate_actual_receipts(self):
        self.assertTrue(self.r['complete'])
        self.assertEqual((self.r['actual_energy_calls'],self.r['actual_field_queries'],self.r['actual_response_solves']),(34,32,40))
        self.assertEqual(self.r['reused_environment_components'],37)
        for mode in ('static','identity','solve'):
            t=next(t for t in self.m['static_tasks']+self.m['response_tasks'] if t['mode']==mode and t['state']=='environment')
            row=read_json(verify(self.r['tasks'][t['task_id']]['result']))
            old=read_json(verify(row['reused_from']))
            self.assertEqual(row['receipt'],old['receipt'])
            self.assertTrue(all(row[k] is False for k in ('energy_attempted','field_attempted','response_attempted')))
            self.assertEqual(parse(verify(row['receipt']['log']).read_text(),t)['energies_kcal_mol'],row['energies_kcal_mol'])

    def test_core_term_removes_old_field_once_and_retains_response_energy(self):
        a=read_json(verify(self.m['sources']['trial_source_audit']))
        for case,c in self.r['variants']['primary']['cases'].items():
            for metal,e in c['endpoints'].items():
                source=a['rows'][case+'_'+metal]
                expected=source['embedded_energy_hartree']*HA_TO_KCAL-source['old_field_coupling_kcal']
                self.assertAlmostEqual(e['components']['DFT_intrinsic_trial_kcal'],expected,places=7)
                self.assertNotIn('DFT_vacuum_kcal',e['components'])
        self.assertIsNone(self.r['reference']);self.assertIsNone(self.r['calibrated_class'])
        self.assertIsNone(self.r['combined_gradient']);self.assertFalse(self.r['baseline_changed'])

    def test_matching_change_algebra_and_explicit_missing_component_rejection(self):
        with tempfile.TemporaryDirectory() as d:
            result=compare(PARENT,REPORT,Path(d)/'comparison.json')
            self.assertEqual(len(result['cases']),4);self.assertEqual(len(result['independent_biological_groups']),2)
            self.assertTrue(all(abs(c['algebra_error_kcal'])<=1e-7 for c in result['cases'].values()))
            damaged=copy.deepcopy(self.r)
            del damaged['variants']['primary']['cases']['GGR_extended']['components_R_kcal']['DFT_intrinsic_trial_kcal']
            p=Path(d)/'explicitly_corrupted_real_report.json';p.write_text(json.dumps(damaged))
            with self.assertRaisesRegex(InvalidArtifact,'missing/incompatible'):
                compare(PARENT,p,Path(d)/'must_not_be_written.json')
            self.assertFalse((Path(d)/'must_not_be_written.json').exists())

    def test_actual_response_requires_unrounded_convergence_evidence(self):
        t=next(t for t in self.m['response_tasks'] if t['state']=='La' and t['variant']=='primary' and t['poleps']==1e-9)
        row=read_json(verify(self.r['tasks'][t['task_id']]['result']))
        text=verify(row['receipt']['log']).read_text()
        self.assertEqual(parse(text,t)['exact_RMS_Debye'],row['exact_RMS_Debye'])
        damaged='\n'.join(line for line in text.splitlines() if not line.strip().startswith('EXACT_FINAL_RMS_DEBYE'))
        with self.assertRaisesRegex(InvalidArtifact,'invalid/nonconverged'):
            parse(damaged,t)


if __name__=='__main__':unittest.main()
