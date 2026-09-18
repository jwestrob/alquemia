"""Real native artifacts only; malformed cases explicitly damage those artifacts."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_native_field_accounting import parse,physical_identity,validate,accepted,task_key,check_runtime_parameters

W=ROOT/'workspaces/mace_omol_20260917/native_field_accounting_v1'
P=ROOT/'workspaces/mace_omol_20260917/tinker_framework_solver_v1'


class ExistingControls(unittest.TestCase):
    def test_static_compatibility_is_physical_not_mask_or_tolerance(self):
        if not (P/'manifest.json').exists():self.skipTest('real framework controls unavailable')
        m=read_json(P/'manifest.json')
        for case in {t['case_id'] for t in m['tasks']}:
            ts={t['variant']:t for t in m['tasks'] if t['case_id']==case}
            params={v:read_json(verify(t['parameters'])) for v,t in ts.items()}
            self.assertEqual(physical_identity(params['all_standard']),physical_identity(params['frozen_tight']))
            self.assertNotEqual(physical_identity(params['all_standard']),physical_identity(params['frozen_rigid']))
            bad=copy.deepcopy(params['all_standard']);bad['atoms'][0]['radius_A']+=.01
            self.assertNotEqual(physical_identity(bad),physical_identity(params['all_standard']))


class NativeAccounting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('accounting preparation unavailable')
        cls.m=validate(W/'manifest.json')

    def test_method_change_invalidates_real_task_key(self):
        t=self.m['tasks'][0];bad=copy.deepcopy(self.m);bad['tolerances']['accounting_kcal_mol']=1.
        self.assertNotEqual(task_key(t,self.m),task_key(t,bad))

    def test_query_and_static_results_and_damaged_records(self):
        items=[(t,accepted(t,W/'tasks'/t['task_id'])) for t in self.m['tasks']]
        recovery=W/'recovery_pinned_v1/recovery.json'
        if recovery.exists():
            recovered=read_json(recovery)['results']
            items=[(t,a or ((read_json(verify(recovered[t['task_id']])),None) if t['task_id'] in recovered else None)) for t,a in items]
        if any(not a for _,a in items):self.skipTest('eighteen real accounting executions incomplete')
        for t,(r,_) in items:
            text=verify(r['receipt']['log']).read_text();p=parse(text,t)
            check_runtime_parameters(p['parameters'],t,64)
            self.assertEqual(p['energies_kcal_mol'],r['energies_kcal_mol'])
            self.assertIsNone(p['numerical_score'])
            self.assertTrue(all(np.all(np.array(v[1:])==0) for v in p['response_eA'].values()))
            with self.assertRaisesRegex(InvalidArtifact,'incomplete'):
                parse(text.replace('ALQUEMIA_ACCOUNTING_COMPLETE','DAMAGED_COMPLETION'),t)
            # This is a corrupted copy of actual output, not simulated successful science.
            damaged='\n'.join(s for s in text.splitlines() if s.split()[:2]!=['RESPONSE','1'])
            with self.assertRaisesRegex(InvalidArtifact,'missing response'):
                parse(damaged,t)

    def test_real_contraction_agrees_with_native_components(self):
        reports=sorted(W.glob('collection_recovery_pinned_*.json'))
        if not reports:self.skipTest('native comparison not yet executed')
        r=read_json(reports[-1]);self.assertTrue(r['complete']);self.assertEqual(len(r['comparisons']),12)
        # A scientific mismatch is retained, never converted into fabricated agreement.
        for p in r['comparisons']:
            observed=p['accounting_residuals_kcal_mol']
            self.assertEqual(p['pass_'],all(abs(v)<=1e-7 for v in observed.values()))
        self.assertIsNone(r['numerical_score'])

    def test_frozen_coordinate_validation_rejects_actual_mutation(self):
        task=next(t for t in self.m['tasks'] if t['variant']=='frozen_rigid')
        params=read_json(verify(task['parameters']));params['inventory'][-1]=64
        check_runtime_parameters(params,task,64)
        params['atoms'][0]['xyz_A'][0]+=.001
        with self.assertRaisesRegex(InvalidArtifact,'runtime parameters differ'):
            check_runtime_parameters(params,task,64)


if __name__=='__main__':unittest.main()
