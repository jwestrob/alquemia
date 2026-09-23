"""Actual source identities and unit algebra; no synthetic scientific endpoints."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import adaptive_completion as completion
from affordable_common import HA_TO_KCAL,read_json,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL


class AdaptiveCompletion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=ROOT/'workspaces/adaptive_completion_20260922/original30_v1/manifest.json'
        cls.m=read_json(cls.p)
        cls.probe=read_json(ROOT/'workspaces/adaptive_completion_20260922/FIRST_STEPS_v1.json')

    def test_declared_all30_and_exact_eight_task_pilot(self):
        self.assertEqual(len(self.m['declared_case_ids']),30)
        self.assertEqual(len(self.m['tasks']),60)
        self.assertEqual(self.m['pilot_task_ids'],[c+'__'+z for c in completion.PILOT for z in ('Ca','La')])
        self.assertEqual(self.m['settings'],completion.SETTINGS)
        for t in self.m['tasks']:
            oldm=read_json(verify(t['prior_angular_manifest']))
            old=next(r for r in oldm['tasks'] if r['task_id']==t['task_id'])
            self.assertEqual({k:v for k,v in t.items() if k not in ('prior_angular_manifest','origin_reuse')},old)

    def test_consistent_positive_unit_scaling_on_actual_gradients(self):
        self.assertGreater(completion.SCALE,0)
        self.assertEqual(completion.SETTINGS['optimizer_ftol']/completion.SCALE,1e-9)
        for t in self.m['tasks']:
            g=np.array(t['origin_reuse']['point']['gradient_kcal_mol_rad'])
            np.testing.assert_allclose((g/EV_TO_KCAL)*completion.SCALE,g/HA_TO_KCAL,atol=1e-15,rtol=1e-14)

    def test_all_real_first_steps_retained_without_displaced_energies(self):
        self.assertEqual(self.probe['denominator'],60)
        self.assertEqual(self.probe['geometry_guard_pass'],60)
        self.assertEqual(self.probe['new_molecular_calls'],0)
        self.assertTrue(all(r['nonzero_objective_evaluations']==0 for r in self.probe['rows']))
        failed=next(r for r in self.probe['rows'] if r['task_id']=='mmol_1770-pqq-la_model__La')
        self.assertTrue(failed['physical_status']['physical_feasible'])
        self.assertTrue(failed['geometry_checks']['pass'])

    def test_exact_origin_reuse_needs_no_model_and_preserves_raw_energy_gradient(self):
        t=next(t for t in self.m['tasks'] if t['task_id']=='mmol_1770-pqq-la_model__La')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'actual_manifest_copy.json';write_new(p,self.m)
            ev=completion.CompletionProposal(p,t,None)
            actual=ev.evaluate(np.zeros(4),'real_origin_replay')
            old=t['origin_reuse']['point']
            for k in ('MACE_eV','gradient_kcal_mol_rad','forces','MACE','coordinate','full_q'):
                self.assertEqual(actual[k],old[k])
            self.assertEqual(xyz(verify(actual['coordinate'])),xyz(verify(t['xyz'])))
            self.assertEqual(ev.evaluate(np.zeros(4),'repeat')['MACE_eV'],old['MACE_eV'])
            self.assertTrue(all(r['scientific_origin_reused'] for r in ev.requests))


if __name__=='__main__':unittest.main()
