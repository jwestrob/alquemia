"""Real pinned fold-source preparation and missing-status checks; no model calls."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import adaptive_completion as completion
import adaptive_completion_folds as folds
from adaptive_force_diagnostic import project
from affordable_common import read_json,verify


class AdaptiveCompletionFolds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=ROOT/'workspaces/adaptive_completion_20260922/primary225_v1'
        cls.m=read_json(cls.root/'manifest.json')
        cls.c=read_json(cls.root/'collection_unrun_v1.json')

    def test_actual_population_replays_full_preflight(self):
        result=folds.validate(self.root/'manifest.json')
        self.assertEqual((result['cases'],result['prepared_cases'],result['tasks'],result['endpoint_denominator']),
                         (225,205,410,450))
        self.assertEqual(result['new_molecular_calls'],0)
        original=read_json(ROOT/'workspaces/adaptive_completion_20260922/original30_v1/manifest.json')
        self.assertEqual(self.m['settings'],original['settings'])
        self.assertEqual(self.m['protocol_id'],original['protocol_id'])
        self.assertIs(completion.optimize,folds.completion.optimize)

    def test_unrun_and_inherited_missing_never_become_success(self):
        self.assertEqual(len(self.c['cases']),225)
        self.assertEqual(len(self.c['endpoints']),450)
        self.assertEqual(self.c['available_candidates'],0)
        self.assertTrue(all(c['status']=='unavailable' for c in self.c['cases']))
        self.assertTrue(all(r['candidate'] is None and r['composite_candidate_energy'] is None
                            for r in self.c['endpoints']))
        unavailable={c['case_id'] for c in self.m['cases'] if c['status']!='prepared'}
        self.assertEqual(len(unavailable),20)
        self.assertTrue(all(t['case_id'] not in unavailable for t in self.m['tasks']))
        self.assertEqual(sum(r['reason']=='not_run' for r in self.c['endpoints']),410)
        self.assertEqual(sum(r['case_id'] in unavailable for r in self.c['endpoints']),40)
        self.assertIn('a8r3s4-pqq-la_model__conditioned_Ca__seed-1_sample-3',unavailable)

    def test_reprojected_four_mode_gradients_use_actual_unchanged_q0_forces(self):
        # All410 source endpoints use their same real forces, not old1/2-mode vectors.
        for t in self.m['tasks']:
            source=t['origin_reuse']['source_point'];point=t['origin_reuse']['point']
            for key in ('coordinate','forces','MACE','MACE_eV','full_q'):
                self.assertEqual(point[key],source[key])
            mapping=read_json(verify(t['mapping']))['context']
            forces=np.load(verify(point['forces']),allow_pickle=False)
            _,_,raw,_,_=project(mapping,point['full_q'],forces)
            np.testing.assert_array_equal(raw[t['active_indices']],point['gradient_kcal_mol_rad'])
            self.assertEqual(len(t['active_indices']),4)
            self.assertEqual(point['active_q_radian'],[0.]*4)


if __name__=='__main__':unittest.main()
