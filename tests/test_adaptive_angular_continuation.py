"""Real source/map continuation checks; no model execution or synthetic energies."""
import ast
import copy
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import adaptive_angular_proposals as adaptive
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_site_kinematics import Kinematics


class AdaptiveAngularContinuation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old_path = ROOT / 'workspaces/adaptive_accommodation_20260922/proposals_v1/manifest.json'
        cls.path = ROOT / 'workspaces/adaptive_accommodation_20260922/remaining26_proposals_v1/manifest.json'
        cls.old = read_json(cls.old_path)
        cls.m = read_json(cls.path)
        cls.parent = read_json(verify(cls.m['source_manifest']))

    def test_fixed_complement_all_52_real_endpoints(self):
        expected = [c['case_id'] for c in self.parent['cases'] if c['case_id'] not in adaptive.PILOT4]
        self.assertEqual(self.m['declared_case_ids'], expected)
        self.assertEqual(len(expected), 26)
        self.assertEqual(self.m['new_optimizer_starts'], 52)
        self.assertEqual({t['case_id'] for t in self.m['tasks']} & set(adaptive.PILOT4), set())
        checked = adaptive.validate(self.path)
        self.assertEqual((checked['cases'], checked['tasks'], checked['molecular_calls']), (26, 52, 0))
        # Crystal 1KB0 retains the parent's same functional-class label scope;
        # source IDs, not that shared label field, identify the fixed split.
        self.assertEqual(sum(c['case_id'].endswith('-pqq-la_model') for c in self.m['cases']), 25)
        self.assertEqual([c['case_id'] for c in self.m['cases'] if not c['case_id'].endswith('-pqq-la_model')], ['1KB0'])

    def test_original_snapshot_and_physics_unchanged(self):
        self.assertEqual(record(self.old_path)['sha256'], 'b5c4218caee58458cc98ee48b0bbcf165bebb88615b08750905de8b66fb5a998')
        self.assertEqual(self.m['settings'], self.old['settings'])
        self.assertEqual(self.m['model'], self.old['model'])
        self.assertEqual(adaptive.validate(self.old_path)['tasks'], 8)
        old_code = Path(verify(self.old['implementation']['adaptive_angular_proposals.py'])).read_text()
        new_code = Path(adaptive.__file__).read_text()
        def functions(text):
            return {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(text).body if isinstance(n, ast.FunctionDef)}
        a, b = functions(old_code), functions(new_code)
        for name in ['full_q', 'constraints', 'physical_status', 'final_geometry', 'optimize']:
            self.assertEqual(a[name], b[name], name)

    def test_default_prepare_replays_original_tasks_and_settings(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / 'default_pilot_replay'
            adaptive.prepare(source=verify(self.old['source_manifest']), diagnostic=verify(self.old['diagnostic']),
                             agreement=verify(self.old['agreement']), output=dest)
            repeated = read_json(dest / 'manifest.json')
            self.assertEqual({k: v for k, v in repeated.items() if k != 'implementation'},
                             {k: v for k, v in self.old.items() if k != 'implementation'})

    def test_remaining_maps_pair_and_small_motion_constraints(self):
        for cid in self.m['declared_case_ids']:
            ca, la = [next(t for t in self.m['tasks'] if t['case_id'] == cid and t['metal'] == z) for z in ('Ca', 'La')]
            self.assertEqual(ca['active_mode_ids'], la['active_mode_ids'])
            self.assertEqual(read_json(verify(ca['mapping'])), read_json(verify(la['mapping'])))
            self.assertEqual(xyz(verify(ca['xyz']))[1:], xyz(verify(la['xyz']))[1:])
            kin = Kinematics(read_json(verify(ca['mapping']))['context'])
            checked = adaptive.final_geometry(kin, ca, np.array([.002, -.001, .0015, -.002]),
                                               [a[0] for a in xyz(verify(ca['xyz']))])
            self.assertTrue(checked['physical_feasible'])
            self.assertTrue(checked['all_other_physical_atoms_fixed'])

    def test_constraint_derivatives_on_fixed_first_middle_last_sources(self):
        q = np.array([.01, -.008, .005, -.003]); step = 1e-6
        for index in [0, 12, 25]:
            cid = self.m['declared_case_ids'][index]
            t = next(t for t in self.m['tasks'] if t['case_id'] == cid and t['metal'] == 'Ca')
            kin = Kinematics(read_json(verify(t['mapping']))['context'])
            _, analytic = adaptive.constraints(kin, t, q)
            numeric = np.column_stack([(adaptive.constraints(kin, t, q + np.eye(4)[j]*step)[0] -
                                        adaptive.constraints(kin, t, q - np.eye(4)[j]*step)[0])/(2*step) for j in range(4)])
            np.testing.assert_allclose(analytic, numeric, atol=1e-8, rtol=1e-6)

    def test_actual_unrun_52_collection_remains_unavailable(self):
        # This real prelaunch collection remains a fixture after a later launch.
        result = read_json(self.path.parent / 'collection_unrun_v1.json')
        self.assertEqual(result['manifest'], record(self.path))
        self.assertEqual((result['case_denominator'], result['endpoint_denominator'], result['available_candidates']), (26, 52, 0))
        self.assertTrue(all(e['reason'] == 'not_run' and e['candidate'] is None and
                            e['composite_candidate_energy'] is None for e in result['endpoints']))
        self.assertTrue(all(c['score'] is None and c['status'] == 'unavailable' for c in result['cases']))

    def test_runtime_gate_still_requires_actual_original_four_pool(self):
        base = ROOT / 'workspaces/nikasha_shared_pool_20260922/pilot_v2'
        completed = base / 'after_solvent_0_1209845.json'
        self.assertEqual(adaptive.pool_gate(self.path, completed), record(completed))
        with self.assertRaisesRegex(InvalidArtifact, 'complete compatible four-case'):
            adaptive.pool_gate(self.path, base / 'after_MACE_1209840.json')
        expanded = ROOT / 'workspaces/adaptive_accommodation_20260922/common_pool_v1'
        # Do not treat an adaptive expansion as the required original protocol.
        manifest = read_json(expanded / 'manifest.json')
        self.assertNotEqual(manifest['protocol_id'], read_json(completed)['protocol_id'])

    def test_corrupted_real_population_and_task_count_are_rejected(self):
        for field, value in [('population', 'pilot4'), ('excluded_case_ids', []), ('new_optimizer_starts', 8)]:
            broken = copy.deepcopy(self.m); broken[field] = value
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / 'explicitly_corrupted_real_manifest.json'
                write_new(path, broken)
                with self.assertRaisesRegex(InvalidArtifact, 'scope differs|exclusions differ'):
                    adaptive.validate(path)


if __name__ == '__main__':
    unittest.main()
