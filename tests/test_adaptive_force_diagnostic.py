"""Geometry/algebra tests on actual archived forces; no molecular evaluations."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from adaptive_force_diagnostic import project, preview


class ArchivedForceDiagnostic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT/'workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json'
        cls.manifest = read_json(path)
        cls.tasks = [t for t in cls.manifest['tasks'] if t['case_id'] == '4MAE']
        cls.points = [read_json(path.parent/'proposals'/t['task_id']/'result.json')['origin'] for t in cls.tasks]
        cls.data = read_json(verify(cls.tasks[0]['mapping']))['context']
        cls.projections = [project(cls.data, p['full_q'], np.load(verify(p['forces']), allow_pickle=False)) for p in cls.points]

    def test_actual_active_derivatives_and_unit_norms(self):
        for task, point, result in zip(self.tasks, self.points, self.projections):
            _, v, raw, normalized, length = result
            np.testing.assert_allclose(raw[task['active_indices']], point['gradient_kcal_mol_rad'], atol=1e-8, rtol=0)
            np.testing.assert_allclose(np.linalg.norm(v, axis=1), 1, atol=1e-12, rtol=0)
            np.testing.assert_allclose(normalized*length, raw, atol=1e-12, rtol=0)

    def test_rigid_axis_rotation_does_not_change_preview(self):
        _, v, _, ca, _ = self.projections[0]
        la = self.projections[1][3]
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        moved = (v.reshape(len(v), -1, 3) @ rotation).reshape(v.shape)
        a = preview(self.data['modes'], v, ca, la)
        b = preview(self.data['modes'], moved, ca, la)
        self.assertEqual([r['id'] for r in a['selected']], [r['id'] for r in b['selected']])
        np.testing.assert_allclose([r['differential_kcal_mol_A'] for r in a['selected']],
                                   [r['differential_kcal_mol_A'] for r in b['selected']], atol=1e-11, rtol=0)

    def test_duplicate_real_coordinate_does_not_add_a_direction(self):
        _, v, _, ca, _ = self.projections[0]; la = self.projections[1][3]
        # Deliberately duplicate an existing real coordinate for an algebra test.
        original = preview(self.data['modes'], v, ca, la, maximum=99)
        i = next(i for i,m in enumerate(self.data['modes']) if m['kind'] == 'sidechain_torsion')
        duplicate = dict(self.data['modes'][i], id='zz_duplicate_for_algebra_test')
        copied = preview(self.data['modes']+[duplicate], np.vstack([v,v[i]]),
                         np.append(ca,ca[i]), np.append(la,la[i]), maximum=99)
        self.assertEqual(len(copied['selected']), len(original['selected']))
        self.assertIn(duplicate['id'], copied['redundant'])

    def test_truncated_real_force_array_is_rejected(self):
        point = self.points[0]; force = np.load(verify(point['forces']), allow_pickle=False)
        with self.assertRaisesRegex(InvalidArtifact, 'force dimensions'):
            project(self.data, point['full_q'], force[:-1])

    def test_seven_actual_composite_contexts_replay_in_common_measure(self):
        from adaptive_composite_diagnostic import analyze
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'comparison.json'
            analyze(ROOT/'workspaces/accommodation_response_20260920/result_v2.json',
                    ROOT/'diagnostics/adaptive_accommodation_20260922/PLAN.md', out)
            result = read_json(out)
            self.assertEqual(result['denominator'], 7)
            self.assertEqual(result['new_molecular_calls'], 0)
            self.assertIsNone(result['affinity_score'])
            self.assertEqual({r['case_id'] for r in result['cases']},
                             {'1H4I','4MAE','1F6S','6IP9','1GLG','2FW0','2FVY'})


if __name__ == '__main__': unittest.main()
