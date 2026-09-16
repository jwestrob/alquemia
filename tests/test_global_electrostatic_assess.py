"""Real-artifact assessor tests; no fabricated scientific output or solver call."""
from pathlib import Path
import copy
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from affordable_tabi import transform_definition, transformed, physical_xyzr
import affordable_tabi as tabi
from global_electrostatic import TABI_COULOMB_KCAL_A
from global_electrostatic_assess import direct_coulomb, collect, assess, report, expected_tasks, _serialized_coordinates


class ArchivedDirectArithmetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / 'workspaces/affordable_challenger_20260915/solver_completion/solver_result.json'
        if not path.is_file():
            raise unittest.SkipTest('real archived charge/geometry arithmetic fixtures unavailable')
        verify({'path': str(path), 'sha256': '542d491159111114e54ba81df220b08f7bbf7d8db619e2a32b970dea3b0ee885'})
        cls.fixtures = []
        for row in read_json(path)['numerical_checks']:
            if row['label'].endswith('/primary'):
                result = row['result']
                manifest = read_json(verify(result['source_manifest']))
                state = read_json(verify(manifest['state']))
                cls.fixtures.append((state, result['components']['direct_core_environment_kcal_mol']))
        if len(cls.fixtures) != 4:
            raise AssertionError('archived four-state fixture set changed')

    def test_direct_sum_reproduces_archived_interaction_after_declared_constant_conversion(self):
        # Only arithmetic: archived CPCM charges are not a new gas/TABI result.
        for state, archived in self.fixtures:
            got = direct_coulomb(state, transform_definition('identity'))
            expected = archived * TABI_COULOMB_KCAL_A / 332.063713299
            self.assertAlmostEqual(got, expected, places=8)

    def test_rigid_transform_preserves_actual_archived_direct_interaction(self):
        for state, _ in self.fixtures:
            reference = direct_coulomb(state, transform_definition('identity'))
            for name in ('translated', 'rotated'):
                transformed = direct_coulomb(state, transform_definition(name))
                self.assertLess(abs(transformed - reference), 1e-8)


class ActualSurfaceCampaign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / 'workspaces/global_electrostatic_20260916/surfaces_v1/campaign_manifest.json'
        if not path.is_file():
            raise unittest.SkipTest('actual frozen surface campaign not yet prepared')
        cls.collection = collect(path)
        cls.assessment = assess(cls.collection)

    def test_actual_schedule_and_state_provenance_are_complete(self):
        self.assertEqual({r['task_id'] for r in self.collection['rows']}, set(expected_tasks()))
        self.assertEqual(len(self.collection['rows']), 25)
        self.assertTrue(all(c['status'] == 'passed' for c in self.collection['state_checks']))
        self.assertEqual(self.collection['quantum']['status'], 'endpoints_complete')

    def test_current_collector_accepts_actual_frozen_adapter_snapshot(self):
        row = next(r for r in self.collection['rows'] if r['task_id'] == '1h4i_qm33_La/primary')
        manifest, _ = tabi._verified_manifest(verify(row['manifest']))
        frozen = manifest['implementation']['affordable_tabi.py']
        self.assertEqual(frozen['sha256'], '50ed8e24668c3777c66f34c888208548332e357769171f248e3bc3dcec0f2454')
        verify(frozen)
        self.assertNotEqual(frozen['sha256'], tabi.LOADED_IMPLEMENTATION['sha256'])

    def test_corrupted_real_collection_does_not_substitute_a_score(self):
        corrupt = copy.deepcopy(self.collection)
        row = next(r for r in corrupt['rows'] if r['task_id'] == '1h4i_qm33_La/primary')
        row.update(status='failed', solver=None, endpoint=None, reason='explicitly corrupted real collection for failure test')
        result = assess(corrupt)
        self.assertEqual(result['gate'], 'failed')
        self.assertFalse(result['stage2_eligible'])
        self.assertIsNone(result['paired_contrasts']['1h4i_qm33/primary']['R_global_kcal_mol'])
        self.assertIsNone(result['calibrated_decision'])

    def test_corrupted_real_collection_rejects_duplicate_rows(self):
        corrupt = copy.deepcopy(self.collection)
        duplicate = copy.deepcopy(corrupt['rows'][0])
        duplicate.update(status='failed', solver=None, endpoint=None,
                         reason='explicit duplicate-row corruption of real collection')
        # Place the failed duplicate first: the old dictionary conversion could
        # silently replace it with the original later row of the same task ID.
        corrupt['rows'].insert(0, duplicate)
        with self.assertRaisesRegex(InvalidArtifact, 'exactly the frozen 25 task rows'):
            assess(corrupt)
        # Also reject a repeated ID that replaces one task without adding rows.
        corrupt = copy.deepcopy(self.collection)
        corrupt['rows'][1] = copy.deepcopy(corrupt['rows'][0])
        with self.assertRaisesRegex(InvalidArtifact, 'frozen task schedule'):
            assess(corrupt)

    def test_all_partition_resolutions_are_required(self):
        names = {c['name'] for c in self.assessment['physical_checks']}
        self.assertTrue({'partition/primary', 'partition/refined', 'partition/tree'} <= names)
        self.assertTrue({'common_mesh/primary/identity', 'common_mesh/refined/identity',
                         'common_mesh/primary/translated', 'common_mesh/primary/rotated',
                         'common_mesh/isolated_qm33/identity'} <= names)

    def test_real_completed_solver_components_are_retained_and_not_double_counted(self):
        if self.collection['counts']['computed'] != 25:
            self.skipTest('25 actual TABI tasks not yet successfully collected; no fake integration result')
        for row in self.collection['rows']:
            if row['endpoint'] is None:
                continue
            e = row['endpoint']
            expected = e['vacuum_energy_hartree'] * 627.509474 + row['direct_coulomb_kcal_mol'] + row['solver']['reaction_field_kJ_mol'] / 4.184
            self.assertEqual(e['total_kcal_mol'], expected)
            self.assertFalse(row['solver']['coulomb_added_to_descriptor'])
        for row in self.assessment['reaction_field_components'].values():
            self.assertEqual(row['status'], 'computed')
            self.assertAlmostEqual(row['reaction_field_cross_kcal_mol'] + row['reaction_field_core_kcal_mol'] +
                                   row['reaction_field_environment_kcal_mol'], row['reaction_field_total_kcal_mol'], places=8)

    def test_actual_report_is_explicit_about_scope_and_missing_values(self):
        text = report(self.assessment)
        self.assertIn('Physical gate:', text)
        self.assertIn('Baseline remains default', text)
        self.assertIn('launches no conditional Stage 2 jobs', text)
        if self.collection['counts']['unavailable']:
            self.assertIn('unavailable', text)
        self.assertFalse(self.assessment['accuracy_trial_executed_by_assessor'])

    def test_real_rotated_input_roundoff_is_bounded_and_coordinate_corruption_fails(self):
        row = next(r for r in self.collection['rows'] if r['task_id'] == '1h4i_qm33_La/rotated')
        manifest = read_json(verify(row['manifest']))
        state = read_json(verify(manifest['state']))
        expected = physical_xyzr(transformed(state['physical_atoms'], manifest['transform']))
        actual = verify(manifest['physical_xyzr']).read_text()
        self.assertLessEqual(_serialized_coordinates(actual, expected, {0, 1, 2}), 1.01e-10)
        # Explicitly corrupted copy of actual input, not a new physical state.
        rows = actual.splitlines(); fields = rows[0].split(); fields[0] = str(float(fields[0]) + .001)
        rows[0] = ' '.join(fields)
        with self.assertRaises(InvalidArtifact):
            _serialized_coordinates('\n'.join(rows) + '\n', expected, {0, 1, 2})


if __name__ == '__main__':
    unittest.main()
