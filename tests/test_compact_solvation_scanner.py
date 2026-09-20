"""Pinned real preparations, partial execution and scalar algebra; no mock energies."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, write_new
import compact_solvation_scanner as scanner

INVENTORY = ROOT / 'diagnostics/compact_solvation_20260920/INVENTORY.json'
COMPARISON = ROOT / 'workspaces/compact_solvation_20260920/full_v1/comparison_v1.json'
MANIFEST = ROOT / 'workspaces/compact_scanner_20260920/pilot_v1/manifest.json'


class Scanner(unittest.TestCase):
    def test_finite_actual_preparations_and_checkpoint(self):
        result = scanner.validate(MANIFEST)
        self.assertEqual((result['native_MACE'], result['native_GFN2']), (8, 16))
        m = read_json(MANIFEST)
        self.assertEqual(m['cases'], list(scanner.CASES))
        mm = read_json(m['mace_manifest']['path'])
        self.assertTrue(all(t['energy_only'] and t['capture_native_readout'] for t in mm['tasks']))
        self.assertEqual(mm['model']['checkpoint']['sha256'], scanner.omol.CHECKPOINT_SHA)

    def test_fixed_boundary_has_no_repeat_tolerance_padding(self):
        bands = read_json(COMPARISON)['calibration']['context']['bands']
        ca = bands['Ca_supported_max_R_model_kcal_mol']
        self.assertEqual(scanner.decision(ca, bands), 'Ca-supported')
        self.assertEqual(scanner.decision(ca + 1e-7, bands), 'inconclusive')
        self.assertEqual(scanner.decision(bands['La_supported_min_R_model_kcal_mol'], bands), 'La-supported')
        self.assertIsNone(scanner.decision(None, bands))

    def test_generic_preparation_cannot_acquire_PQQ_bands_by_label_scope(self):
        cases = read_json(INVENTORY)['cases']
        pqq = next(c for c in cases if c['case_id'] == '1H4I')
        alpha = copy.deepcopy(next(c for c in cases if c['case_id'] == '1F6S'))
        self.assertTrue(scanner.pqq_decision_compatible(pqq))
        alpha['label_scope'] = 'PQQ_prediction'  # explicitly corrupted metadata on a real preparation
        self.assertFalse(scanner.pqq_decision_compatible(alpha))

    def test_gpu_venv_invocation_preserved_despite_python_symlink(self):
        m = read_json(MANIFEST)
        self.assertIn('/venv/bin/python', m['gpu_python_invocation'])
        self.assertNotEqual(m['gpu_python_invocation'], m['gpu_python']['path'])
        self.assertEqual(record(m['gpu_python_invocation']), m['gpu_python'])

    def test_explicit_corruption_of_real_state_rejected(self):
        m = read_json(MANIFEST)
        with tempfile.TemporaryDirectory() as directory:
            mm = read_json(m['mace_manifest']['path'])
            mm['tasks'][0]['charge'] += 1  # identified corrupted real fixture, not scientific evidence
            mp = Path(directory) / 'corrupted_real_mace_manifest.json'; write_new(mp, mm)
            m['mace_manifest'] = record(mp)
            path = Path(directory) / 'corrupted_real_scanner_manifest.json'; write_new(path, m)
            with self.assertRaises(InvalidArtifact):
                scanner.validate(path)

    def test_explicit_pair_needs_no_old_energy_and_partial_stays_null(self):
        inv = read_json(INVENTORY)
        inv['cases'] = [copy.deepcopy(next(c for c in inv['cases'] if c['case_id'] == '1H4I'))]
        for endpoint in inv['cases'][0]['representations']['context']['endpoints'].values():
            for key in ('native_MACE_energy_eV', 'native_MACE_receipt', 'source_manifest', 'source_task_id'):
                endpoint.pop(key)
        inv['cases'][0]['representations'].pop('core')
        m = read_json(MANIFEST)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'actual_1H4I_without_archived_energies.json'; write_new(path, inv)
            out = Path(directory) / 'prepared'
            result = scanner.prepare_pairs(path, COMPARISON, m['agreement']['path'], out,
                                           m['cpu_python_invocation'], m['gpu_python_invocation'])
            self.assertEqual((result['native_MACE'], result['native_GFN2']), (2, 4))
            result = scanner.collect(out / 'manifest.json', out / 'unrun_collection.json')
            self.assertEqual(result['available'], 0)
            self.assertIsNone(result['rows'][0]['composite_R_model_kcal_mol'])
            self.assertIsNone(result['rows'][0]['archived_composite_R_model_kcal_mol'])
            self.assertIsNone(result['rows'][0]['published_PQQ_decision'])

    def test_actual_scalar_success_does_not_hide_failed_MPI_components(self):
        actual = read_json(ROOT / 'workspaces/compact_scanner_20260920/pilot_v1/result_1203299.json')
        self.assertEqual(actual['available'], 0)
        for row in actual['rows']:
            self.assertIsNone(row['composite_R_model_kcal_mol'])
            self.assertIsNone(row['published_PQQ_decision'])
            for native in row['native_MACE'].values():
                self.assertEqual(native['status'], 'available')
                self.assertEqual(native['archived_difference_kcal_mol'], 0)
                self.assertEqual(native['analytic_gradient_status'], 'not_requested')
            for low in row['low_level_endpoints'].values():
                self.assertEqual(low['status'], 'unavailable')
                self.assertIsNone(low['delta_solv_hartree'])
                self.assertTrue(all(e['available_artifacts'] for e in low['endpoints'].values()))

    def test_complete_actual_integrated_repeat_and_unit_algebra(self):
        path = ROOT / 'workspaces/compact_scanner_20260920/pilot_v2/result_1203319.json'
        actual = read_json(path)
        self.assertEqual((actual['available'], actual['case_denominator']), (4, 4))
        for row in actual['rows']:
            native = {z: e['energy_eV'] for z, e in row['native_MACE'].items()}
            low = row['low_level_endpoints']
            recomputed = scanner.mix_pair(native, {z: e['vacuum_hartree'] for z, e in low.items()},
                                          {z: e['alpb_hartree'] for z, e in low.items()})
            self.assertEqual(recomputed['composite_R_model_kcal_mol'], row['composite_R_model_kcal_mol'])
            self.assertTrue(row['repeat_pass'])
            self.assertTrue(all(e['archived_difference_kcal_mol'] == 0 for e in row['native_MACE'].values()))
            self.assertTrue(all(e['normal_termination_confirmed'] for metal in low.values() for e in metal['endpoints'].values()))
        self.assertGreater(actual['alpha_minus_GGR_model_kcal_mol'], 0)
        self.assertEqual([r['published_PQQ_decision'] for r in actual['rows']], ['Ca-supported', 'Ca-supported', None, None])


if __name__ == '__main__':
    unittest.main()
