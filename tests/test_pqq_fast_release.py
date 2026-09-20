"""Real PQQ source reconstruction and finite release scope, without dummy science."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, write_new, xyz
import pqq_fast_prepare as prep
import pqq_standard as standard
from compact_solvation_scanner import TOLERANCES


class FastSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = read_json(ROOT / 'diagnostics/pqq_fast_release_20260920/SOURCES.json')
        cls.fixed = prep.helpers(cls.source['config'])[0]

    def test_actual_full_panel_reconstructs_same_context_from_protonated_source(self):
        actual = read_json(ROOT / 'workspaces/pqq_fast_release_20260920/core_replay_v2/RESULT.json')
        self.assertEqual((actual['passed'], len(actual['rows'])), (28, 28))
        self.assertEqual(actual['new_scientific_calls'], 0)
        self.assertFalse(actual['fresh_protonation'])
        for row in actual['rows']:
            for comparison in row['comparisons'].values():
                self.assertTrue(comparison['exact'])
                self.assertEqual(comparison['max_displacement_A'], 0)

    def test_pinned_source_coverage_and_roles(self):
        cases = self.source['cases']
        self.assertEqual(len(cases), 28)
        self.assertEqual(sum(c['role'] == 'calibration' for c in cases), 25)
        self.assertEqual({c['case_id'] for c in cases if c['role'] != 'calibration'}, {'1H4I', '4MAE', '1KB0'})
        for c in cases:
            verify(c['source_structure']); verify(c['archived_core'])
            self.assertEqual(set(c['roles']), set(self.fixed.ROLE_ORDER))

    def test_real_core_rebuild_rejects_corrupted_donor_identity(self):
        source = copy.deepcopy(self.source['cases'][0])
        parent = read_json(verify(source['archived_core']))
        source['roles']['anchor_glutamate']['resname'] = 'VAL'  # explicitly corrupted real selector
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises((InvalidArtifact, ValueError)):
                prep.build_core(source, verify(source['archived_protonated']), verify(parent['protonation_manifest']),
                                verify(source['archived_normalized']), Path(d), self.source['config'], self.fixed)

    def test_actual_fresh_source_reproduction_includes_all_hydrogens_and_charges(self):
        actual = read_json(ROOT / 'workspaces/pqq_fast_release_20260920/full_v1/source_preparation/preparation.json')
        self.assertEqual((actual['supported'], actual['denominator']), (28, 28))
        self.assertFalse(actual['preparation_reused'])
        self.assertGreater(actual['fresh_preparation_seconds'], 0)
        for case in actual['cases']:
            self.assertTrue(case['normalized_bytes_equal_archive'])
            self.assertFalse(case['preparation_reused'])
            for comparison in case['context_comparison'].values():
                self.assertTrue(comparison['exact'])
                self.assertTrue(comparison['charge_identical'])
                self.assertEqual(comparison['changed_atom_indices'], [])

    def test_auto_PQQ_and_explicit_DFT_select_different_backends(self):
        self.assertEqual(standard.choose_backend(self.source), 'fast_pqq')
        self.assertEqual(standard.choose_backend(self.source, 'fast-pqq'), 'fast_pqq')
        self.assertEqual(standard.choose_backend(self.source, 'dft-reference'), 'dft_reference')

    def test_actual_existing_water_request_keeps_baseline_route(self):
        request = read_json(ROOT / 'workspaces/water_promotion_20260919/request_v1.json')
        self.assertEqual(standard.choose_backend(request), 'baseline_water')
        self.assertEqual(standard.choose_backend(request, 'dft-reference'), 'baseline_water')
        with self.assertRaises(InvalidArtifact): standard.choose_backend(request, 'fast-pqq')

    def test_unknown_source_protocol_is_not_silent_DFT_fallback(self):
        corrupted = copy.deepcopy(self.source)
        corrupted['schema_version'] = 'unsupported_source_protocol'  # corrupted real request
        with self.assertRaises(InvalidArtifact): standard.choose_backend(corrupted)

    def test_invalid_case_path_rejected_before_source_preparation(self):
        corrupted = copy.deepcopy(self.source)
        corrupted['cases'][0]['case_id'] = '../outside'  # corrupted real manifest
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'request.json'; write_new(path, corrupted)
            with self.assertRaises(InvalidArtifact): prep.prepare(path, Path(d) / 'output')
            self.assertFalse((Path(d) / 'output').exists())

    def test_full_fresh_classification_preserves_exact_band_edges(self):
        actual = read_json(ROOT / 'workspaces/pqq_fast_release_20260920/full_v1/RESULT_1203463.json')
        self.assertTrue(actual['promotion_gate'])
        self.assertEqual((actual['calibration_correct'], actual['transfer_correct']), (25, 3))
        self.assertFalse(actual['threshold_refitted'])
        self.assertFalse(actual['source_reused']); self.assertFalse(actual['energy_reused'])
        for row in actual['rows']:
            self.assertTrue(row['classification_pass'])
            self.assertLess(abs(row['score_difference_model_kcal_mol']),
                            TOLERANCES['composite_repeat_model_kcal_mol'])

    def test_explicit_DFT_manifest_retains_original_core_and_recipe(self):
        source = ROOT / 'workspaces/pqq_fast_release_20260920/full_v1/source_preparation/preparation.json'
        with tempfile.TemporaryDirectory() as d:
            result = standard.prepare_dft(source, Path(d) / 'dft', standard.DEFAULT_RELEASE)
            self.assertEqual((result['status'], result['tasks']), ('dry_run_pass', 56))
            manifest = read_json(verify(result['manifest']))
            prepared = read_json(source)
            for task in manifest['tasks']:
                case = next(c for c in prepared['cases'] if c['case_id'] == task['case'])
                old = case['core']['endpoints'][task['metal']]
                self.assertEqual(xyz(verify(task['xyz'])), xyz(verify(old['xyz'])))
                self.assertEqual(verify(task['input']).read_bytes(), verify(old['input']).read_bytes())

    def test_release_gate_and_implementation_are_checked(self):
        r = standard.release(standard.DEFAULT_RELEASE)
        corrupted = copy.deepcopy(r); corrupted['promoted'] = False
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'release.json'; write_new(path, corrupted)
            with self.assertRaises(InvalidArtifact): standard.release(path)

    def test_literal_standard_CLI_used_fresh_source_and_components(self):
        run = ROOT / 'workspaces/pqq_fast_release_20260920/standard_1H4I_v1'
        result = read_json(run / 'result_1203729.json')
        receipt = read_json(run / 'execution.json')
        self.assertEqual(result['available'], 1)
        self.assertEqual(result['rows'][0]['published_PQQ_decision'], 'Ca-supported')
        self.assertTrue(receipt['source_preparation_included'])
        self.assertFalse(receipt['folding_included'])
        self.assertGreater(receipt['source_to_score_seconds'], 0)
        request = read_json(verify(read_json(run / 'plan.json')['request']))
        self.assertFalse(any(k.startswith('archived_') for k in request['cases'][0]))
        reference = next(c for c in self.source['cases'] if c['case_id'] == '1H4I')
        fresh = read_json(run / 'source_preparation/preparation.json')['cases'][0]
        for metal in ('Ca', 'La'):
            self.assertEqual(xyz(verify(fresh['representations']['context']['endpoints'][metal]['xyz'])),
                             xyz(verify(reference['archived_context']['endpoints'][metal]['xyz'])))

    def test_prediction_request_needs_no_label_and_retains_unknown_domain(self):
        request = read_json(ROOT / 'diagnostics/pqq_fast_release_20260920/examples/1H4I_source.json')
        # Real source, but remove its entry from a test-only domain registry.
        # This tests metadata dispatch without inventing a structure or energy.
        r = standard.release(standard.DEFAULT_RELEASE)
        r['reference_source_sha256'] = []
        case = request['cases'][0]
        case.pop('expected_class'); case['role'] = 'prediction'; case['label_scope'] = 'PQQ_prediction'
        with tempfile.TemporaryDirectory() as d:
            d = Path(d); rp = d / 'release.json'; qp = d / 'request.json'
            write_new(rp, r); write_new(qp, request)
            result = standard.prepare(qp, d / 'run', release_path=rp)
            plan = standard.checked_plan(verify(result['plan']))
            self.assertEqual(plan['backend'], 'fast_pqq')
            self.assertEqual(plan['source_domain_status']['1H4I'], 'unvalidated_input_domain')


if __name__ == '__main__': unittest.main()
