"""Strict three-fold checks on real sources and archived scalar values."""
import copy
from pathlib import Path
import statistics
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, write_new
import pqq_ensemble as en


class Ensemble(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.request = read_json(ROOT / 'diagnostics/pqq_ensemble_20260920/REQUEST.json')
        cls.group = en.declared_groups(cls.request)[0]
        comparison = read_json(ROOT / 'workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json')
        cls.rows = {r['case_id']: {**r, 'source_domain_status': 'unvalidated_input_domain'}
                    for r in comparison['rows'] if r['representation'] == 'context' and r['case_id'] in cls.group['members']}
        ref = read_json(ROOT / 'workspaces/compact_solvation_20260920/full_v1/comparison_v1.json')
        cls.bands = ref['calibration']['context']['bands']

    def test_real_sources_match_actual_sequence_roles_and_La(self):
        self.assertEqual(self.group['source_identity']['sequence_length'], 623)
        self.assertEqual(len(set(c['source_structure']['sha256'] for c in self.request['cases'])), 3)
        self.assertEqual([c['case_id'].rsplit('-', 1)[-1] for c in self.request['cases']], ['1', '2', '3'])

    def test_corrupted_membership_conditioning_and_roles_fail(self):
        for fault in ('duplicate', 'conditioning', 'role', 'protein'):
            bad = copy.deepcopy(self.request)
            if fault == 'duplicate': bad['ensemble']['groups'][0]['members'][1] = bad['cases'][0]['case_id']
            if fault == 'conditioning': bad['cases'][0]['source_conditioning_metal'] = 'Ca'
            if fault == 'role': bad['cases'][0]['roles']['anchor_asparagine']['resname'] = 'ALA'
            if fault == 'protein': bad['cases'][0]['root_case_id'] = 'different_protein'
            with self.assertRaises(InvalidArtifact): en.declared_groups(bad)

    def test_corrupted_real_source_sequence_fails(self):
        import gemmi
        bad = copy.deepcopy(self.request)
        st = gemmi.read_structure(bad['cases'][0]['source_structure']['path'])
        # Explicitly corrupted copy of the real first residue; no scientific run.
        first = st[0][0][0]; first.name = 'GLY' if first.name != 'GLY' else 'ALA'
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'corrupted_sequence.cif'
            st.make_mmcif_document().write_file(str(path))
            bad['cases'][0]['source_structure'] = record(path)
            with self.assertRaisesRegex(InvalidArtifact, 'sequence'): en.declared_groups(bad)

    def test_actual_raw_values_median_and_spread(self):
        result = en.aggregate(self.group, self.rows, self.bands)
        values = [self.rows[c]['composite_R_model_kcal_mol'] for c in self.group['members']]
        self.assertEqual(result['median_R_model_kcal_mol'], statistics.median(values))
        self.assertEqual(result['components']['composite_R_model_kcal_mol']['spread'], max(values) - min(values))
        self.assertEqual(result['developmental_frozen_band_transfer'], 'Ca-supported')
        self.assertEqual(result['source_domain_status'], ['unvalidated_input_domain'])

    def test_extra_real_protein_chain_in_corrupted_copy_changes_assembly(self):
        import gemmi
        bad = copy.deepcopy(self.request)
        st = gemmi.read_structure(bad['cases'][0]['source_structure']['path'])
        # Explicitly corrupted fixture: duplicate a real chain under a new ID.
        extra = st[0][0].clone(); extra.name = 'Z'; st[0].add_chain(extra)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'corrupted_extra_chain.cif'; st.make_mmcif_document().write_file(str(path))
            bad['cases'][0]['source_structure'] = record(path)
            with self.assertRaisesRegex(InvalidArtifact, 'sequence'): en.declared_groups(bad)

    def test_missing_real_score_makes_entire_descriptor_unavailable(self):
        bad = copy.deepcopy(self.rows)
        bad[self.group['members'][0]]['composite_R_model_kcal_mol'] = None  # removed actual value
        result = en.aggregate(self.group, bad, self.bands)
        self.assertEqual(result['status'], 'unavailable')
        self.assertIsNone(result['median_R_model_kcal_mol'])
        self.assertIsNone(result['components']['composite_R_model_kcal_mol']['spread'])
        self.assertEqual(result['available_members'], 2)
        self.assertEqual(result['required_members'], 3)

    def test_existing_standard_result_cannot_gain_posthoc_members(self):
        actual = ROOT / 'workspaces/pqq_fast_release_20260920/standard_1H4I_v1'
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(InvalidArtifact, 'original source request'):
                en.collect(actual / 'plan.json', Path(d) / 'result.json', actual / 'result_1203729.json')

    def test_fresh_literal_standard_execution_and_replay(self):
        run = ROOT / 'workspaces/pqq_ensemble_20260920/pilot_v1'
        plan = read_json(run / 'plan.json'); execution = read_json(run / 'execution.json')
        self.assertEqual((plan['new_MACE_calls'], plan['new_GFN2_calls'], plan['new_DFT_calls']), (6, 12, 0))
        self.assertTrue(execution['source_preparation_included'])
        self.assertFalse(execution['folding_included'])
        actual = read_json(run / 'ensemble_1204055.json')
        self.assertEqual(actual['status'], 'complete')
        self.assertFalse(actual['default_changed']); self.assertFalse(actual['thermal_ensemble'])
        self.assertFalse(actual['threshold_refitted']); self.assertIsNone(actual['probabilities'])
        group = actual['groups'][0]
        self.assertEqual(group['developmental_frozen_band_transfer'], 'Ca-supported')
        self.assertGreater(group['components']['composite_R_model_kcal_mol']['spread'], 20)
        for row in group['members']:
            self.assertEqual(row['source_domain_status'], 'unvalidated_input_domain')
            self.assertLess(abs(row['composite_R_model_kcal_mol'] - self.rows[row['case_id']]['composite_R_model_kcal_mol']), .02)
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'replay.json'
            en.collect(run / 'plan.json', out, run / 'result_1204055.json')
            self.assertEqual(read_json(out)['groups'], actual['groups'])

    def test_unexecuted_or_missing_members_do_not_gain_a_median(self):
        run = ROOT / 'workspaces/pqq_ensemble_20260920/pilot_v1'
        actual = read_json(run / 'pre_execution_unavailable.json')
        self.assertEqual(actual['groups'][0]['available_members'], 0)
        self.assertIsNone(actual['groups'][0]['median_R_model_kcal_mol'])
        # A removed row from the actual collection is explicitly partial input.
        bad = read_json(run / 'result_1204055.json'); bad['rows'].pop()
        with tempfile.TemporaryDirectory() as d:
            d = Path(d); cp = d / 'missing_row.json'; op = d / 'summary.json'; write_new(cp, bad)
            en.collect(run / 'plan.json', op, cp)
            group = read_json(op)['groups'][0]
            self.assertEqual(group['available_members'], 2)
            self.assertIsNone(group['median_R_model_kcal_mol'])

    def test_corrupted_real_result_bands_score_and_domain_rejected(self):
        run = ROOT / 'workspaces/pqq_ensemble_20260920/pilot_v1'
        for field in ('bands', 'score', 'domain'):
            bad = read_json(run / 'result_1204055.json')
            if field == 'bands': bad['published_context_bands']['Ca_supported_max_R_model_kcal_mol'] += 1
            if field == 'score': bad['rows'][0]['composite_R_model_kcal_mol'] += 1
            if field == 'domain': bad['rows'][0]['source_domain_status'] = 'validated_everywhere'
            with tempfile.TemporaryDirectory() as d:
                d = Path(d); cp = d / 'corrupted.json'; write_new(cp, bad)
                with self.assertRaises(InvalidArtifact): en.collect(run / 'plan.json', d / 'output.json', cp)


if __name__ == '__main__': unittest.main()
