"""Real archived endpoint algebra and declared fold-pool accounting; no fake science."""
from __future__ import annotations
import copy
from pathlib import Path
import sys
import statistics
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify
import accommodation_folds_compare as comparison


class FoldComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = read_json(ROOT / 'diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json')
        cls.reference = read_json(ROOT / 'workspaces/compact_solvation_20260920/full_v1/comparison_v1.json')
        cls.inventory_path = ROOT / 'diagnostics/compact_solvation_20260920/INVENTORY.json'
        cls.inventory = read_json(cls.inventory_path)
        cls.collection = read_json(ROOT / 'workspaces/compact_solvation_20260920/full_v1/collection_1203165.json')

    def test_actual_source_scope_and_primary_membership(self):
        groups = comparison.source_groups(self.source)
        self.assertEqual(len(groups), 25)
        self.assertEqual(sum(c['canonical_coordinate_match'] for c in self.source['cases']), 25)
        self.assertEqual(sum(c['primary_evaluation_pool'] for c in self.source['cases']), 225)
        for members in groups.values():
            self.assertEqual(sum(c['source_conditioning_metal'] == 'La' and not c['canonical_coordinate_match'] for c in members), 4)
            self.assertEqual(sum(c['source_conditioning_metal'] == 'Ca' for c in members), 5)

    def test_corrupted_real_replay_assignment_rejected(self):
        source = copy.deepcopy(self.source)
        source['cases'][0]['canonical_coordinate_match'] = True
        with self.assertRaises(InvalidArtifact): comparison.source_groups(source)

    def test_old_bands_are_exact_and_edges_not_padded(self):
        bands = comparison.bands(self.reference)
        self.assertEqual(bands['core']['composite']['Ca_max'], -405480.04029496165)
        self.assertEqual(bands['context']['composite']['La_min'], -405459.1113819997)
        for b in bands.values():
            for band in b.values():
                self.assertEqual(comparison.decision(band['Ca_max'], band), 'Ca-supported')
                self.assertEqual(comparison.decision(band['La_min'], band), 'La-supported')
                self.assertEqual(comparison.decision((band['Ca_max'] + band['La_min']) / 2, band), 'inconclusive')
                self.assertEqual(comparison.decision(None, band), 'unavailable')

    def test_actual_reused_receipts_keep_original_task_identity(self):
        row = next(r for r in self.collection['rows'] if (r['case_id'], r['representation'], r['metal']) == ('1H4I', 'context', 'Ca'))
        source = next(c for c in self.inventory['cases'] if c['case_id'] == '1H4I')['representations']['context']['endpoints']['Ca']
        reuse = comparison.checked_scoped_transfer(row, source, record(self.inventory_path), self.collection['manifest'])
        self.assertEqual(set(reuse), {'vacuum', 'alpb'})
        for medium in reuse:
            self.assertNotEqual(reuse[medium]['original_manifest'], self.collection['manifest'])
            self.assertEqual(reuse[medium]['original_task_id'], row['endpoints'][medium]['task_id'])

    def test_actual_fresh_receipt_and_corrupted_transfer_sign(self):
        row = next(r for r in self.collection['rows'] if (r['case_id'], r['representation'], r['metal']) == ('4MAE', 'core', 'La'))
        source = next(c for c in self.inventory['cases'] if c['case_id'] == '4MAE')['representations']['core']['endpoints']['La']
        reuse = comparison.checked_scoped_transfer(row, source, record(self.inventory_path), self.collection['manifest'])
        self.assertEqual(reuse, {})
        corrupted = copy.deepcopy(row); corrupted['delta_solv_hartree'] *= -1
        with self.assertRaises(InvalidArtifact):
            comparison.checked_scoped_transfer(corrupted, source, record(self.inventory_path), self.collection['manifest'])

    def test_real_canonical_values_do_not_fill_missing_new_folds(self):
        frozen = comparison.bands(self.reference); rows = []
        old = {(r['case_id'], r['representation']): r for r in self.reference['rows']}
        for c in self.source['cases']:
            for rep in comparison.REPS:
                row = {k: c[k] for k in ('case_id', 'root_case_id', 'source_conditioning_metal', 'canonical_coordinate_match',
                        'primary_evaluation_pool', 'biological_group', 'expected_class')}
                row.update(representation=rep, state_audit={'status': 'unavailable'}, **{k: None for k in comparison.COMPONENTS})
                if c['canonical_coordinate_match']:
                    # Exact established canonical coordinate mapping; actual archived energy only.
                    actual = old[(c['root_case_id'], rep)]
                    row.update({k: actual[k] for k in comparison.COMPONENTS})
                for method, field in [('native', comparison.COMPONENTS[0]), ('composite', comparison.COMPONENTS[-1])]:
                    row[method + '_decision'] = comparison.decision(row[field], frozen[rep][method])
                rows.append(row)
        summaries, counts = comparison.summarize(self.source, rows, frozen)
        for rep in comparison.REPS:
            totals = counts[rep]
            self.assertEqual(totals['single_sources']['canonical25_replay']['composite']['correct'], 25)
            self.assertEqual(totals['single_sources']['primary225']['composite']['unavailable'], 225)
            for group in totals['ensemble_descriptors'].values():
                self.assertEqual(group['composite']['unavailable'], 25)
        for protein in summaries:
            self.assertIsNone(protein['equal_mean_of_arm_medians']['composite_R_model_kcal_mol'])
            self.assertIsNone(protein['Ca5_median_minus_La4_median']['composite_R_model_kcal_mol'])

    def test_actual_new_preparation_retains_state_and_context_membership(self):
        directory = ROOT / 'workspaces/accommodation_goal_20260920/folds_v1/source_preparation/cases'
        pin = next(directory.glob('*/preparation_result.json'))
        prepared = read_json(pin)
        audit = comparison.state_audit(prepared)
        self.assertEqual(audit['status'], 'available')
        self.assertTrue(audit['core_source_coordinates_unchanged'])
        self.assertTrue(audit['water_inventory_unchanged'])
        self.assertGreater(audit['context_atoms'], 0)
        self.assertEqual(audit['core_state']['charges']['La'] - audit['core_state']['charges']['Ca'], 1)

    def test_actual_warm_native_complete_pool_and_missing_member(self):
        actual = read_json(ROOT / 'workspaces/accommodation_goal_20260920/folds_v1/native_comparison_v1.json')
        frozen = comparison.bands(self.reference)
        group = comparison.source_groups(self.source)['a0a3f2yly8-pqq-la_model']
        members = [m for m in group if m['source_conditioning_metal'] == 'La' and not m['canonical_coordinate_match']]
        index = {r['case_id']: r for r in actual['rows'] if r['representation'] == 'core'}
        result = comparison.pool_summary(members, index, frozen['core'])
        values = [index[m['case_id']]['native_R_model_kcal_mol'] for m in members]
        self.assertEqual(result['components']['native_R_model_kcal_mol']['median'], statistics.median(values))
        self.assertEqual(result['components']['native_R_model_kcal_mol']['range'], max(values) - min(values))
        self.assertEqual(result['components']['composite_R_model_kcal_mol']['status'], 'unavailable')
        corrupted = copy.deepcopy(index); corrupted[members[0]['case_id']]['native_R_model_kcal_mol'] = None
        missing = comparison.pool_summary(members, corrupted, frozen['core'])
        self.assertIsNone(missing['components']['native_R_model_kcal_mol']['median'])
        self.assertEqual(missing['native_decision'], 'unavailable')

    def test_actual_failed_vacuum_retains_valid_ALPB_and_null_transfer(self):
        directory = ROOT / 'workspaces/accommodation_goal_20260920/folds_v1'
        col = read_json(directory / 'solvent_shards_v2/shard_3/collection_after_failure_1203800.json')
        inv_path = directory / 'INVENTORY.json'; inv = read_json(inv_path)
        failed = next(r for r in col['rows'] if r['status'] != 'complete')
        self.assertEqual(failed['case_id'], 'a8r3s4-pqq-la_model__conditioned_Ca__seed-1_sample-3')
        self.assertEqual(failed['endpoints']['alpb']['status'], 'complete')
        self.assertIsNone(failed['vacuum_hartree']); self.assertIsNone(failed['delta_solv_hartree'])
        case = next(c for c in inv['cases'] if c['case_id'] == failed['case_id'])
        endpoint = case['representations'][failed['representation']]['endpoints'][failed['metal']]
        self.assertEqual(comparison.checked_scoped_transfer(failed, endpoint, record(inv_path), col['manifest']), {})

    def test_complete_actual_result_keeps_denominators_and_La4_gain_distinct(self):
        actual = read_json(ROOT / 'workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json')
        self.assertEqual(len(actual['rows']), 500)
        self.assertEqual(sum(r['composite_status'] == 'available' for r in actual['rows']), 464)
        self.assertEqual(sum(len(r['endpoints']) for r in actual['explicit_original_receipt_reuse']), 196)
        for rep in comparison.REPS:
            self.assertEqual(actual['counts'][rep]['single_sources']['canonical25_replay']['composite']['correct'], 25)
        context = actual['counts']['context']
        self.assertEqual(context['single_sources']['primary225']['composite']['denominator'], 225)
        self.assertEqual(context['single_sources']['primary225']['composite']['wrong'], 2)
        self.assertEqual(context['ensemble_descriptors']['La4']['native']['wrong'], 1)
        self.assertEqual(context['ensemble_descriptors']['La4']['composite']['correct'], 23)
        self.assertEqual(context['ensemble_descriptors']['balanced']['composite']['unavailable'], 5)
        q9 = next(p for p in actual['protein_summaries'] if p['representation'] == 'context'
                  and p['root_case_id'] == 'q9z4j7-pqq-la_model')
        self.assertEqual((q9['pools']['La4']['native_decision'], q9['pools']['La4']['composite_decision']),
                         ('La-supported', 'Ca-supported'))
        for row in actual['rows']:
            if row['preparation_status'] == 'prepared' and row['composite_status'] != 'available':
                self.assertIsNotNone(row['GFN2_ALPB_R_kcal_mol'])
                self.assertIsNone(row['GFN2_vacuum_R_kcal_mol'])
                self.assertIsNone(row['solvation_delta_R_kcal_mol'])


if __name__ == '__main__': unittest.main()
