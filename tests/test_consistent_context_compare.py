"""Actual union-context component/reference checks; no molecular calls."""
from pathlib import Path
import copy
import json
import statistics
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json
from consistent_context_compare import collect, reference, transfer_compare
from compact_solvation_compare import mix_pair

BASE = ROOT/'workspaces/consistent_context_20260922'


class ActualUnionComponents(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.saved = read_json(BASE/'calibration28_v1/collection_final_v2.json')

    def test_actual_collection_complete_status_and_component_algebra(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'actual.json'
            collect(BASE/'calibration28_v1', path)
            current = read_json(path)
        self.assertEqual((current['complete'], current['denominator']), (28, 28))
        old = {r['case_id']: r for r in self.saved['rows']}
        for row in current['rows']:
            self.assertEqual(row['status'], 'complete')
            self.assertTrue(all(e['status'] == 'complete' for e in row['native_endpoints'].values()))
            self.assertEqual(row['composite_R_model_kcal_mol'], old[row['case_id']]['composite_R_model_kcal_mol'])
            actual = mix_pair({z: row['native_endpoints'][z]['native_MACE_energy_eV'] for z in ('Ca', 'La')},
                    {z: row['solvent_endpoints'][z]['vacuum']['energy_hartree'] for z in ('Ca', 'La')},
                    {z: row['solvent_endpoints'][z]['alpb']['energy_hartree'] for z in ('Ca', 'La')})
            for key, value in actual.items(): self.assertEqual(row[key], value)

    def test_frozen_reference_uses_exact_canonical_only_values(self):
        old = read_json(BASE/'calibration28_v1/REFERENCE_v2.json')
        with tempfile.TemporaryDirectory() as d:
            new = reference(BASE/'calibration28_v1/collection_final_v2.json', Path(d)/'reference.json')
        for key in ('bands', 'class_extrema', 'gap_model_kcal_mol', 'canonical_case_ids', 'status'):
            self.assertEqual(new[key], old[key])
        self.assertEqual(new['available'], 25)
        self.assertFalse(new['transfers_used'])
        self.assertEqual(new['ordered_crossclass_pairs'], new['crossclass_pairs'])

    def test_corrupted_real_missing_endpoint_cannot_make_reference(self):
        damaged = copy.deepcopy(self.saved)
        row = next(r for r in damaged['rows'] if r['role'] == 'calibration')
        row['status'] = 'unavailable'; row['composite_R_model_kcal_mol'] = None
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'explicitly_corrupted_actual_collection.json'; path.write_text(json.dumps(damaged))
            ref = reference(path, Path(d)/'reference.json')
        self.assertEqual(ref['status'], 'unavailable_incomplete_canonical')
        self.assertIsNone(ref['bands'])
        self.assertEqual(ref['available'], 24)


class ActualUnionTransfer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.output = Path(cls.temporary.name)/'comparison.json'
        transfer_compare(BASE/'pilot36_v1/collection_final_v1.json', BASE/'calibration28_v1/REFERENCE_v2.json',
              ROOT/'workspaces/nikasha_recovery_20260922/proposal_comparison.json',
              ROOT/'diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json', cls.output)
        cls.result = read_json(cls.output)

    @classmethod
    def tearDownClass(cls): cls.temporary.cleanup()

    def test_same_coverage_and_original_baseline_values(self):
        rows = self.result['rows']; self.assertEqual(len(rows), 36)
        old = {r['case_id']: r for r in read_json(ROOT/'workspaces/nikasha_recovery_20260922/proposal_comparison.json')['rows']}
        for r in rows:
            for method, actual in old[r['case_id']]['methods'].items(): self.assertEqual(r['methods'][method], actual)
        counts = self.result['counts']['single_sources']['all']['context_union']
        self.assertEqual(counts, {'denominator': 36, 'correct': 19, 'wrong': 1, 'inconclusive': 5, 'unavailable': 11})
        self.assertEqual(sum(v for k, v in counts.items() if k != 'denominator'), 36)

    def test_strict_medians_missing_members_and_balanced_weights(self):
        idx = {r['case_id']: r['methods']['context_union']['R'] for r in self.result['rows']}
        self.assertEqual((len(self.result['pools']), len(self.result['triples'])), (12, 16))
        for row in self.result['pools']+self.result['triples']:
            v = row['methods']['context_union']; missing = [c for c in row['members'] if idx[c] is None]
            self.assertEqual(set(v['missing_members']), set(missing))
            if missing: self.assertIsNone(v['R'])
            elif row.get('pool') != 'balanced':
                self.assertEqual(v['R'], statistics.median(idx[c] for c in row['members']))
            else:
                arms = [p['methods']['context_union']['R'] for p in self.result['pools']
                        if p['root_case_id'] == row['root_case_id'] and p['pool'] in ('La4', 'Ca5')]
                self.assertEqual(v['R'], sum(arms)/2)

    def test_corrupted_real_reference_pin_rejected(self):
        ref = read_json(BASE/'calibration28_v1/REFERENCE_v2.json')
        ref['bands']['Ca_max'] += 1.0  # Explicit corruption of a real frozen reference.
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'corrupted_reference.json'; path.write_text(json.dumps(ref))
            with self.assertRaisesRegex(InvalidArtifact, 'another reference'):
                transfer_compare(BASE/'pilot36_v1/collection_final_v1.json', path,
                    ROOT/'workspaces/nikasha_recovery_20260922/proposal_comparison.json',
                    ROOT/'diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json', Path(d)/'out.json')

    def test_actual_full_manifest_retains_scope_and_reuses_pilot_receipts(self):
        ready = read_json(BASE/'primary225_v1/READY.json')
        self.assertEqual((ready['stage_denominator'], ready['prepared_cases']), (225, 208))
        self.assertEqual((ready['MACE_calls'], ready['MACE_reused'], ready['GFN2_calls'], ready['GFN2_reused']),
                         (188, 228, 378, 454))
        native = read_json(BASE/'primary225_v1/mace/manifest.json')
        solvent = read_json(BASE/'primary225_v1/solvent/manifest.json')
        prior = read_json(BASE/'pilot36_v1/collection_final_v1.json')
        inv = read_json(BASE/'primary225_v1/SCORE_INPUTS.json')
        self.assertEqual(len(inv['cases'])+len(inv['failures']), 225)
        for row in prior['rows']:
            if row['status'] != 'complete': continue
            for metal in ('Ca', 'La'):
                tid = row['case_id']+'__context__'+metal
                self.assertEqual(native['reused'][tid]['receipt'], row['native_endpoints'][metal]['native_MACE_receipt'])
                for medium in ('vacuum', 'alpb'):
                    pin = solvent['reused'][tid+'__'+medium+'__native']; old = row['solvent_endpoints'][metal][medium]
                    for key in ('receipt', 'manifest', 'output', 'task_id', 'energy_hartree'):
                        self.assertEqual(pin[key], old[key])


if __name__ == '__main__': unittest.main()
