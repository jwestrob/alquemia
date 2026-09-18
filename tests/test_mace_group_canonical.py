"""Real archived parser/calibration tests and pinned full-protein source checks."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json
from mace_group_canonical import prepared, validate, calibrate, classify, parent_results
from mace_group_transfer import task
from mace_hybrid import EV_TO_KCAL

WORK = ROOT / 'workspaces/mace_group_canonical_20260918'
PREP = WORK / 'prepared_v1/preparation.json'
MANIFEST = WORK / 'model_v2/manifest.json'


@unittest.skipUnless(PREP.exists(), 'pinned real full-protein preparations unavailable')
class CanonicalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p, cls.cfg, cls.data, cls.cases = prepared(PREP)

    def test_all_25_real_pairs_preserve_geometry_groups_and_charge(self):
        self.assertEqual(len(self.cases), 25)
        self.assertEqual(sum(c['case']['expected_class'] == 'La' for c in self.cases.values()), 11)
        for name, row in self.cases.items():
            ca, a = task(name, 'Ca', 'primary', row)
            la, b = task(name, 'La', 'primary', row)
            np.testing.assert_array_equal(np.array(a, dtype=object)[:, 1:], np.array(b, dtype=object)[:, 1:])
            self.assertEqual(ca['atom_group_indices'], la['atom_group_indices'])
            delta = np.array(la['group_charges_e']) - ca['group_charges_e']
            self.assertEqual(np.flatnonzero(delta).tolist(), [row['primary']['selected_group_index']])
            self.assertEqual(delta.sum(), 1)
            self.assertEqual(la['charge'] - ca['charge'], 1)
            self.assertEqual(sum(ca['group_charges_e']), ca['charge'])
            self.assertEqual(ca['state']['all_electron_count'] % 2, 0)
            self.assertEqual(la['state']['all_electron_count'] % 2, 0)

    def test_calibration_algebra_on_archived_baseline_values_only(self):
        # Parser/algebra fixture, NOT new MACE evidence or a MACE calibration.
        scores = {r['case_id']: {'evaluation_role': r['evaluation_role'], 'expected_class': r['expected_class'],
                                 'R_model_kcal': r['baseline']['published_R_kcal_mol']} for r in self.data['rows']}
        c = calibrate(scores, True)
        self.assertEqual(c['status'], 'available_for_research')
        self.assertAlmostEqual(c['observed_class_gap_model_kcal'], 8.602932002686895, places=7)
        for s in scores.values():
            self.assertEqual(classify(s['R_model_kcal'], c), s['expected_class'])
        # Withholding a real case must disable the entire calibration.
        missing = copy.deepcopy(scores); missing[next(iter(self.cfg['cases']))]['R_model_kcal'] = None
        self.assertEqual(calibrate(missing, True)['status'], 'unavailable_incomplete')
        self.assertIsNone(classify(scores['1H4I']['R_model_kcal'], calibrate(missing, True)))
        self.assertEqual(calibrate(scores, False)['status'], 'unavailable_numerical_failure')

    def test_actual_failed_omol_panel_does_not_acquire_a_threshold(self):
        # Different archived model: exercise the refusal rule using real overlapping scores.
        old = read_json(ROOT / 'diagnostics/mace_omol_20260917/INTACT_CANONICAL_RESULT.json')
        scores = {r['case_id']: {'evaluation_role': r['evaluation_role'], 'expected_class': r['expected_class'],
                                'R_model_kcal': r['R_coord_kcal_mol']} for r in old['scores']}
        c = calibrate(scores, True)
        self.assertEqual(c['status'], 'unavailable_class_overlap')
        self.assertAlmostEqual(c['observed_class_gap_model_kcal'], old['calibration_gap_kcal_mol'], places=9)
        self.assertIsNone(c['Ca_max_R_model_kcal']); self.assertIsNone(c['La_min_R_model_kcal'])
        self.assertIsNone(classify(scores['1H4I']['R_model_kcal'], c))

    def test_reused_crystal_receipts_and_sign_match_actual_previous_result(self):
        parent, report, rows = parent_results(self.cfg)
        self.assertEqual(len(rows), 4)
        for name in ('1H4I', '4MAE'):
            ca = rows[name + '_Ca_primary']['energy_eV']; la = rows[name + '_La_primary']['energy_eV']
            r = (ca - la) * EV_TO_KCAL
            self.assertEqual(r, report['scores']['PQQ_' + name]['R_model_kcal'])
        self.assertFalse(report['representation_checks_pass'])

    @unittest.skipUnless(MANIFEST.exists(), 'actual finite calibration manifest unavailable')
    def test_manifest_and_explicit_corrupted_real_cache_and_model_rejected(self):
        self.assertEqual(validate(MANIFEST)['tasks'], 50)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'corrupted_real_manifest.json'
            m = read_json(MANIFEST); m['tasks'][0]['cache_key'] = 'corrupted_fixture'
            path.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):
                validate(path)
            m = read_json(MANIFEST); m['model']['chemical_group_policy'] += '_corrupted_fixture'
            path.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):
                validate(path)


@unittest.skipUnless((WORK / 'partial_report_v1/result.json').exists(), 'actual partial execution report unavailable')
class ActualPartialTests(unittest.TestCase):
    def test_actual_incomplete_execution_preserves_missing_scores_and_no_bands(self):
        d = read_json(WORK / 'partial_report_v1/result.json')
        self.assertEqual(d['status'], 'incomplete')
        self.assertEqual(d['calibration']['status'], 'unavailable_incomplete')
        self.assertEqual(len(d['scores']), 28)
        self.assertEqual(len(d['calibration_pairwise']), 154)
        self.assertFalse(d['representation_checks_pass'])
        self.assertIsNone(d['scores']['1KB0']['R_model_kcal'])
        self.assertIn('unsupported peptide connection', d['scores']['1KB0']['status'])
        endpoints = {**d['collection']['rows'], **d['collection']['reused_rows']}
        missing = 0
        for name, score in d['scores'].items():
            ca, la = (endpoints.get(name + '_' + metal + '_primary') for metal in ('Ca', 'La'))
            if ca and la and ca['status'] == la['status'] == 'computed':
                self.assertEqual(score['R_model_kcal'], (ca['energy_eV'] - la['energy_eV']) * EV_TO_KCAL)
            else:
                missing += 1
                self.assertIsNone(score['R_model_kcal'])
            self.assertIsNone(score['research_class']); self.assertIsNone(score['production_class'])
        self.assertGreater(missing, 1)


@unittest.skipUnless((WORK / 'report_v1/result.json').exists(), 'completed real canonical report unavailable')
class CompletedCanonicalTests(unittest.TestCase):
    def test_complete_actual_algebra_and_failed_calibration_are_preserved(self):
        d = read_json(WORK / 'report_v1/result.json')
        self.assertEqual(d['status'], 'complete'); self.assertTrue(d['numerical_checks_pass'])
        self.assertEqual(len(d['checks']), 54); self.assertTrue(all(c['pass'] for c in d['checks']))
        self.assertEqual(d['calibration_scored_count'], 25)
        self.assertEqual(d['pairwise_scored_count'], 154)
        self.assertEqual(d['retrospective_transfer_scored'], 2)
        self.assertEqual(d['retrospective_transfer_denominator'], 3)
        self.assertEqual(d['retrospective_transfer_classified'], 0)
        self.assertEqual(d['calibration']['status'], 'unavailable_class_overlap')
        self.assertLess(d['calibration']['observed_class_gap_model_kcal'], 0)
        self.assertFalse(d['representation_checks_pass'])
        rows = {**d['collection']['rows'], **d['collection']['reused_rows']}
        for name, score in d['scores'].items():
            self.assertIsNone(score['research_class']); self.assertIsNone(score['production_class'])
            if name == '1KB0':
                self.assertIsNone(score['R_model_kcal']); continue
            expected = (rows[name + '_Ca_primary']['energy_eV'] - rows[name + '_La_primary']['energy_eV']) * EV_TO_KCAL
            self.assertEqual(score['R_model_kcal'], expected)
        pairs = d['calibration_pairwise']
        for pair in pairs:
            margin = d['scores'][pair['La']]['R_model_kcal'] - d['scores'][pair['Ca']]['R_model_kcal']
            self.assertEqual(pair['margin_model_kcal'], margin)
            self.assertEqual(pair['pass'], margin > .02)
        self.assertEqual(d['pairwise_pass_count'], sum(p['pass'] for p in pairs))

    def test_current_report_replays_actual_complete_scientific_values(self):
        from mace_group_canonical import report
        original = read_json(WORK / 'report_v1/result.json')
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'replay'; report(MANIFEST, out)
            replay = read_json(out / 'result.json')
        for key in ('scores', 'checks', 'calibration', 'calibration_pairwise', 'grouping_checks',
                    'numerical_checks_pass', 'representation_checks_pass', 'pairwise_pass_count'):
            self.assertEqual(replay[key], original[key])


if __name__ == '__main__':
    unittest.main()
