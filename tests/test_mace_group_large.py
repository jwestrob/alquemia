"""Pinned real-state and execution-gate checks for the larger grouped model."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, xyz
from mace_group_large import sources, tasks, validate, native_gate, report
from mace_hybrid import EV_TO_KCAL

WORK = ROOT / 'workspaces/mace_group_large_20260918'
CONFIG = WORK / 'config.json'
MANIFEST = WORK / 'model_v1/manifest.json'


@unittest.skipUnless(MANIFEST.exists(), 'real finite large-checkpoint manifest unavailable')
class LargeGroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg, cls.small, cls.large, cls.previous, cls.native = sources(CONFIG)
        cls.m = read_json(MANIFEST)

    def test_all_55_scientific_inputs_are_exactly_the_existing_medium_panel(self):
        self.assertEqual(len(self.m['tasks']), 57)
        for new, old in zip(self.m['tasks'][2:], self.small['tasks']):
            self.assertEqual({k: v for k, v in new.items() if k != 'cache_key'},
                             {k: v for k, v in old.items() if k != 'cache_key'})
            self.assertNotEqual(new['cache_key'], old['cache_key'])
        self.assertEqual(self.m['model']['checkpoint'], self.large['model']['checkpoint'])
        self.assertEqual(self.m['model']['pair_kernel'], self.large['model']['pair_kernel'])
        for key in ('chemical_group_policy', 'charge_group_adapter', 'native_polar_source'):
            self.assertEqual(self.m['model'][key], self.small['model'][key])

    def test_native_controls_preserve_actual_archived_large_coordinate_order(self):
        for t in self.m['tasks'][:2]:
            old = next(r for r in self.large['tasks'] if r['task_id'] == 'GGR_1GLG_' + t['metal'] + '_primary')
            self.assertEqual(t['xyz'], old['xyz']); self.assertEqual(t['charge'], old['charge'])
            self.assertEqual(t['preparation'], old['preparation'])
            self.assertEqual(t['atom_group_indices'], [0] * len(xyz(verify(t['xyz']))))
            self.assertEqual(t['group_charges_e'], [t['charge']])

    def test_missing_native_execution_blocks_the_gate(self):
        # Real manifest copied to an empty execution location; no fabricated outputs.
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'real_manifest_without_execution.json'
            path.write_text(json.dumps(self.m))
            gate = native_gate(path)
            self.assertEqual(gate['status'], 'unavailable_or_failed')
            self.assertEqual(len(gate['checks']), 6)
            self.assertTrue(all(c['error'] is None and not c['pass'] for c in gate['checks']))

    def test_actual_manifest_and_explicit_corrupted_checkpoint_rejected(self):
        self.assertEqual(validate(MANIFEST)['tasks'], 57)
        broken = copy.deepcopy(self.m); broken['model']['checkpoint'] = self.small['model']['checkpoint']
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'corrupted_real_manifest.json'; path.write_text(json.dumps(broken))
            with self.assertRaises(InvalidArtifact):
                validate(path)


@unittest.skipUnless((WORK / 'report_v1/checkpoint_comparison.json').exists(),
                     'completed real large-checkpoint result unavailable')
class CompletedLargeTests(unittest.TestCase):
    def test_actual_scores_native_checks_and_checkpoint_comparisons(self):
        d = read_json(WORK / 'report_v1/result.json')
        paired = read_json(WORK / 'report_v1/checkpoint_comparison.json')
        old = read_json(verify(paired['medium_result']))
        self.assertEqual(d['status'], 'complete')
        self.assertEqual(len(d['collection']['rows']), 57)
        self.assertEqual(len(d['scores']), 34)
        self.assertEqual(len(d['checks']), 69)
        self.assertEqual(len(d['comparisons']), 22)
        self.assertEqual(d['collection']['native_identity']['status'], 'pass')
        self.assertEqual(d['numerical_checks_pass'], all(c['pass'] for c in d['checks']))
        self.assertEqual(d['representation_checks_pass'], all(c['pass'] for c in d['grouping_checks']))
        rows = d['collection']['rows']
        for s in d['scores'].values():
            self.assertEqual(s['R_model_kcal'],
                             (rows[s['Ca_task']]['energy_eV'] - rows[s['La_task']]['energy_eV']) * EV_TO_KCAL)
            self.assertIsNone(s['calibrated_class'])
        for new, previous, p in zip(d['comparisons'], old['comparisons'], paired['paired']):
            for key in ('higher_expected', 'lower_expected', 'stratum'):
                self.assertEqual(new[key], previous[key])
            self.assertEqual(p['large_minus_medium_margin_model_kcal'],
                             new['margin_model_kcal'] - previous['margin_model_kcal'])
            self.assertEqual(new['pass'], new['margin_model_kcal'] > .02)
        passed = sum(c['pass'] for c in d['comparisons'])
        self.assertEqual(d['raw_pass_count'], passed)
        self.assertEqual(paired['declared_directional_improvement'],
                         passed > 17 and all(c['pass'] for c in d['comparisons'][:7]))
        self.assertEqual(d['qualified_pass_count'], passed if d['numerical_checks_pass'] and d['representation_checks_pass'] else 0)
        self.assertEqual([r['site'] for r in d['aequorin_ordered']],
                         ['AEQ_1SL8_EF1', 'AEQ_1SL8_EF3', 'AEQ_1SL8_EF4'])
        self.assertTrue(all(r['site_label'] is None for r in d['aequorin_ordered']))

    def test_report_replays_actual_complete_scientific_values(self):
        original = read_json(WORK / 'report_v1/result.json')
        original_pair = read_json(WORK / 'report_v1/checkpoint_comparison.json')
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'replay'; report(MANIFEST, out)
            replay = read_json(out / 'result.json')
            paired = read_json(out / 'checkpoint_comparison.json')
        for key in ('scores', 'checks', 'grouping_checks', 'comparisons', 'raw_pass_count',
                    'qualified_pass_count', 'aequorin_ordered', 'Khoury_means'):
            self.assertEqual(replay[key], original[key])
        self.assertEqual(paired['paired'], original_pair['paired'])
        self.assertEqual(paired['declared_directional_improvement'], original_pair['declared_directional_improvement'])


if __name__ == '__main__':
    unittest.main()
