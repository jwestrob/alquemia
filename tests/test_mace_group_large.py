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
from mace_group_large import sources, tasks, validate, native_gate

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


if __name__ == '__main__':
    unittest.main()
