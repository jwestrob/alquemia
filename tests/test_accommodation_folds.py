"""Actual pinned source and cross-host replay checks; no invented molecules."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from accommodation_folds import adapt_ca, atoms, reconcile
from affordable_common import InvalidArtifact, read_json, verify
import gemmi


class FoldInputs(unittest.TestCase):
    def test_real_ca_adapter_preserves_every_atom(self):
        m = read_json(ROOT / 'diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json')
        source = next(c for c in m['cases'] if c['source_conditioning_metal'] == 'Ca')
        with tempfile.TemporaryDirectory() as d:
            result = adapt_ca(source, Path(d) / 'alias')
            original = gemmi.read_structure(str(verify(source['source_structure'])))
            actual = gemmi.read_structure(str(verify(result['source_structure'])))
            self.assertEqual(atoms(original), atoms(actual))
            self.assertEqual(actual[0]['B'][0].name, 'LA')
            self.assertEqual(actual[0]['B'][0][0].element.name, 'Ca')

    def test_corrupted_real_ca_selector_is_rejected(self):
        m = read_json(ROOT / 'diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json')
        source = copy.deepcopy(next(c for c in m['cases'] if c['source_conditioning_metal'] == 'Ca'))
        source['raw_source_metal']['atom'] = 'CORRUPTED_REAL_SELECTOR'
        with tempfile.TemporaryDirectory() as d, self.assertRaises(InvalidArtifact):
            adapt_ca(source, Path(d) / 'alias')

    def test_actual_last_bit_replays_preserve_failures_and_coordinates(self):
        source = ROOT / 'workspaces/accommodation_goal_20260920/folds_v1/source_preparation/preparation.json'
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'reconciled.json'; reconcile(source, path)
            old, new = read_json(source), read_json(path)
            self.assertEqual(len(old['cases']), 250)
            self.assertEqual(len(new['cases']), 250)
            self.assertEqual(new['numeric_reconciliation']['case_ids'].__len__(), 7)
            for a, b in zip(old['cases'], new['cases']):
                if a['status'] == 'unsupported': self.assertEqual(a, b)
                if 'representations' in a:
                    self.assertEqual(a['representations'], b['representations'])
                    self.assertEqual(a['core'], b['core'])
                    self.assertEqual(a['context_comparison'], b['context_comparison'])


if __name__ == '__main__': unittest.main()
