"""Admission/recovery tests using pinned real response-trace manifests."""
from pathlib import Path
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json
from mace_response_trace import validate, core_gate
from mace_hybrid import collect

MANIFEST = ROOT / 'workspaces/mace_response_trace_20260916/medium_v1/manifest.json'


@unittest.skipUnless(MANIFEST.exists(), 'requires real prepared trace manifest')
class ResponseTraceTests(unittest.TestCase):
    def test_trace_admission_and_physical_model_change_rejected(self):
        m = read_json(MANIFEST)
        validate(m)
        m['model']['external_field'] = [1., 0., 0.]
        with self.assertRaisesRegex(InvalidArtifact, 'physical model'):
            validate(m)

    def test_corrupted_real_charge_and_tolerances_rejected(self):
        m = read_json(MANIFEST)
        m['tasks'][0]['charge'] += 1
        with self.assertRaisesRegex(InvalidArtifact, 'physical input'):
            validate(m)
        m = read_json(MANIFEST)
        m['trace_tolerances']['density'] = 1.
        with self.assertRaisesRegex(InvalidArtifact, 'tolerances'):
            validate(m)

    def test_missing_trace_execution_cannot_pass_gate_or_supply_score(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / 'manifest.json'
            p.write_bytes(MANIFEST.read_bytes())
            c = collect(p)
            self.assertEqual(c['status'], 'incomplete')
            self.assertEqual(len(c['rows']), 6)
            self.assertEqual(core_gate(p)['status'], 'failed')
            self.assertIsNone(c['direct'])
            self.assertIsNone(c['S_kcal_mol'])


if __name__ == '__main__':
    unittest.main()
