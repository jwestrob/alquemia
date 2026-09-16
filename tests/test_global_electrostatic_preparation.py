"""Preparation invariants on the four real, archived 1H4I endpoints."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from global_electrostatic import prepare_partition, input_state, dry_run, validate_vacuum_output

PILOT = ROOT / 'workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json'
ENV = ROOT / 'workspaces/affordable_challenger_20260915/environment'
AGREEMENT = ROOT / 'diagnostics/global_electrostatic_20260916/AGREEMENT.md'


class ActualPartitionInputs(unittest.TestCase):
    def test_actual_vacuum_run_ignores_smd_author_credit(self):
        p = ROOT / 'workspaces/global_electrostatic_20260916/partition_tasks_v1/manifest.json'
        if not p.exists():
            self.skipTest('approved vacuum integration has not run in this checkout')
        t = next(t for t in read_json(p)['tasks'] if t['task_id'] == '1h4i_qm33_Ca')
        op = Path(t['output_path'])
        if not op.exists() or 'ORCA TERMINATED NORMALLY' not in op.read_text():
            self.skipTest('approved real vacuum endpoint is incomplete')
        self.assertIn('smd solvation model', op.read_text())
        self.assertTrue(validate_vacuum_output(t, op)['vacuum'])

    def test_exact_coordinates_charges_and_explicit_solvent_change(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'prepared'
            m = prepare_partition(PILOT, ENV, out, AGREEMENT)
            self.assertEqual(len(m['tasks']), 4)
            self.assertEqual(dry_run(out / 'manifest.json')['status'], 'dry_run_pass')
            for t in m['tasks']:
                self.assertEqual(verify(t['xyz']).read_bytes(), verify(t['source_xyz']).read_bytes())
                self.assertEqual(input_state(verify(t['input']))['charge'], t['charge'])
                self.assertNotIn('CPCM', verify(t['input']).read_text())
                with self.assertRaises(InvalidArtifact):
                    validate_vacuum_output(t, t['baseline_output']['path'])
            self.assertEqual({(t['case'], t['metal']): t['charge'] for t in m['tasks']},
                             {('1h4i_qm33', 'La'): -1, ('1h4i_qm33', 'Ca'): -2,
                              ('1h4i_qm36', 'La'): -2, ('1h4i_qm36', 'Ca'): -3})
            with self.assertRaises(InvalidArtifact):
                prepare_partition(PILOT, ENV, out, AGREEMENT)

    def test_unexpected_hamiltonian_rejected_on_corrupted_real_input(self):
        t = next(t for t in read_json(PILOT)['tasks'] if t['case'] == '1h4i_qm33')
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'explicitly_corrupted_real_input.inp'
            p.write_text(verify(t['input']).read_text().replace('r2SCAN-3c', 'B3LYP'))
            with self.assertRaises(InvalidArtifact):
                input_state(p, cpcm=True)


if __name__ == '__main__':
    unittest.main()
