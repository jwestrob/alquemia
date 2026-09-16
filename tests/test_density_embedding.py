"""Pinned real-density and prepared-state checks; no solver execution."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from density_embedding import parse_potential


class RealDensityInputs(unittest.TestCase):
    def test_existing_quantum_potential_is_parsed_in_actual_probe_order(self):
        p = ROOT / 'workspaces/global_electrostatic_20260916/esp_v1/1h4i_qm33_La/quality.json'
        if not p.exists(): self.skipTest('pinned real ESP fixture unavailable')
        verify({'path':str(p),'sha256':'50e5a8438281d29370b6224233a91f3b431a7ffdbb2a6704f1e2caa0aaf14de7'})
        q = read_json(p); points = np.loadtxt(verify(q['points']),skiprows=1)
        potential = parse_potential(verify(q['actual_quantum_potential']),points)
        self.assertEqual(len(potential),q['point_count'])
        self.assertTrue(np.isfinite(potential).all())
        # Explicitly corrupted coordinate copy; no fabricated potential.
        bad = points.copy(); bad[0,0] += .001
        with self.assertRaises(InvalidArtifact): parse_potential(verify(q['actual_quantum_potential']),bad)

    def test_prepared_field_is_exactly_the_existing_environment(self):
        p = ROOT / 'workspaces/density_embedding_20260916/prepared_v1/manifest.json'
        if not p.exists(): self.skipTest('actual prepared field manifest unavailable')
        m = read_json(p)
        self.assertEqual(len(m['tasks']),4)
        for task in m['tasks']:
            state = read_json(verify(task['source_state']))
            actual = np.loadtxt(verify(task['pointcharges']),skiprows=1)
            expected = np.array([[a['charge_e'],*a['xyz_A']] for a in state['environment_atoms']])
            np.testing.assert_array_equal(actual,expected)
            self.assertEqual(task['xyz']['sha256'],task['vacuum_task']['xyz']['sha256'])
            text = verify(task['input']).read_text()
            self.assertIn('DoEQ false',text)
            self.assertNotIn('CPCM',text)
            self.assertNotIn('Opt',text)

    def test_paired_external_charge_files_are_identical(self):
        p = ROOT / 'workspaces/density_embedding_20260916/prepared_v1/manifest.json'
        if not p.exists(): self.skipTest('actual prepared field manifest unavailable')
        tasks = read_json(p)['tasks']
        for case in ('1h4i_qm33','1h4i_qm36'):
            hashes = {t['pointcharges']['sha256'] for t in tasks if t['case']==case}
            self.assertEqual(len(hashes),1)


if __name__ == '__main__': unittest.main()
