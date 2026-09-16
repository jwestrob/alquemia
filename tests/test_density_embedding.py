"""Pinned real-density and prepared-state checks; no solver execution."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from density_embedding import parse_potential, parse_chelpg


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

    def test_actual_recovered_density_couplings_close_by_residue(self):
        p = ROOT / 'workspaces/density_embedding_20260916/potential_retry_v1/potentials_1199983.json'
        if not p.exists(): self.skipTest('actual recovered vpot results unavailable')
        result = read_json(p)
        self.assertEqual(result['status'],'complete')
        self.assertEqual(len(result['rows']),4)
        for row in result['rows']:
            receipt = read_json(verify(row['execution_receipt']))
            verify(receipt['task']['density_info'])
            self.assertEqual(receipt['returncode'],0)
            self.assertAlmostEqual(sum(v['density_kcal_mol'] for v in row['per_residue'].values()),
                                   row['direct_density_kcal_mol'],places=9)
            self.assertAlmostEqual(sum(v['mbis_kcal_mol'] for v in row['per_residue'].values()),
                                   row['direct_mbis_atomic_units_kcal_mol'],places=9)

    def test_real_chelpg_charge_order_closure_and_native_controls(self):
        p=ROOT/'workspaces/density_embedding_20260916/chelpg_v1/manifest.json'
        if not p.exists():self.skipTest('real CHELPG task manifest unavailable')
        for task in read_json(p)['tasks']:
            receipt_path=Path(task['directory'])/'execution.json'
            if not receipt_path.exists():self.skipTest('actual CHELPG calculation unavailable')
            receipt=read_json(receipt_path);state=read_json(verify(task['source_potential_task']['state']))
            self.assertEqual(receipt['returncode'],0)
            charge=parse_chelpg(verify(receipt['log']).read_text(),[a['element'] for a in state['core_atoms']],state['core_total_charge_e'])
            self.assertLess(abs(sum(charge)-state['core_total_charge_e']),5e-5)

    def test_corrupted_real_chelpg_settings_are_rejected(self):
        p=ROOT/'workspaces/density_embedding_20260916/chelpg_v1/1h4i_qm33_La/execution.json'
        if not p.exists():self.skipTest('real CHELPG output unavailable')
        receipt=read_json(p);task=receipt['task'];state=read_json(verify(task['source_potential_task']['state']))
        text=verify(receipt['log']).read_text()
        # Corrupt the echoed setting of the real output, not a generated energy.
        self.assertIn('0.300000',text)
        with self.assertRaises(InvalidArtifact):
            parse_chelpg(text.replace('0.300000','0.600000'),[a['element'] for a in state['core_atoms']],state['core_total_charge_e'])


if __name__ == '__main__': unittest.main()
