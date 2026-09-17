"""Checks on saved scientific artifacts; these tests run no model inference."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from mace_output_audit import load_pair

AUDIT = ROOT / 'workspaces/mace_large_20260916/output_audit_v1/result.json'


@unittest.skipUnless(AUDIT.exists(), 'requires completed real medium/large outputs')
class SavedOutputTests(unittest.TestCase):
    def test_direct_gradient_sign_and_physical_mapping(self):
        r = read_json(AUDIT)
        for name, ref in r['collections'].items():
            _, atoms, pair = load_pair(verify(ref))
            gradient = np.load(verify(r['models'][name]['gradient']))
            # Floating subtraction followed by addition is not bitwise invertible.
            scale = np.maximum(1., np.abs(gradient) + np.abs(pair['Ca']['forces']))
            error = np.abs(gradient + pair['Ca']['forces'] - pair['La']['forces'])
            self.assertTrue(np.all(error <= 4 * np.finfo(float).eps * scale))
            self.assertEqual(len(atoms), r['atoms'])
            self.assertIsNone(r['hybrid_gradient'])
            self.assertIsNone(r['relaxation_correction'])

    def test_all_atoms_charge_and_component_accounting(self):
        r = read_json(AUDIT)
        for m in r['models'].values():
            self.assertEqual(sum(s['atoms'] for s in m['shells']), r['atoms'])
            self.assertAlmostEqual(sum(s['Ca_minus_La_charge_sum_e'] for s in m['shells']), -1., places=9)
            self.assertAlmostEqual(sum(m['R_components_kcal_mol'].values()), m['raw_R_kcal_mol'], places=8)
            gradient = np.load(verify(m['gradient']))
            self.assertAlmostEqual(sum(s['direct_R_gradient_squared_norm_eV2_A2'] for s in m['shells']), float(np.square(gradient).sum()), places=8)

    def test_corrupted_real_coordinate_order_rejected(self):
        ref = read_json(AUDIT)['collections']['large']
        from mace_output_audit import xyz
        def reversed_coordinates(path):
            return list(reversed(xyz(path)))
        with patch('mace_output_audit.xyz', side_effect=reversed_coordinates):
            with self.assertRaisesRegex(InvalidArtifact, 'coordinate/element'):
                load_pair(verify(ref))


if __name__ == '__main__':
    unittest.main()
