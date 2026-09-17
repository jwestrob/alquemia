"""Actual parameter artifacts and deliberately corrupted real preparations.

No synthetic molecular data and no energy/force evaluations.
"""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, write_new
from mace_amoeba_capability import topology, validate

W = ROOT/'workspaces/mace_omol_20260917/amoeba_capability_v1'


class AmoebaCapability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'result.json').exists():
            raise unittest.SkipTest('real executed AMOEBA framework preparations unavailable')
        cls.result = read_json(W/'result.json')
        cls.prep = read_json(verify(cls.result['cases'][0]['source_preparation']))

    def test_executed_frameworks_preserve_physical_coordinates_and_charge(self):
        self.assertEqual(self.result['framework_passes'], 3)
        self.assertFalse(self.result['full_hybrid_qualified'])
        for row in self.result['cases']:
            prep = read_json(verify(row['source_preparation']))
            mapping = read_json(verify(row['mapping']))
            verify(row['framework'])
            self.assertEqual(mapping['physical_atoms'], prep['physical_atoms'])
            physical = {a['id']: a for a in prep['physical_atoms']}
            np.testing.assert_array_equal(mapping['positions_A'],
                [physical[i]['xyz_A'] for i in mapping['system_atom_ids']])
            self.assertEqual(set(mapping['system_atom_ids'])|{'metal'}, set(physical))
            self.assertEqual(row['unparameterized_atoms'], 1)
            self.assertEqual(row['new_energy_calls'], 0)
            self.assertEqual(row['new_force_calls'], 0)
            self.assertAlmostEqual(row['protein_charge_e'], prep['protein_charge_e'], places=5)
            for param in mapping['parameters']:
                self.assertAlmostEqual(param['multipole_md_units'][0], param['gk_md_units'][0], places=12)

    def test_duplicate_real_atom_rejected(self):
        prep = copy.deepcopy(self.prep)
        prep['physical_atoms'].append(copy.deepcopy(prep['physical_atoms'][0]))
        with self.assertRaisesRegex(InvalidArtifact, 'duplicate physical identity'):
            topology(prep)

    def test_modified_real_heavy_coordinate_rejected(self):
        prep = copy.deepcopy(self.prep)
        atom = next(a for a in prep['physical_atoms'] if a['element'] == 'C')
        atom['xyz_A'][0] += 0.1  # explicitly corrupted fixture, not a displacement experiment
        with self.assertRaisesRegex(InvalidArtifact, 'paired endpoint coordinates'):
            topology(prep)

    def test_unsupported_cofactor_not_silently_dropped(self):
        prep = copy.deepcopy(self.prep)
        prep['physical_atoms'][0]['kind'] = 'unrecognized_cofactor'
        with self.assertRaisesRegex(InvalidArtifact, 'unsupported cofactor'):
            topology(prep)

    def test_missing_real_bond_rejected(self):
        prep = copy.deepcopy(self.prep)
        prep['preparation_details']['bonds'].pop()
        with self.assertRaisesRegex(InvalidArtifact, 'connectivity/disulfides'):
            topology(prep)

    def test_manifest_change_invalidates_cache(self):
        self.assertEqual(validate(W/'manifest.json')['tasks'], 3)
        manifest = read_json(W/'manifest.json')
        manifest['settings']['polarization'] = 'direct'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'corrupted_manifest.json'
            write_new(path, manifest)
            with self.assertRaisesRegex(InvalidArtifact, 'manifest key differs'):
                validate(path)


if __name__ == '__main__':
    unittest.main()
