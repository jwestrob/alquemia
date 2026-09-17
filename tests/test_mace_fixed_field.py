"""Real archived fixed-field components; no fabricated energies or inference."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from mace_fixed_field import archive_equal, contrast, coordinate_check, direct_coupling
from mace_hybrid import rotation
RESULT = ROOT/'workspaces/mace_fixed_field_20260916/screen_v1/result.json'


@unittest.skipUnless(RESULT.exists(), 'requires completed real fixed-field screen')
class FixedFieldTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = read_json(RESULT)
        cls.fit = read_json(verify(cls.result['charge_fit']))

    def test_components_reproduce_recorded_partition_and_full_cancellation(self):
        for model in self.result['models'].values():
            pairs = {size: contrast(model['rows'][f'1h4i_{size}_Ca'], model['rows'][f'1h4i_{size}_La'])
                     for size in ('qm33', 'qm36')}
            delta = pairs['qm36']['R_vacuum_screen_kcal_mol']-pairs['qm33']['R_vacuum_screen_kcal_mol']
            self.assertAlmostEqual(delta, model['partition']['R_vacuum_screen_kcal_mol'], places=9)
            self.assertEqual(model['partition']['short_full_R_kcal_mol'], 0.)
            self.assertIsNone(pairs['qm33']['solution_score'])
            self.assertIsNone(pairs['qm36']['calibrated_class'])
        self.assertTrue(self.result['accounting_checks_pass'])
        self.assertFalse(self.result['primary_partition_gate_pass'])

    def test_actual_fitted_charges_replay_coupling_and_rigid_invariance(self):
        r = rotation()
        for row in self.fit['rows']:
            source = self.result['models']['medium']['rows'][row['task_id']]['source_state']
            state = read_json(verify(source)); value = direct_coupling(state, row['charge_e'])
            self.assertAlmostEqual(value, row['direct_fit_kcal_mol'], places=8)
            transformed = copy.deepcopy(state)
            for key in ('core_atoms', 'environment_atoms'):
                for atom in transformed[key]:
                    atom['xyz_A'] = (r @ np.array(atom['xyz_A'])+np.array([4., -7., 2.])).tolist()
            self.assertAlmostEqual(value, direct_coupling(transformed, row['charge_e']), places=8)

    def test_corrupted_actual_atom_identity_cannot_pass_mapping(self):
        name = next(iter(self.result['models']['medium']['rows']))
        state = read_json(verify(self.result['models']['medium']['rows'][name]['source_state']))
        atoms = [(a['element'], *a['xyz_A']) for a in state['core_atoms']]
        coordinate_check(atoms, state['core_atoms'], state['metal'])
        corrupted = list(atoms); corrupted[1] = ('H', *corrupted[1][1:])
        with self.assertRaises(InvalidArtifact):
            coordinate_check(corrupted, state['core_atoms'], state['metal'])

    def test_replay_allows_roundoff_but_rejects_altered_real_coupling_or_metadata(self):
        actual = self.fit['rows'][0]
        rounded = copy.deepcopy(actual); rounded['direct_fit_kcal_mol'] += 1e-12
        self.assertTrue(archive_equal(actual, rounded))
        corrupted = copy.deepcopy(actual); corrupted['direct_fit_kcal_mol'] += .01
        self.assertFalse(archive_equal(actual, corrupted))
        corrupted = copy.deepcopy(actual); corrupted['task_id'] = 'corrupted_fixture'
        self.assertFalse(archive_equal(actual, corrupted))


if __name__ == '__main__':
    unittest.main()
