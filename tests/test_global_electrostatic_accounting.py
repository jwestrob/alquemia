"""Read-only checks of real archived ownership and arithmetic, never new science.

The arithmetic fixtures deliberately come from the older CPCM/APBS protocol.
Their numbers exercise signs and units only; they are not vacuum/TABI scores.
No structure, charge distribution, energy, or solver success is fabricated.
"""
from pathlib import Path
import math
import sys
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import HA_TO_KCAL, InvalidArtifact, energy, read_json, record, verify, xyz
from affordable_environment import physical_boundary_key
from global_electrostatic import endpoint_components, paired_components

BASE = ROOT / 'workspaces/affordable_challenger_20260915'
PINS = {
    '1h4i_qm33_La': 'b6597141c4143f854b23f2d9cfdcfce0e4c973cee08b597c81550d72c5fa9c8a',
    '1h4i_qm33_Ca': '25f7b78bca9e000b263ce2514f6999cb248f79912f1e0b7636507e3a288df952',
    '1h4i_qm36_La': '12e7496aece17805dc0173f714d4f9b1beb1aee535632551b3591331585d6c4b',
    '1h4i_qm36_Ca': 'f2f2d39cc9109383814a399678412ca15687f8926e983934a256b0b98ebcd971',
}
BOUNDARY_KEY = '4c518676f0b57cd3d30b135734e226e49b98627f365908f668cd7e27c7ecb1c4'


class ArchivedBoundaryOwnership(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (BASE / 'environment').is_dir():
            raise unittest.SkipTest('real archived qm33/qm36 environment fixtures unavailable')
        cls.states = {}
        for name, sha in PINS.items():
            path = BASE / 'environment' / name / 'skeleton.json'
            verify({'path': str(path), 'sha256': sha})
            cls.states[name] = read_json(path)

    def test_pinned_sources_and_existing_coordinates(self):
        for state in self.states.values():
            verify(state['source'])
            boundary = read_json(verify(state['boundary_mapping']))
            verify(boundary['forcefield'])
            self.assertEqual(boundary['source'], state['source'])
            coordinates = xyz(verify(state['source_endpoint_xyz']))
            self.assertEqual(len(coordinates), len(state['core_atoms']))
            for row, atom in zip(coordinates, state['core_atoms']):
                self.assertEqual(row[0], atom['element'])
                self.assertEqual(list(row[1:]), atom['xyz_A'])

    def test_common_source_cavity_and_declared_chain(self):
        serializations = []
        for state in self.states.values():
            self.assertEqual(physical_boundary_key(state['physical_atoms']), BOUNDARY_KEY)
            self.assertEqual(len(state['physical_atoms']), 9141)
            self.assertEqual(state['assembly'], 'deposited_catalytic_chain_A')
            self.assertEqual(state['explicit_waters'], [])
            boundary = read_json(verify(state['boundary_mapping']))
            self.assertEqual(boundary['explicitly_excluded_other_chains'], list('BCDEF'))
            self.assertEqual(boundary['source_protein_atom_count'], 9113)
            self.assertFalse(any(a['id'].startswith('cap/') for a in state['physical_atoms']))
            serializations.append(sorted((a['id'], *(f'{v:.10f}' for v in a['xyz_A']),
                                           f"{a['radius_A']:.6f}")
                                          for a in state['physical_atoms']))
        self.assertTrue(all(s == serializations[0] for s in serializations))

    def test_paired_state_and_whole_system_charge(self):
        for n, core_count, env_count, la_charge, env_charge in (
                (33, 47, 9094, -1, -7), (36, 54, 9087, -2, -6)):
            la, ca = (self.states[f'1h4i_qm{n}_{m}'] for m in ('La', 'Ca'))
            for key in ('source', 'assembly', 'microstate', 'explicit_waters',
                        'environment_atoms', 'physical_atoms'):
                self.assertEqual(la[key], ca[key])
            for state, charge in ((la, la_charge), (ca, la_charge - 1)):
                self.assertEqual(len(state['core_atoms']), core_count)
                self.assertEqual(len(state['environment_atoms']), env_count)
                self.assertEqual(state['core_total_charge_e'], charge)
                self.assertAlmostEqual(sum(a['charge_e'] for a in state['environment_atoms']), env_charge)
                self.assertAlmostEqual(charge + state['expected_environment_charge_e'],
                                       -8 if state['metal'] == 'La' else -9)
            for a, b in zip(la['core_atoms'], ca['core_atoms']):
                for key in ('id', 'kind', 'xyz_A', 'radius_A'):
                    self.assertEqual(a[key], b[key])

    def test_source_ownership_and_residue_charge_closure_against_ff(self):
        for state in self.states.values():
            boundary = read_json(verify(state['boundary_mapping']))
            templates = {r.attrib['name']: {a.attrib['name']: float(a.attrib['charge'])
                                          for a in r.findall('Atom')}
                         for r in ET.parse(verify(boundary['forcefield'])).findall('./Residues/Residue')}
            physical = {a['id']: a for a in state['physical_atoms']}
            core = {a['id']: a for a in state['core_atoms']}
            env = {a['id']: a for a in state['environment_atoms']}
            self.assertEqual(len(physical), len(state['physical_atoms']))
            self.assertEqual(len(core), len(state['core_atoms']))
            self.assertEqual(len(env), len(state['environment_atoms']))
            self.assertFalse(set(core) & set(env))
            expected_cavity_only = set()
            for ledger in boundary['ledgers']:
                fragment = ledger['fragment_id'].split(':')[1]
                residue_name, residue_number = fragment[:3], fragment[3:]
                prefix = f'A/{residue_number}/ /'
                names = templates[residue_name]
                removed_ids = {state['physical_atoms'][i]['id'] for i in ledger['removed_source_indices']}
                source_ids = {state['physical_atoms'][i]['id'] for i in ledger['source_qm_indices']}
                self.assertEqual(removed_ids - source_ids, {prefix + 'CA'})
                self.assertTrue(source_ids <= set(core))
                self.assertFalse(removed_ids & set(env))
                expected_cavity_only.add(prefix + 'CA')
                recipients = {state['physical_atoms'][i]['id'] for i in ledger['recipients']}
                self.assertEqual(recipients, {prefix + 'N', prefix + 'C'})
                for atom_id in recipients:
                    self.assertAlmostEqual(env[atom_id]['charge_e'],
                                           names[atom_id.rsplit('/', 1)[1]] + ledger['each_increment_e'])
                remaining_charge = sum(a['charge_e'] for i, a in env.items() if i.startswith(prefix))
                self.assertAlmostEqual(remaining_charge + ledger['formal_charge'], sum(names.values()))
            self.assertEqual(set(physical) - set(core) - set(env), expected_cavity_only)

    def test_partition_moves_only_recorded_asp303_ownership(self):
        small = self.states['1h4i_qm33_La']; large = self.states['1h4i_qm36_La']
        se = {a['id']: a for a in small['environment_atoms']}
        le = {a['id']: a for a in large['environment_atoms']}
        removed = {f'A/303/ /{name}' for name in ('CA', 'CB', 'CG', 'HB2', 'HB3', 'OD1', 'OD2')}
        self.assertEqual(set(se) - set(le), removed)
        self.assertFalse(set(le) - set(se))
        changes = {i: le[i]['charge_e'] - se[i]['charge_e']
                   for i in set(se) & set(le) if le[i]['charge_e'] != se[i]['charge_e']}
        self.assertEqual(set(changes), {'A/303/ /N', 'A/303/ /C'})
        for difference in changes.values():
            self.assertAlmostEqual(difference, .09)
        sc = {a['id'] for a in small['core_atoms']}
        lc = {a['id'] for a in large['core_atoms']}
        self.assertEqual(lc - sc, (removed - {'A/303/ /CA'}) | {'cap/A:ASP303'})

    def test_caps_are_inside_source_cavity_and_not_on_environment_charges(self):
        expected = {'cap/A:GLU177': (0.0561118804045424, 1.4374210582560665),
                    'cap/A:ASN261': (0.05442971700191357, 1.453939387007239),
                    'cap/A:ASP303': (0.05756479663910424, 1.3668652637231622)}
        for state in self.states.values():
            for atom in state['core_atoms']:
                closest = min(math.dist(atom['xyz_A'], e['xyz_A']) for e in state['environment_atoms'])
                self.assertGreaterEqual(closest, 1.0)
                if atom['kind'] != 'cap':
                    continue
                margin = max(p['radius_A'] - atom['radius_A'] - math.dist(atom['xyz_A'], p['xyz_A'])
                             for p in state['physical_atoms'])
                self.assertGreater(margin, 0.)
                self.assertAlmostEqual(margin, expected[atom['id']][0])
                self.assertAlmostEqual(closest, expected[atom['id']][1])


class ArchivedValueArithmeticOnly(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = BASE / 'solver_completion/solver_result.json'
        if not path.is_file():
            raise unittest.SkipTest('real archived APBS arithmetic fixtures unavailable')
        verify({'path': str(path), 'sha256': '542d491159111114e54ba81df220b08f7bbf7d8db619e2a32b970dea3b0ee885'})
        results = {r['label'].split('/')[0]: r['result'] for r in read_json(path)['numerical_checks']
                   if r['label'].endswith('/primary')}
        pilot = BASE / 'pilot/pilot_manifest.json'
        verify({'path': str(pilot), 'sha256': '91ab53de327806f3e042791957f72c34ace68a9368eef0f99601efdfa4dc987a'})
        cls.values = {}
        for task in read_json(pilot)['tasks']:
            if task['task_id'] not in results:
                continue
            result = results[task['task_id']]
            verify(result['source_manifest'])
            charging = result['charging_energies_kJ_mol']
            cls.values[task['task_id']] = (
                energy(task['output_path']), result['components']['direct_core_environment_kcal_mol'],
                charging['target_total'] - charging['homogeneous_total'])

    def test_actual_archived_values_convert_once_without_old_counterterms(self):
        for values in self.values.values():
            got = endpoint_components(*values)
            e, direct, rf = values
            self.assertEqual(got['reaction_field_kcal_mol'], rf / 4.184)
            self.assertEqual(got['total_kcal_mol'], e * HA_TO_KCAL + direct + rf / 4.184)
            self.assertNotIn('reaction_field_reference_core_kcal_mol', got)
            self.assertNotIn('reaction_field_environment_kcal_mol', got)

    def test_real_values_obey_ca_minus_la_and_no_reference_is_invented(self):
        for n in (33, 36):
            ca, la = (self.values[f'1h4i_qm{n}_{m}'] for m in ('Ca', 'La'))
            result = paired_components(endpoint_components(*ca), endpoint_components(*la))
            expected = (ca[0] - la[0]) * HA_TO_KCAL + ca[1] - la[1] + (ca[2] - la[2]) / 4.184
            self.assertAlmostEqual(result['R_global_kcal_mol'], expected, places=8)
            self.assertIsNone(result['S_global_kcal_mol'])
            self.assertEqual(result['reference_status'], 'unavailable')
            reversed_result = paired_components(endpoint_components(*la), endpoint_components(*ca))
            self.assertAlmostEqual(reversed_result['R_global_kcal_mol'], -result['R_global_kcal_mol'])

    def test_corrupted_real_components_remain_missing_or_fail(self):
        # Explicitly corrupt real scalar inputs; never create fake successful solver output.
        values = self.values['1h4i_qm33_La']
        complete = endpoint_components(*values)
        for i in range(3):
            missing = list(values); missing[i] = None
            got = endpoint_components(*missing)
            self.assertIsNone(got['total_kcal_mol'])
            self.assertEqual(paired_components(complete, got)['status'], 'unavailable')
            for bad in (math.nan, math.inf):
                corrupt = list(values); corrupt[i] = bad
                with self.assertRaises(InvalidArtifact):
                    endpoint_components(*corrupt)


if __name__ == '__main__':
    unittest.main()
