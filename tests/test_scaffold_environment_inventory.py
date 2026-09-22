"""Actual source/topology fixtures; no molecular energy or force evaluations."""
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import openmm as mm
from openmm import app, unit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
import scaffold_environment_inventory as inventory


class RealScaffoldInventory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = read_json(ROOT/'workspaces/scaffold_environment_20260922/inventory_v2/INVENTORY.json')
        cls.config = read_json(verify(cls.result['config']))
        cls.manifest = read_json(verify(cls.config['source_manifest']))

    def test_exact_source_atoms_protons_and_explicit_terminal_completion(self):
        self.assertEqual([r['case_id'] for r in self.result['rows']], list(inventory.CASES))
        for row in self.result['rows']:
            pdb = app.PDBFile(str(verify(row['source'])))
            pos = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
            actual = {inventory.atom_id(a): (a, pos[a.index]) for a in pdb.topology.atoms() if a.residue.name in inventory.STANDARD}
            saved = read_json(verify(row['artifacts']['parent_atoms.json']))
            state = read_json(verify(row['artifacts']['source_proton_inventory.json']))
            self.assertFalse(state['source_H_moved'])
            self.assertFalse(state['global_archived_H_normalization_reused'])
            originals = [a for a in saved if a['source_or_archived_completion'] == 'source']
            self.assertEqual({a['source_id'] for a in originals}, set(actual))
            for a in originals:
                source, p = actual[a['source_id']]
                self.assertEqual(a['element'], source.element.symbol)
                np.testing.assert_array_equal(a['xyz_A'], p)
            self.assertEqual(len(saved)-len(originals), int(row['case_id'] in ('1H4I', '4MAE')))

    def test_actual_paired_context_and_every_cap_have_real_source_bond(self):
        for row in self.result['rows']:
            maps = [read_json(verify(row['endpoints'][z]['mapping'])) for z in ('Ca', 'La')]
            self.assertEqual(maps[0], maps[1])
            mapping = read_json(verify(row['artifacts']['context_mapping.json']))
            atoms = read_json(verify(row['artifacts']['parent_atoms.json']))
            real = np.asarray([a['xyz_A'] for a in atoms])
            actual = np.asarray([a[1:] for a in xyz(verify(row['endpoints']['Ca']['xyz']))])
            np.testing.assert_array_equal(actual, maps[0]['context']['core_positions_A'])
            for link in mapping['source_context_atoms']:
                # OpenMM's nm->Angstrom conversion has last-bit floating error;
                # use the inventory's predeclared 1e-12 A mapping criterion.
                np.testing.assert_allclose(actual[link['context_index']], real[link['parent_index']], atol=1e-12, rtol=0)
            caps = sorted(tuple(sorted((c['retained_parent_index'], c['omitted_parent_index']))) for c in mapping['caps'])
            self.assertEqual(caps, [tuple(p) for p in mapping['actual_cut_bonds']])
            self.assertGreater(len(row['local_capped_unmatched_residues']), 0)

    def test_serialized_parent_parameters_and_partition_are_complete(self):
        for row in self.result['rows']:
            system = mm.XmlSerializer.deserialize(verify(row['artifacts']['parent_system.xml']).read_text())
            self.assertEqual(system.getNumParticles(), row['parent_atom_count'])
            self.assertEqual(system.getNumConstraints(), 0)
            terms = read_json(verify(row['artifacts']['term_supports.json']))
            local = set(read_json(verify(row['artifacts']['context_mapping.json']))['local_parent_indices'])
            for force in terms:
                self.assertEqual(force['type'], type(system.getForce(force['index'])).__name__)
                self.assertEqual(len(force['supports']), len(force['partitions']))
                for support, part in zip(force['supports'], force['partitions']):
                    self.assertEqual(part, inventory.partition(support, local))
            nb = next(r for r in row['parameterized_force_supports'] if r['type'] == 'NonbondedForce')
            n = row['parent_atom_count']
            self.assertEqual(sum(nb['ordinary_pair_counts'].values())+nb['support_count'], n*(n-1)//2)
            self.assertAlmostEqual(nb['charge_e'], round(nb['charge_e']), places=8)
            self.assertFalse(row['absolute_or_relaxation_score_available'])

    def test_corrupted_real_mapping_missing_cap_is_rejected(self):
        manifest = deepcopy(self.manifest)
        tasks = [t for t in manifest['tasks'] if t['case_id'] == '1H4I']
        mapping = read_json(verify(tasks[0]['mapping']))
        links = mapping['context']['core_links']
        links.remove(next(l for l in links if l[1] == 'cap'))
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); path = tmp/'explicitly_corrupted_real_mapping.json'
            write_new(path, mapping)
            for task in tasks: task['mapping'] = record(path)
            with self.assertRaisesRegex(InvalidArtifact, 'real cut bonds differ from cap coverage'):
                inventory.inspect_case('1H4I', manifest, self.config, tmp/'result')

    def test_no_silent_unknown_cofactor_omission(self):
        source = verify(self.result['rows'][0]['source'])
        lines = source.read_text().splitlines(keepends=True); changes = 0
        for i, line in enumerate(lines):
            if line.startswith(('ATOM  ', 'HETATM')) and line[17:20] == 'PQQ':
                lines[i] = line[:17]+'UNK'+line[20:]; changes += 1
        self.assertGreater(changes, 0)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'explicitly_corrupted_real_PQQ.pdb'; path.write_text(''.join(lines))
            with self.assertRaisesRegex(InvalidArtifact, 'unsupported source residues'):
                inventory.parent_model(path, verify(self.config['forcefield']), None)


if __name__ == '__main__':
    unittest.main()
