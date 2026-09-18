"""Source-backed multisite grouping tests; no fabricated molecular evaluations."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from mace_group_transfer import prepared, family_cases, specs, task, validate
from mace_site_groups import build

PREP = ROOT / 'workspaces/mace_group_transfer_20260918/prepared_v2/preparation.json'
MANIFEST = PREP.parent.parent / 'model_v3/manifest.json'


@unittest.skipUnless(PREP.exists(), 'pinned real multisite source preparations unavailable')
class TransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p, cls.cfg, cls.cases = prepared(PREP)

    def test_all_real_sites_have_same_all_calcium_groups_and_paired_sources(self):
        self.assertEqual(len(self.cases), 34)
        self.assertEqual(len(specs(self.cfg)), 55)
        for names in family_cases(self.cfg).values():
            first = self.cases[names[0]]['primary']
            for name in names:
                row = self.cases[name]; g = row['primary']
                self.assertEqual(g['all_Ca_atoms'], first['all_Ca_atoms'])
                self.assertEqual(g['atom_group_indices'], first['atom_group_indices'])
                delta = np.array(g['endpoint_group_charges_e']['La']) - g['endpoint_group_charges_e']['Ca']
                self.assertEqual(np.flatnonzero(delta).tolist(), [g['selected_group_index']])
                self.assertEqual(delta.sum(), 1)
                ca, ca_atoms = task(name, 'Ca', 'primary', row)
                la, la_atoms = task(name, 'La', 'primary', row)
                np.testing.assert_array_equal(np.array(ca_atoms, dtype=object)[:, 1:], np.array(la_atoms, dtype=object)[:, 1:])
                changed = [i for i, (a, b) in enumerate(zip(ca_atoms, la_atoms)) if a[0] != b[0]]
                self.assertEqual(changed, [g['selected_metal_index']])
                self.assertEqual(la['charge'] - ca['charge'], 1)
                self.assertEqual(ca['state']['all_electron_count'] % 2, 0)
                self.assertEqual(la['state']['all_electron_count'] % 2, 0)

    def test_actual_amide_and_shared_metal_groups_use_source_bonds(self):
        row = self.cases['GGR_1GLG']; p = read_json(verify(row['case']['preparation']))
        g = row['primary']; assignment = dict(zip(g['atom_ids'], g['atom_group_indices']))
        adjacency = {}
        for b in p['preparation_details']['bonds']:
            a, c = b['atom_a_id'], b['atom_b_id']
            adjacency.setdefault(a, set()).add(c); adjacency.setdefault(c, set()).add(a)
        amides = 0
        for c in g['contacts']:
            if c['donor_type'] != 'backbone_O':
                continue
            carbon = next(a for a in adjacency[c['donor_id']] if a.endswith('/C'))
            for n in adjacency[carbon]:
                if n.endswith('/N'):
                    amides += 1
                    self.assertEqual(assignment[n], assignment[c['metal_id']])
        self.assertGreater(amides, 0)
        rtx = self.cases[family_cases(self.cfg)['RTX'][0]]['primary']
        metal_groups = [rtx['atom_group_indices'][i] for i in rtx['all_metal_indices']]
        self.assertLess(len(set(metal_groups)), len(metal_groups))

    def test_explicit_corrupted_real_charge_duplicate_and_water_fail(self):
        row = self.cases['PARV_4CPV_CD']; p = read_json(verify(row['case']['preparation']))
        bonds = [tuple(sorted((b['atom_a_id'], b['atom_b_id']))) for b in p['preparation_details']['bonds']]
        broken = copy.deepcopy(p); broken['endpoints']['La']['charge'] += 1
        with self.assertRaises(InvalidArtifact):
            build(broken, row['protein_ledger'], bonds)
        broken = copy.deepcopy(p); broken['physical_atoms'].append(copy.deepcopy(broken['physical_atoms'][0]))
        with self.assertRaises(InvalidArtifact):
            build(broken, row['protein_ledger'], bonds)
        broken = copy.deepcopy(p)
        index = next(i for i, a in enumerate(broken['physical_atoms']) if a['kind'] == 'retained_site_water' and a['element'] == 'H')
        broken['physical_atoms'].pop(index)
        with self.assertRaises(InvalidArtifact):
            build(broken, row['protein_ledger'], bonds)

    @unittest.skipUnless(MANIFEST.exists(), 'actual finite transfer manifest unavailable')
    def test_actual_manifest_and_rejected_charge_mutation(self):
        self.assertEqual(validate(MANIFEST)['tasks'], 55)
        m = read_json(MANIFEST); m['tasks'][0]['group_charges_e'][0] += 1
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'corrupted_real_manifest.json'; path.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):
                validate(path)

    def test_actual_permutation_preserves_physical_group_state(self):
        name = family_cases(self.cfg)['RTX'][0]; row = self.cases[name]
        a, xyz = task(name, 'La', 'primary', row)
        b, permuted = task(name, 'La', 'permute', row)
        self.assertEqual(xyz, permuted[::-1])
        self.assertEqual(a['atom_group_indices'], b['atom_group_indices'][::-1])
        self.assertEqual(a['source_atom_ids'], b['source_atom_ids'][::-1])
        self.assertEqual(a['charge'], b['charge'])
        self.assertEqual(a['selected_metal_index'], len(xyz) - 1 - b['selected_metal_index'])

    def test_pure_loader_matches_installed_source_definitions_on_real_atoms(self):
        from coordination_policy import donor_type as original
        from mace_site_groups import donor_type as loaded
        for row in self.cases.values():
            p = read_json(verify(row['case']['preparation']))
            for a in p['physical_atoms']:
                if 'resname' in a:
                    args = (a['resname'], a['id'].rsplit('/', 1)[-1], a['element'])
                    self.assertEqual(loaded(*args), original(*args))


if __name__ == '__main__':
    unittest.main()
