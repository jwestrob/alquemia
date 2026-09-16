"""Stage-B integrity tests use only pinned real GGR structures and corrupt copies."""
from pathlib import Path
import sys
import tempfile
import unittest

import gemmi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, paired, xyz
from ggr_structure_prepare import inspect_source


class GGRSourceIntegrity(unittest.TestCase):
    def setUp(self):
        self.sources = ROOT / 'workspaces/ggr_mechanism_20260915/stage_b_sources'
        self.reference = ROOT / 'diagnostics/laca_benchmark_expansion_20260915/sources/pdb/1GLG.cif'
        self.repair = ROOT / 'workspaces/affordable_challenger_20260915/verified_repairs/ggr_1glg_GGR/repair_manifest.json'
        if not self.repair.is_file() or any(not (self.sources / (sid+'.cif')).is_file() for sid in ('2FW0','2FVY')):
            self.skipTest('pinned real GGR source/repair fixtures unavailable')
        self.topology = verify(read_json(self.repair)['topology_definition'])

    def inspect(self, source):
        return inspect_source(source, self.reference, self.repair, self.topology)

    def test_real_sources_recover_frozen_composition_without_forced_selection(self):
        for sid in ('2FW0', '2FVY'):
            row = self.inspect(self.sources / (sid+'.cif'))
            self.assertEqual(row['status'], 'source_integrity_pass')
            self.assertEqual(row['sequence_length'], 309)
            self.assertEqual(row['missing_sequence_positions'], [1,307,308,309])
            self.assertEqual(row['actual_donor_sequence_positions'], [134,136,138,140,142,205])
            self.assertEqual(len(row['typed_qualification']), 7)
            self.assertEqual(row['explicit_water_inventory'], [])
            self.assertEqual(row['local_missing_heavy_atoms'], [])
            self.assertFalse(any(x['resname']=='HIS' for x in row['nonstandard_species']))

    def test_radiation_and_crystal_solutes_remain_visible(self):
        row = self.inspect(self.sources / '2FVY.cif')
        self.assertIn('Glu149', row['radiation_damage_caveat'])
        glu149 = next(r for r in row['alternate_conformers'] if r['resnum']==149)
        self.assertFalse(glu149['in_local_representation_or_cap_neighbours'])
        self.assertEqual({a['altloc'] for a in glu149['selected_atoms']} - {''}, {'A'})
        self.assertEqual({x['resname'] for x in row['nonstandard_species']}, {'BGC','ACT','CO2','GOL'})
        self.assertEqual({x['resname'] for x in self.inspect(self.sources/'2FW0.cif')['nonstandard_species']}, {'NA','CIT','MLA'})

    def test_sequence_mapping_does_not_assume_author_numbering(self):
        with tempfile.TemporaryDirectory() as td:
            st = gemmi.read_structure(str(self.sources/'2FW0.cif'))
            for residue in st[0]['A']:
                if residue.label_seq is not None:
                    residue.seqid = gemmi.SeqId(residue.seqid.num+1000, residue.seqid.icode)
            altered = Path(td)/'real_2FW0_author_identifiers_renumbered.cif'
            st.make_mmcif_document().write_file(str(altered))
            row = self.inspect(altered)
            self.assertEqual(row['status'], 'source_integrity_pass')
            self.assertEqual(row['actual_donor_sequence_positions'], [134,136,138,140,142,205])
            self.assertEqual({x['resnum'] for x in row['typed_inclusion']}, {1134,1136,1138,1140,1142,1205})

    def test_missing_real_amide_nitrogen_is_unsupported(self):
        with tempfile.TemporaryDirectory() as td:
            st = gemmi.read_structure(str(self.sources/'2FW0.cif'))
            residue = next(r for r in st[0]['A'] if r.label_seq==141)
            residue.remove_atom('N', '\x00', gemmi.Element('N'))
            altered = Path(td)/'explicitly_corrupted_real_2FW0_missing_amide_N.cif'
            st.make_mmcif_document().write_file(str(altered))
            with self.assertRaisesRegex(InvalidArtifact, 'amide N'):
                self.inspect(altered)

    def test_prepared_pairs_share_one_frozen_source_per_structure(self):
        inventory = ROOT/'workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/inventory.json'
        if not inventory.is_file():
            self.skipTest('actual Stage B preparations unavailable; no simulated success')
        data = read_json(inventory)
        self.assertEqual(data['status'], 'prepared')
        self.assertEqual(len(data['preparations']), 4)
        sources = {}
        for row in data['preparations']:
            manifest = read_json(verify(row['manifest']))
            endpoints = manifest['outputs']
            result = paired(verify(endpoints['La']['xyz']), verify(endpoints['Ca']['xyz']),
                            endpoints['La']['charge'], endpoints['Ca']['charge'])
            self.assertEqual(result['atom_count'], 52 if row['representation']=='amide_v3' else 58)
            self.assertEqual((endpoints['La']['charge'], endpoints['Ca']['charge']), (0,-1))
            self.assertEqual(manifest['explicit_water_inventory'], [])
            self.assertEqual(row['source_coordinate_max_displacement_A'], 0.0)
            self.assertEqual(manifest['source_structure'], row['source_protonation'])
            sources.setdefault(row['pdb_id'], []).append(row['source_protonation'])
        for reused in sources.values():
            self.assertEqual(reused[0], reused[1])

    def test_all_prepared_source_heavy_atoms_match_deposited_selected_coordinates(self):
        inventory = ROOT/'workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/inventory.json'
        if not inventory.is_file():
            self.skipTest('actual Stage B preparations unavailable; no simulated success')
        data = read_json(inventory)
        cases = {r['pdb_id']: r for r in data['cases']}
        for row in data['preparations']:
            selection = read_json(verify(cases[row['pdb_id']]['source_selection']))['atoms']
            original = {(a['chain'],a['resnum'],a['insertion_code'],a['atom']):a
                        for a in selection if a['element'] not in ('H','D')}
            manifest = read_json(verify(row['manifest']))
            coords = xyz(verify(manifest['outputs']['La']['xyz']))
            for atom in manifest['atom_graph']['source_to_qm']:
                if atom['kind']!='source' or atom['source']['element'] in ('H','D'):
                    continue
                s = atom['source']
                key = (s['chain'],s['resnum'],s['insertion_code'],s['atom'])
                self.assertEqual(list(coords[atom['qm_index']][1:]), original[key]['xyz_A'])
        for case in cases.values():
            self.assertEqual(case['protonation']['ph'], 7.0)
            self.assertFalse(case['protonation']['add_missing_residues'])
            self.assertEqual(case['protonation']['repaired_missing_atom_count'], 0)
            self.assertEqual(case['heavy_coordinate_audit']['maximum_displacement_A'], 0.0)
            self.assertEqual(case['heavy_coordinate_audit']['added_heavy_atom_identifiers'], [['A',306,'','OXT']])


if __name__ == '__main__':
    unittest.main()
