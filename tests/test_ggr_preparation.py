"""Real-structure preparation checks; these do not execute ORCA."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify, xyz
from affordable_peptide import SourceGraph
from ggr_preparation import prepare, expand_selection, _source_selection

FIXTURES = {
    'ggr': ('workspaces/affordable_challenger_20260915/verified_repairs/ggr_1glg_GGR/repair_manifest.json', 58),
    'aequorin': ('workspaces/affordable_challenger_20260915/verified_repairs/aequorin_1sl8_EF3/repair_manifest.json', 49),
    'alpha_ca': ('workspaces/benchmark_set_20260915/prepared/alacta_1f6s_v2/strong_site/amide_v3/repair_manifest.json', 52),
    'alpha_la': ('workspaces/benchmark_set_20260915/prepared/alacta_6ip9_v2/strong_site/amide_v3/repair_manifest.json', 55),
}


class RealBoundaryPreparations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.results = {name: prepare(ROOT / path, cls.root / name)
                       for name, (path, _) in FIXTURES.items()}
        cls.connected = prepare(ROOT / FIXTURES['ggr'][0], cls.root / 'connected', 'connected_segment')

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_all_real_pair_counts_charge_coordinates_water_and_parity(self):
        for name, result in self.results.items():
            parent = read_json(ROOT / FIXTURES[name][0])
            self.assertEqual(result['paired_invariants']['atom_count'], FIXTURES[name][1])
            self.assertEqual(result['source_coordinate_max_displacement_A'], 0.)
            self.assertEqual(result['coordination'], parent['coordination'])
            self.assertEqual(result['explicit_water_inventory'], parent['explicit_water_inventory'])
            self.assertEqual(result['charge_ledger'], parent['charge_ledger'])
            self.assertEqual(result['ligand_formal_charge'], -3)
            self.assertEqual(result['decision'], 'uncalibrated_protocol')
            self.assertEqual(result['absolute_reference_status'], 'not_registered')
            la, ca = (xyz(verify(result['outputs'][m]['xyz'])) for m in ('La', 'Ca'))
            self.assertEqual(la[1:], ca[1:])
            self.assertEqual(la[0][1:], ca[0][1:])
            for metal, charge in (('La', 0), ('Ca', -1)):
                self.assertEqual(result['outputs'][metal]['charge'], charge)
                self.assertEqual(result['outputs'][metal]['all_electron_count_for_parity'] % 2, 0)

    def test_ggr_nma_retains_actual_anchors_and_moves_caps_outward(self):
        ggr = self.results['ggr']
        source = {(a['source']['resnum'], a['source']['atom'])
                  for a in ggr['atom_graph']['source_to_qm'] if a['kind'] == 'source'}
        for key in ((140, 'CA'), (140, 'C'), (140, 'O'), (141, 'N'), (141, 'H'), (141, 'CA')):
            self.assertIn(key, source)
        cuts = {(a['retained']['resnum'], a['retained']['atom'], a['omitted']['atom'])
                for a in ggr['atom_graph']['cut_bonds_and_caps']}
        self.assertNotIn((140, 'C', 'CA'), cuts)
        self.assertNotIn((141, 'N', 'CA'), cuts)
        for key in ((140, 'CA', 'N'), (140, 'CA', 'CB'), (141, 'CA', 'C'), (141, 'CA', 'CB')):
            self.assertIn(key, cuts)

    def test_connected_native_segment_merges_donors_without_charge_import(self):
        m = self.connected
        self.assertEqual(m['paired_invariants']['atom_count'], 111)
        self.assertEqual(m['source_coordinate_max_displacement_A'], 0.)
        self.assertEqual(m['ligand_formal_charge'], -3)
        self.assertEqual([r['resnum'] for r in m['extension']['complete_native_residues']], list(range(138, 143)))
        atoms = [(a['source']['resnum'], a['source']['atom'])
                 for a in m['atom_graph']['source_to_qm'] if a['kind'] == 'source']
        self.assertEqual(len(atoms), len(set(atoms)))
        self.assertNotIn((137, 'NZ'), atoms)
        self.assertIn((140, 'NE2'), atoms)
        self.assertIn((141, 'CD1'), atoms)
        self.assertEqual(len(m['atom_graph']['cut_bonds_and_caps']), 7)
        edges = {(a['a']['resnum'], a['a']['atom'], a['b']['resnum'], a['b']['atom'])
                 for a in m['atom_graph']['retained_bonds']}
        for n in range(137, 143):
            self.assertIn((n, 'C', n+1, 'N'), edges)

    def test_repeat_preparation_is_coordinate_and_cache_identical(self):
        repeat = prepare(ROOT / FIXTURES['ggr'][0], self.root / 'repeat')
        first = self.results['ggr']
        self.assertEqual(first['cache_key'], repeat['cache_key'])
        for metal in ('La', 'Ca'):
            self.assertEqual(first['outputs'][metal]['xyz']['sha256'], repeat['outputs'][metal]['xyz']['sha256'])
        with self.assertRaises(InvalidArtifact):
            prepare(ROOT / FIXTURES['ggr'][0], self.root / 'repeat')

    def test_preserved_code_matches_recorded_live_bytes(self):
        for result in (*self.results.values(), self.connected):
            for row in read_json(verify(result['implementation_snapshot'])).values():
                self.assertEqual(record(verify(row['historical_live_source_at_preparation']))['sha256'],
                                 record(verify(row['preserved_copy']))['sha256'])

    def test_missing_anchor_in_corrupted_real_structure_fails(self):
        m = read_json(ROOT / FIXTURES['ggr'][0])
        source = verify(m['source_structure'])
        destination = self.root / 'explicitly_corrupted_missing_Ile141_CA.pdb'
        destination.write_text(''.join(line for line in source.read_text().splitlines(keepends=True)
                                      if not (line.startswith('ATOM  ') and line[21] == 'A'
                                              and line[22:26].strip() == '141'
                                              and line[12:16].strip() == 'CA')))
        graph = SourceGraph(destination, verify(m['topology_definition']))
        selected, _ = _source_selection(graph, m)
        with self.assertRaises(InvalidArtifact):
            expand_selection(graph, selected, m['charge_ledger'], 'alpha_caps')

    def test_no_charge_changing_aequorin_segment_is_silently_used(self):
        with self.assertRaisesRegex(InvalidArtifact, 'outside approved GGR chemistry'):
            prepare(ROOT / FIXTURES['aequorin'][0], self.root / 'unsupported_aequorin_connected', 'connected_segment')
        self.assertFalse((self.root / 'unsupported_aequorin_connected').exists())

    def test_pqq_inputs_remain_pinned_and_are_not_accepted_by_new_policy(self):
        release = ROOT / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json'
        with self.assertRaisesRegex(InvalidArtifact, 'PQQ stays unchanged'):
            prepare(release, self.root / 'unsupported_pqq')
        for row in read_json(release)['scores']:
            for metal in ('La', 'Ca'):
                for kind in ('input', 'xyz'):
                    verify(row['artifacts'][metal][kind])
        self.assertFalse((self.root / 'unsupported_pqq').exists())


if __name__ == '__main__':
    unittest.main()
