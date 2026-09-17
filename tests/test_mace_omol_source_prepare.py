"""Real archived source replay and rejection of unsafe molecular filtering."""
from pathlib import Path
import copy
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, xyz
from mace_omol_source_prepare import inspect_source
from mace_omol_prepared import audit_preparation

WORK = ROOT/'workspaces/mace_omol_20260917'


class RealSourcePreparationTests(unittest.TestCase):
    def setUp(self):
        try:
            import openmm, gemmi
        except ImportError:
            self.skipTest('requires recorded OpenMM/gemmi preparation driver')

    def test_real_ggr_exact_physical_replay(self):
        new = WORK/'ggr_source_bridge_1glg_v2/source_preparation.json'
        old = ROOT/'workspaces/mace_global_benchmark_20260916/prepared_v1/GGR_1GLG/preparation.json'
        if not new.exists() or not old.exists():
            self.skipTest('requires actual source-bridge and original physical preparations')
        self.assertEqual(audit_preparation(new)['status'], 'pass')
        a, b = read_json(new), read_json(old)
        for key in ('physical_atoms', 'preparation_details', 'protein_charge_e', 'cofactor_charge_e',
                    'explicit_waters', 'water_H_moves', 'microstate'):
            self.assertEqual(a[key], b[key])
        for metal in ('La', 'Ca'):
            self.assertEqual(xyz(a['endpoints'][metal]['xyz']['path']), xyz(b['endpoints'][metal]['xyz']['path']))

    def test_real_ggr_cannot_hide_undeclared_sugar_or_change_site(self):
        p = WORK/'ggr_structure_sources_v1/1glg_source_row.json'
        if not p.exists(): self.skipTest('requires real GGR source row')
        row = read_json(p)
        with self.assertRaisesRegex(InvalidArtifact, 'explicit supported exclusion'):
            inspect_source(row, ['HOH'])
        damaged = copy.deepcopy(row); damaged['selected_site']['xyz_A'][0] += 1.
        with self.assertRaisesRegex(InvalidArtifact, 'actual source metal'):
            inspect_source(damaged, ['GAL', 'HOH'])

    def test_real_1kb0_modified_polymer_cannot_be_filtered_to_standard_chain(self):
        p = WORK/'intact_panel_prepared_v2/preparation_manifest.json'
        if not p.exists(): self.skipTest('requires pinned original 1KB0 source records')
        panel = read_json(p); inventory = read_json(panel['source_inventory']['path'])
        source = next(r for r in inventory['rows'] if r['case_id']=='1KB0')
        pm = source['preparation_provenance']['protonation_manifest']
        norm = read_json(source['preparation_provenance']['normalization_manifest']['path'])
        row = {'case_id':'PQQ_1KB0', 'source_structure':read_json(pm['path'])['output'],
               'source_protonation_manifest':pm, 'raw_source':norm['source'],
               'preparation_manifest':source['source_manifest'], 'endpoints':source['endpoints'],
               'selected_site':source['assembly']['selected_site']}
        with self.assertRaisesRegex(InvalidArtifact, 'unsupported source polymer residue: TRO512'):
            inspect_source(row, ['HOH'])


if __name__ == '__main__':
    unittest.main()
