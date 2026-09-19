"""Real archived structure/energy checks; no fabricated scientific outputs."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import read_json, verify, xyz, InvalidArtifact
from hydration_scanner import transfer, contrasts
from mace_hybrid import EV_TO_KCAL


class ScannerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ip = ROOT/'workspaces/site_classifier_20260918/inventory_v2/inventory.json'
        cp = ROOT/'workspaces/hydration_network_20260918/core_transfer_v1/manifest.json'
        if not ip.exists() or not cp.exists():
            raise unittest.SkipTest('real pinned alpha/GGR fixtures unavailable')
        cls.inventory = read_json(ip); cls.core = read_json(cp)

    def test_source_transfer_preserves_scaffold_and_water_shape(self):
        for t in self.core['tasks']:
            row = next(r for r in self.inventory['rows'] if r['case_id'] == 'ALPHA_'+t['case'])
            p = read_json(verify(row['preparation'])); parent = read_json(verify(t['parent']))
            original = xyz(verify(p['endpoints'][t['metal']]['xyz'])); core = xyz(verify(t['xyz']))
            result, mapping = transfer(p, parent, core, t['metal'])
            moving = {a['whole_index'] for a in mapping}
            self.assertEqual(len(moving), 2*len(p['explicit_waters']))
            for i in range(len(result)):
                if i not in moving: self.assertEqual(result[i], original[i])
            for a in mapping: self.assertEqual(result[a['whole_index']], core[a['core_index']])

    def test_corrupted_real_water_mapping_rejected(self):
        t = self.core['tasks'][0]
        r = next(r for r in self.inventory['rows'] if r['case_id'] == 'ALPHA_'+t['case'])
        p = read_json(verify(r['preparation'])); parent = read_json(verify(t['parent']))
        corrupted = copy.deepcopy(parent)
        corrupted['atom_graph']['source_to_qm'].pop()
        with self.assertRaises(InvalidArtifact): transfer(p, corrupted, xyz(verify(t['xyz'])), t['metal'])
        coords = xyz(verify(t['xyz'])); oi = next(a['source_qm_index'] for a in p['physical_atoms'] if a['kind'] == 'retained_site_water' and a['element'] == 'O')
        changed = list(coords); o = coords[oi]; changed[oi] = (o[0], o[1]+.1, o[2], o[3])
        with self.assertRaises(InvalidArtifact): transfer(p, parent, changed, t['metal'])

    def test_archived_energy_algebra_and_missing_reference(self):
        d = read_json(ROOT/'workspaces/mace_omol_20260917/charge_ablation_report_v1/result.json')
        for name in ('ALPHA_1F6S', 'ALPHA_6IP9', 'GGR_1GLG'):
            row = d['scores'][name]
            inv = next(r for r in self.inventory['rows'] if r['case_id'] == name)
            energies = {(m, pos): row['endpoints'][m][pos]['energy_eV'] for m in ('Ca', 'La') for pos in ('bound', 'detached')}
            c = contrasts(energies, inv['disconnected'])
            self.assertAlmostEqual(c['interaction_R_model_kcal'], row['R_mask_model_kcal'], places=7)
            self.assertLess(abs(c['component_closure_error_model_kcal']), 1e-7)
            bound = {k: v for k, v in energies.items() if k[1] == 'bound'}
            with self.assertRaises(InvalidArtifact): contrasts(bound, inv['disconnected'])
            total = contrasts(bound, inv['disconnected'], common_geometry=True)
            self.assertAlmostEqual(total['interaction_R_model_kcal'], c['interaction_R_model_kcal'], places=6)

    def test_actual_scientific_integration(self):
        path = ROOT/'workspaces/hydration_scanner_20260919/mace_v2/collection_job_1202084.json'
        if not path.exists():
            self.skipTest('actual allocated MACE integration has not been collected')
        r = read_json(path)
        from hydration_scanner import collect
        self.assertEqual(r, collect(verify(r['manifest'])))
        self.assertEqual(r['status'], 'complete')
        self.assertEqual(len(r['rows']), 12)
        self.assertEqual(len(r['dry_identities']), 30)
        self.assertEqual(sum(c['MACE_normalized'] > 0 for c in r['comparisons']), 2)
        self.assertEqual(sum(c['MACE_contextual'] > 0 for c in r['comparisons']), 6)
        for c in r['cases']:
            self.assertLess(abs(c['MACE_contextual']['component_closure_error_model_kcal']), 1e-7)


if __name__ == '__main__': unittest.main()
