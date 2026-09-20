"""Real prepared-fold population and existing physical-map reuse checks."""
from pathlib import Path
import sys
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import accommodation_fold_mappings as fm
from affordable_common import read_json


class FoldMappings(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.preparation = ROOT/'workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json'
        cls.old = ROOT/'workspaces/accommodation_reference_geometry_20260920/mappings_v1/result.json'

    def test_primary_population_and_exact_reuse_overlap(self):
        _, rows, old = fm.population(self.preparation, self.old)
        self.assertEqual(len(rows), 225)
        self.assertEqual(sum(r['status'] == 'prepared' for r in rows), 208)
        self.assertEqual(sum(r['case_id'] in old for r in rows), 57)
        self.assertEqual(sum(r['case_id'] not in old and r['status'] == 'prepared' for r in rows), 151)
        self.assertEqual(sum(r['status'] != 'prepared' for r in rows), 17)
        self.assertTrue(all(not r['source']['canonical_coordinate_match'] for r in rows))

    def test_actual_extra_Asp_and_nonacidic_homolog_mapping_reuse(self):
        _, rows, old = fm.population(self.preparation, self.old)
        for resname in ('ASP', 'ALA', 'SER', 'THR'):
            c = next(r for r in rows if r['case_id'] in old and
                     r['source']['roles']['extra_acidic_ligand_homolog']['resname'] == resname)
            result = fm.supported_case(c, old[c['case_id']])
            self.assertEqual(result['origins'], c['representations']['context']['endpoints'])
            self.assertEqual(result['maps']['Ca'], old[c['case_id']]['mapping'])
            self.assertEqual(result['maps']['Ca'], result['maps']['La'])
            self.assertEqual('extra_acidic_ligand_homolog' in result['modes'], resname == 'ASP')
            self.assertTrue(result['modes']['anchor_glutamate'].endswith('/chi3'))
            if resname == 'ASP': self.assertTrue(result['modes']['extra_acidic_ligand_homolog'].endswith('/chi2'))
            self.assertLessEqual(result['map_origin_max_abs_difference_A'], 1e-12)
            self.assertEqual(result['source_conditioning_metal'], c['source']['source_conditioning_metal'])

    def test_completed_mapping_design_retains_all_sources_and_exact_origins(self):
        _, selected, _ = fm.population(self.preparation, self.old)
        design = read_json(ROOT/'workspaces/accommodation_nonlinear_20260920/fold_maps_v1/design.json')
        source = {r['case_id']: r for r in selected}
        self.assertEqual([r['case_id'] for r in design['rows']], [r['case_id'] for r in selected])
        self.assertEqual(design['status_counts'], {'supported': 208, 'preparation_unavailable': 17})
        self.assertEqual(design['maps_reused'], 57)
        self.assertEqual(design['maps_new_attempted'], 151)
        self.assertEqual(design['new_molecular_calls'], 0)
        for c in design['cases']:
            original = source[c['case_id']]
            self.assertEqual(c['origins'], original['representations']['context']['endpoints'])
            self.assertEqual(c['preparation'], original['representations']['context']['preparation'])
            self.assertEqual(c['source'], original['core'])
            self.assertEqual('extra_acidic_ligand_homolog' in c['modes'], c['extra_homolog_residue'] == 'ASP')
            self.assertFalse(c['canonical_coordinate_match'])
        for r in design['rows']:
            if r['status'] == 'preparation_unavailable':
                self.assertEqual(r['preparation_reason'], source[r['case_id']]['reason'])
                self.assertIsNone(r['mapping_check'])


if __name__ == '__main__': unittest.main()
