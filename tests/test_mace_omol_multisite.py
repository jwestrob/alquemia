"""Pinned real multisite fixtures; no generated structures or model outputs."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, xyz
from mace_hybrid import check_atoms
from mace_omol_multisite import audit


class RealMultisitePreparation(unittest.TestCase):
    def setUp(self):
        self.root = ROOT/'workspaces/mace_omol_20260917/multisite_prepared_v2'
        if not (self.root/'PARV_4CPV_CD/preparation.json').exists():
            self.skipTest('requires actual multisite source preparations and OpenMM driver')

    def test_real_acetyl_background_ions_waters_and_common_ca_geometry(self):
        for prefix, sites, waters, backgrounds in [('PARV_4CPV', ['CD','EF'],1,1),
                                                   ('AEQ_1SL8',['EF1','EF3','EF4'],3,2)]:
            ca_systems = []
            for site in sites:
                path = self.root/(prefix+'_'+site)/'preparation.json'; p = read_json(path)
                self.assertEqual(audit(path)['status'], 'pass')
                self.assertEqual(len(p['background_metals']), backgrounds)
                self.assertEqual(len(p['explicit_waters']), waters)
                self.assertEqual(p['endpoints']['La']['charge']-p['endpoints']['Ca']['charge'],1)
                ca_systems.append(sorted(xyz(p['endpoints']['Ca']['xyz']['path'])))
                if prefix.startswith('PARV'):
                    self.assertEqual(sum(a.get('resname')=='ACE' for a in p['physical_atoms']),6)
                    self.assertIn(['A/0//C','A/1//N'],p['source_covalent_connections'])
            self.assertTrue(all(s==ca_systems[0] for s in ca_systems))

    def test_real_multisite_state_requires_explicit_background_mapping(self):
        p = read_json(self.root/'AEQ_1SL8_EF1/preparation.json')
        e = p['endpoints']['La']; coords = xyz(e['xyz']['path'])
        with self.assertRaisesRegex(InvalidArtifact,'exactly one selected metal'):
            check_atoms(coords,e['charge'])
        indices = tuple(i for i,a in enumerate(p['physical_atoms']) if a['kind']=='background_metal')
        self.assertEqual(check_atoms(coords,e['charge'],background_calcium_indices=indices),e['state'])
        with self.assertRaisesRegex(InvalidArtifact,'background calcium'):
            check_atoms(coords,e['charge'],background_calcium_indices=indices+(len(coords)-1,))

    def test_corrupted_real_preparation_cannot_drop_background_ion(self):
        p = copy.deepcopy(read_json(self.root/'PARV_4CPV_CD/preparation.json'))
        p['physical_atoms'] = [a for a in p['physical_atoms'] if a['kind']!='background_metal']
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'corrupted_real_preparation.json'; path.write_text(json.dumps(p))
            with self.assertRaisesRegex(InvalidArtifact,'full source replay'): audit(path)

    def test_actual_multisite_scores_keep_failed_supporting_gate_and_unresolved_labels(self):
        from mace_omol_multisite_report import report, ORDER
        w = self.root.parent
        if not (w/'multisite_comparison_v1/result.json').exists():
            self.skipTest('requires actual completed ten-forward multisite experiment')
        with tempfile.TemporaryDirectory() as tmp:
            r = report([w/'multisite_reports_v1'/n/'result.json' for n in ORDER],
                       w/'ggr_structure_report_v1/result.json',
                       ROOT/'diagnostics/mace_omol_20260917/MULTISITE_PANEL_PLAN.md',Path(tmp)/'report')
        self.assertTrue(r['numerical_gate_pass'])
        self.assertFalse(r['parvalbumin_supporting_all_case_gate_pass'])
        self.assertEqual(sum(x['pass'] for x in r['parvalbumin_supporting_contrasts']),4)
        self.assertEqual(len(r['aequorin_ordered_vector_model_kcal']),3)
        self.assertIsNone(r['aequorin_site_resolved_direction'])
        self.assertEqual(r['gold_same_assay_site_resolved_groups_added'],0)
        self.assertEqual(r['model_cost']['new_forwards'],10)


if __name__=='__main__': unittest.main()
