"""Real archived canonical inputs and parser/calibration algebra, without inference."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_canonical_run import inventory,validate
from mace_canonical_report import fit_bands,decision,contrast_components
W=ROOT/'workspaces/mace_canonical_20260916'
I=W/'audit_v2/inventory.json'


@unittest.skipUnless(I.exists(),'requires exact archived canonical source inventory')
class CanonicalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inventory=inventory(I)
        cls.release=read_json(verify(cls.inventory['sources']['calibration']))
        cls.rows=[{'case_id':r['case_id'],'expected_class':r['expected_class'],'evaluation_role':r['evaluation_role'],
                   'status':'computed','R_kcal_mol':r['baseline']['published_R_kcal_mol']} for r in cls.inventory['rows']]

    def test_every_original_case_and_sequence_overlap_is_retained(self):
        self.assertEqual(len(self.inventory['rows']),28)
        self.assertEqual(sum(r['evaluation_role']=='calibration' for r in self.inventory['rows']),25)
        groups={r['case_id']:r['sequence_accession_group'] for r in self.inventory['rows']}
        self.assertEqual(groups['1H4I'],'p16027');self.assertEqual(groups['4MAE'],'i0jwn7')
        self.assertTrue(all(not r['prospectively_blind'] for r in self.inventory['rows']))
        self.assertTrue(all(r['explicit_water_inventory']==[] for r in self.inventory['rows']))

    def test_actual_archive_reproduces_published_band_extrema_and_transfer_calls(self):
        bands=fit_bands(self.rows);published=self.release['calibration']
        self.assertEqual(bands['U_max_Ca_kcal_mol'],published['U_max_Ca_kcal_mol'])
        self.assertEqual(bands['L_min_La_kcal_mol'],published['L_min_La_kcal_mol'])
        for row in self.rows:
            self.assertEqual(decision(row['R_kcal_mol'],bands['decision_bands']),row['expected_class']+'-supported')
        self.assertEqual(bands,fit_bands([r for r in self.rows if r['evaluation_role']=='calibration']))

    def test_missing_real_endpoint_does_not_fit_a_smaller_calibration(self):
        rows=copy.deepcopy(self.rows);rows[0].update(status='unavailable',R_kcal_mol=None)
        bands=fit_bands(rows)
        self.assertEqual(bands['valid_calibration_count'],24)
        self.assertIsNone(bands['decision_bands'])
        self.assertEqual(decision(rows[-1]['R_kcal_mol'],bands['decision_bands']),'unavailable_calibration')
        self.assertEqual(decision(None,bands['decision_bands']),'unavailable_endpoint')

    def test_exact_method_matching_reuses_only_four_medium_crystal_endpoints(self):
        medium=read_json(W/'mace_v2/medium/manifest.json');large=read_json(W/'mace_v2/large/manifest.json')
        self.assertEqual(set(medium['reused']),{'1H4I_La','1H4I_Ca','4MAE_La','4MAE_Ca'})
        self.assertEqual(len(medium['tasks']),52);self.assertEqual(len(large['tasks']),56)
        self.assertEqual(large['reused'],{})
        self.assertEqual(validate(W/'mace_v2/medium/manifest.json')['status'],'pass')
        self.assertEqual(validate(W/'mace_v2/large/manifest.json')['status'],'pass')

    def test_corrupted_real_task_cannot_change_charge_or_calibration_assignment(self):
        m=read_json(W/'mace_v2/medium/manifest.json')
        with tempfile.TemporaryDirectory() as d:
            for field,value in [('charge',m['tasks'][0]['charge']+1),('evaluation_role','retrospective_structural_transfer')]:
                corrupted=copy.deepcopy(m);corrupted['tasks'][0][field]=value
                path=Path(d)/('corrupted_'+field+'.json');path.write_text(json.dumps(corrupted))
                with self.assertRaises(InvalidArtifact):validate(path)

    @unittest.skipUnless((W/'report_v1/result.json').exists(),'requires completed canonical scientific results')
    def test_actual_failed_calibration_cannot_release_transfer_decisions(self):
        result=read_json(W/'report_v1/result.json')
        for model in result['models'].values():
            bands=fit_bands(model['rows'])
            self.assertEqual(bands['status'],'calibration_separation_failed')
            self.assertIsNone(bands['decision_bands'])
            for row in model['rows']:
                if row['evaluation_role']!='calibration':
                    self.assertEqual(decision(row['R_kcal_mol'],bands['decision_bands']),'unavailable_calibration')

    def test_actual_MACE_GB_archive_replays_sign_and_single_unit_conversion(self):
        m=read_json(ROOT/'workspaces/mace_local_correction_20260916/mace_v2/collection_job_1200717.json')
        g=read_json(ROOT/'workspaces/mace_local_correction_20260916/gb_v1/collection_job_1200719.json')
        old=read_json(ROOT/'diagnostics/mace_local_correction_20260916/result.json')
        for name in ('PQQ_1H4I','PQQ_4MAE'):
            endpoint={}
            for metal in ('La','Ca'):
                key=f'{name}_archived_{metal}';mr=m['rows'][key];gr=g['rows'][key]
                endpoint[metal]={'MACE_energy_eV':mr['energy_eV'],'GB_reaction_kcal_mol':gr['GB_reaction_kcal_mol'],
                                 'MACE_components_eV':mr['energy_components_eV']}
            current=contrast_components(endpoint);reference=old['scores'][name]['states']['archived']
            self.assertAlmostEqual(current['R_kcal_mol'],reference['R_kcal_mol'],places=8)
            self.assertAlmostEqual(current['R_GB_kcal_mol'],reference['GB_Ca_minus_La_kcal_mol'],places=10)
            self.assertAlmostEqual(current['R_vacuum_kcal_mol'],reference['MACE_vacuum_R_kcal_mol'],places=8)


if __name__=='__main__':unittest.main()
