"""Actual source receipts and algebra; no invented GFN2 scientific output."""
from pathlib import Path
import sys
import json
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, verify
from compact_solvation_compare import comparison, mix_pair, native_endpoint


class CompactSources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'diagnostics/compact_solvation_20260920/INVENTORY.json'
        cls.inv=read_json(cls.path)

    def test_all_132_actual_native_sources_and_fixed_denominators(self):
        self.assertEqual(len(self.inv['cases']),33)
        self.assertEqual(sum(c['role']=='calibration' for c in self.inv['cases']),25)
        self.assertEqual(sum(c['role'] not in ('calibration','direct_direction_development') for c in self.inv['cases']),3)
        n=0
        for case in self.inv['cases']:
            for rep in case['representations'].values():
                for source in rep['endpoints'].values():
                    actual=read_json(verify(source['native_MACE_receipt']))
                    self.assertEqual(native_endpoint(actual,self.inv['model']),source);n+=1
        self.assertEqual(n,132)
        self.assertEqual(self.inv['comparison_definitions']['alpha_GGR_biological_comparisons'],1)

    def test_replays_existing_alpha_ggr_native_context_gain(self):
        cases={c['case_id']:c for c in self.inv['cases']}
        counts={}
        for representation in ('core','context'):
            margins=[cases[p['La_like']]['representations'][representation]['native_R_model_kcal_mol']-
                     cases[p['Ca_like']]['representations'][representation]['native_R_model_kcal_mol']
                     for p in self.inv['comparison_definitions']['alpha_GGR_pairs']]
            counts[representation]=sum(v>0 for v in margins)
        self.assertEqual(counts,{'core':2,'context':6})

    def test_no_low_level_outputs_retains_native_and_null_composite(self):
        result=comparison(self.path,[],'native')
        self.assertEqual(result['case_representation_denominator'],66)
        self.assertEqual(result['available_composite_case_representations'],0)
        self.assertTrue(all(r['native_status']=='available' and r['composite_R_model_kcal_mol'] is None for r in result['rows']))
        self.assertTrue(all(v['bands'] is None for v in result['calibration'].values()))
        self.assertTrue(all(p['composite_margin_model_kcal_mol'] is None for p in result['alpha_GGR_matrix']))

    def test_identity_unit_algebra_uses_actual_gfn_values(self):
        case=next(c for c in self.inv['cases'] if c['case_id']=='1H4I')
        mace={m:case['representations']['context']['endpoints'][m]['native_MACE_energy_eV'] for m in ('Ca','La')}
        quantum=read_json(ROOT/'workspaces/compact_solvation_20260920/pilot_v1/collection_1203145.json')
        ha={m:next(r['vacuum_hartree'] for r in quantum['rows'] if r['case_id']=='1H4I' and r['metal']==m) for m in ('Ca','La')}
        result=mix_pair(mace,ha,ha)
        self.assertEqual(result['solvation_delta_R_kcal_mol'],0)
        self.assertEqual(result['composite_R_model_kcal_mol'],case['representations']['context']['native_R_model_kcal_mol'])
        self.assertEqual(result['GFN2_vacuum_R_kcal_mol'],(ha['Ca']-ha['La'])*HA_TO_KCAL)
        with self.assertRaises(InvalidArtifact):mix_pair(mace,ha,{'Ca':None,'La':ha['La']})

    def test_explicit_corruption_of_actual_native_energy_is_rejected(self):
        source=self.inv['cases'][0]['representations']['core']['endpoints']['Ca']
        actual=read_json(verify(source['native_MACE_receipt']))
        actual['energy_eV']+=1 # identified corrupted real receipt copy, never science evidence
        with self.assertRaises(InvalidArtifact):native_endpoint(actual,self.inv['model'])

    def test_actual_eight_call_gfn_pilot_partial_coverage_and_algebra(self):
        path=ROOT/'workspaces/compact_solvation_20260920/pilot_v1/collection_1203145.json'
        result=comparison(self.path,[path],'native');actual=read_json(path)
        self.assertEqual(result['available_composite_case_representations'],2)
        self.assertEqual(sum(r['composite_R_model_kcal_mol'] is None for r in result['rows']),64)
        self.assertIsNone(result['calibration']['context']['bands'])
        for case in ('1H4I','q9z4j7-pqq-la_model'):
            eps={m:next(r for r in actual['rows'] if r['case_id']==case and r['metal']==m) for m in ('Ca','La')}
            row=next(r for r in result['rows'] if r['case_id']==case and r['representation']=='context')
            expected=((eps['Ca']['alpb_hartree']-eps['Ca']['vacuum_hartree'])-(eps['La']['alpb_hartree']-eps['La']['vacuum_hartree']))*HA_TO_KCAL
            self.assertEqual(row['solvation_delta_R_kcal_mol'],expected)
            self.assertEqual(row['composite_R_model_kcal_mol'],row['native_R_model_kcal_mol']+expected)
            src=next(c for c in self.inv['cases'] if c['case_id']==case)['representations']['context']
            energies={m:src['endpoints'][m]['native_MACE_energy_eV'] for m in ('Ca','La')}
            vacuum={m:eps[m]['vacuum_hartree'] for m in ('Ca','La')}
            alpb={m:eps[m]['alpb_hartree'] for m in ('Ca','La')}
            # Direct algebra on actual GFN values; reversed transfer reverses sign.
            self.assertEqual(mix_pair(energies,alpb,vacuum)['solvation_delta_R_kcal_mol'],-expected)
            self.assertEqual(mix_pair(energies,vacuum,vacuum)['solvation_delta_R_kcal_mol'],0)

    def test_corrupted_real_gfn_collection_is_rejected(self):
        actual=read_json(ROOT/'workspaces/compact_solvation_20260920/pilot_v1/collection_1203145.json')
        actual['rows'][0]['vacuum_hartree']+=0.001 # explicitly corrupted real fixture copy
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'corrupted_real_collection.json';path.write_text(json.dumps(actual))
            with self.assertRaises(InvalidArtifact):comparison(self.path,[path],'native')

    def test_actual_audited_primary_and_failed_numerical_check_stay_separate(self):
        primary=ROOT/'workspaces/compact_solvation_20260920/pilot_v1/audited_collection_v2.json'
        failed=ROOT/'workspaces/compact_solvation_20260920/numerical_v1/failed_collection_v1.json'
        result=comparison(self.path,[primary],'native',[failed])
        self.assertEqual(result['available_composite_case_representations'],2)
        self.assertEqual(result['numerical_qualification'],'numerical_crosscheck_unavailable')
        self.assertEqual(len(result['numerical_check']['rows']),4)
        self.assertIsNone(result['numerical_check']['score_correction_differences_kcal_mol'])
        self.assertTrue(all(e['energy_hartree'] is None for row in result['numerical_check']['rows'] for e in row['endpoints'].values()))
        self.assertTrue(all(r['composite_R_model_kcal_mol'] is None for r in result['rows'] if r['composite_status']=='unavailable'))
        with self.assertRaises(InvalidArtifact):comparison(self.path,[failed],'ordinary_tight',[failed])


if __name__=='__main__':unittest.main()
