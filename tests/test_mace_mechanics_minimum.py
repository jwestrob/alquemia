"""Selection and geometry integrity on the actual, frozen response predictions."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_mechanics_minimum import eligible,validate,ELIGIBLE
A=ROOT/'workspaces/mace_mechanics_20260916/assessment_v1/result.json'
P=ROOT/'workspaces/mace_mechanics_20260916/minimum_v1/preparation.json'


@unittest.skipUnless(A.exists() and P.exists(),'requires actual completed mechanics assessment and fixed minima')
class MinimumTests(unittest.TestCase):
    def test_all_and_only_eligible_predictions_are_prepared(self):
        a,selected=eligible(A);p=read_json(P)
        self.assertEqual(set(selected),{'ALPHA_1F6S_Ca','GGR_connected_Ca','GGR_extended_Ca'})
        self.assertEqual(set(selected),set(p['minima']))
        self.assertEqual(p['new_DFT_calls'],3);self.assertEqual(p['new_short_calls'],6)
        for key,row in p['minima'].items():
            self.assertEqual(row['prediction'],selected[key]['prediction'])
            self.assertLessEqual(row['physical_validation']['maximum_heavy_displacement_A'],.05)
        self.assertTrue(all(a['rows'][c]['stationary_points']['La']['predicted_relaxation_kcal_mol'] is None for c in a['rows']))

    def test_ineligible_real_prediction_cannot_be_promoted_by_a_status_flag(self):
        a=copy.deepcopy(read_json(A))
        a['rows']['ALPHA_1F6S']['stationary_points']['La']['status']=ELIGIBLE
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'corrupted_actual_assessment.json';f.write_text(json.dumps(a))
            with self.assertRaises(InvalidArtifact):eligible(f)

    def test_qm_caps_remain_fixed_and_source_coordinates_match_full_system(self):
        p=read_json(P);prep=read_json(verify(p['original_preparation']))
        for key,row in p['minima'].items():
            case=read_json(verify(prep['cases'][row['case_id']]))
            repair=read_json(verify(case['source_preparation']))
            core=xyz(verify(row['endpoints']['core']['xyz']));full=xyz(verify(row['endpoints']['full']['xyz']))
            center=xyz(verify(case['grids']['center']['endpoints'][row['metal']]['xyz']))
            physical=read_json(verify(case['physical_atoms']))
            by_id={atom['source_key']:full[i][1:] for i,atom in enumerate(physical)}
            from affordable_response import source_key
            for atom in repair['atom_graph']['source_to_qm']:
                i=atom['qm_index']
                if atom['kind']=='cap':self.assertEqual(core[i],center[i])
                elif atom['kind']=='source':
                    np.testing.assert_allclose(core[i][1:],by_id[source_key(atom['source'])],atol=1e-9,rtol=0)
            self.assertEqual(core[0][0],row['metal'])

    def test_successful_Ca_validation_cannot_fill_missing_La_or_release_a_score(self):
        path=ROOT/'workspaces/mace_mechanics_20260916/minimum_validation_v1/result.json'
        if not path.exists():self.skipTest('requires actual independently executed DFT validation')
        result=read_json(path)
        self.assertEqual(len(result['rows']),3)
        self.assertTrue(all(row['pass'] for row in result['rows'].values()))
        for row in result['scores'].values():
            self.assertIsNone(row['endpoint_response_kcal_mol']['La'])
            self.assertIsNone(row['response_correction_kcal_mol'])
            self.assertIsNone(row['corrected_R_kcal_mol'])
            self.assertIsNone(row['calibrated_class'])
        self.assertEqual(result['partition_check']['status'],'unavailable')
        self.assertFalse(result['predictive_improvement_claimed'])


if __name__=='__main__':unittest.main()
