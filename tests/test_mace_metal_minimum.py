"""Actual completed curvature grids and fixed native-DFT validation inputs."""
from pathlib import Path
import sys,tempfile,json,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz,InvalidArtifact
from mace_metal_minimum import eligible,validate
from mace_metal_response import prediction
ASSESS=ROOT/'workspaces/mace_metal_response_20260918/assessment_v1/result.json'
PREP=ROOT/'workspaces/mace_metal_response_20260918/minimum_v2/preparation.json'

@unittest.skipUnless(PREP.exists(),'requires actual completed 3D grids and fixed validation preparation')
class MetalMinimumTests(unittest.TestCase):
    def test_actual_gradients_curvature_and_sign_algebra(self):
        a,p,selected=eligible(ASSESS);self.assertEqual(len(selected),4)
        for name,row in a['rows'].items():
            for metal,state in row.items():
                g=np.array(state['DFT_gradient_kcal_mol_A'])+state['J_gradient_kcal_mol_A']
                pred=prediction(g,state['relative_grid_energies_kcal_mol'])
                self.assertEqual(pred,state['prediction'])
                k=sum(np.array(v['fine']) for v in pred['matrices_kcal_mol_A2'].values())
                self.assertGreater(np.linalg.eigvalsh(k).min(),0)
                expected=-.5*g@np.linalg.solve(k,g)
                self.assertAlmostEqual(expected,pred['exploratory_quadratic_change_kcal_mol'],places=10)
            if name.startswith('GGR'):
                expected=row['Ca']['prediction']['predicted_response_kcal_mol']-row['La']['prediction']['predicted_response_kcal_mol']
                self.assertEqual(expected,a['scores'][name]['unvalidated_predicted_delta_R_kcal_mol'])
            else:
                self.assertIsNone(a['scores'][name]['unvalidated_predicted_delta_R_kcal_mol'])

    def test_actual_fixed_quantum_inputs_keep_scaffold_and_analytic_gradients(self):
        a,p,selected=eligible(ASSESS);q=read_json(PREP.parent/'quantum/manifest.json')
        self.assertEqual({t['task_id'] for t in q['tasks']},set(selected));self.assertEqual(len(q['tasks']),4)
        for t in q['tasks']:
            self.assertTrue(t['case_id'].startswith('GGR'))
            text=verify(t['input']).read_text();self.assertIn('r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF EnGrad',text)
            self.assertNotIn('NumGrad',text);self.assertEqual(t['task_type'],'analytic_gradient')
            base=xyz(verify(p['centers'][t['case_id']][t['metal']]['core_state']['xyz']))
            moved=xyz(verify(t['xyz']));self.assertEqual(base[1:],moved[1:])
            self.assertLessEqual(np.linalg.norm(np.array(moved[0][1:])-base[0][1:]),.20)
        self.assertEqual(validate(PREP.parent/'short/manifest.json')['new_short_calls'],8)

    def test_actual_partial_result_keeps_partition_dependent_scores_unavailable(self):
        path=ROOT/'workspaces/mace_metal_response_20260918/minimum_partial_report_v1/result.json'
        if not path.exists():self.skipTest('requires actual partial native DFT report')
        result=read_json(path);self.assertEqual(result['status'],'incomplete')
        self.assertIsNone(result['partition_check']['pass'])
        self.assertIsNotNone(result['scores']['GGR_extended']['endpoint_validated_delta_R_kcal_mol'])
        self.assertTrue(all(row['validated_R_kcal_mol'] is None for row in result['scores'].values()))

    def test_actual_native_DFT_result_and_partition_algebra(self):
        path=ROOT/'workspaces/mace_metal_response_20260918/minimum_report_v2/result.json'
        if not path.exists():self.skipTest('requires completed actual native DFT validation')
        result=read_json(path);verify(result['implementation']);self.assertEqual(result['status'],'complete')
        self.assertEqual(len(result['rows']),4);self.assertTrue(all(row['pass'] for row in result['rows'].values()))
        for name in ('GGR_extended','GGR_connected'):
            expected=result['rows'][name+'_Ca']['actual_energy_change_kcal_mol']-result['rows'][name+'_La']['actual_energy_change_kcal_mol']
            self.assertEqual(result['scores'][name]['validated_delta_R_kcal_mol'],expected)
        delta=result['scores']['GGR_connected']['validated_delta_R_kcal_mol']-result['scores']['GGR_extended']['validated_delta_R_kcal_mol']
        self.assertEqual(result['partition_check']['difference_kcal_mol'],delta)
        self.assertTrue(result['partition_check']['pass'])
        self.assertTrue(all(result['scores'][name]['validated_R_kcal_mol'] is None for name in ('ALPHA_1F6S','ALPHA_6IP9')))

    def test_status_only_promotion_of_real_outside_prediction_rejected(self):
        # Explicitly corrupted copy of a real result; never a scientific output.
        a=read_json(ASSESS);a['rows']['ALPHA_1F6S']['La']['prediction']['status']='eligible_for_DFT_validation'
        with tempfile.TemporaryDirectory() as td:
            damaged=Path(td)/'corrupted_actual_prediction.json';damaged.write_text(json.dumps(a))
            with self.assertRaises(InvalidArtifact):eligible(damaged)

if __name__=='__main__':unittest.main()
