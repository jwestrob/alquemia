"""Actual matched-H response grids, gradients, and independently prepared points."""
from pathlib import Path
import copy,json,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz,InvalidArtifact
from mace_hybrid import EV_TO_KCAL
from mace_omol_hybrid_minimum import inputs,validate,validate_quantum,report
from mace_metal_response import matrix,grid
BASE=ROOT/'workspaces/mace_omol_hybrid_response_20260918'
ASSESS=BASE/'assessment_v1/result.json'
PREP=BASE/'minimum_v1/preparation.json'

@unittest.skipUnless(ASSESS.exists(),'requires all actual completed matched-H response grids')
class HybridMinimumTests(unittest.TestCase):
    def test_same_hybrid_gradient_and_whole_curvature_from_real_outputs(self):
        a,initial,p,r,gm,states,selected=inputs(ASSESS);self.assertEqual(len(selected),8)
        rows={}
        for ref in a['collections']:rows.update(read_json(verify(ref))['rows'])
        for key,item in initial['reused_GGR_whole_gradients'].items():rows['full__'+key+'__center']=item['result']
        for name,ends in a['rows'].items():
            for metal,pred in ends.items():
                s=states[name][metal];gid=s['global_id'];full=rows['full__'+gid+'_'+metal+'__center'];core=rows['core__'+name+'_'+metal+'__center']
                context=(np.load(verify(full['gradient']))[s['full']['metal_index']]-np.load(verify(core['gradient']))[0])*EV_TO_KCAL
                np.testing.assert_allclose(pred['gradient_kcal_mol_A'],np.array(s['DFT_gradient_kcal_mol_A'])+context,atol=1e-10,rtol=0)
                values={point:(rows['full__'+gid+'_'+metal+'__'+point]['energy_eV']-full['energy_eV'])*EV_TO_KCAL for point in grid()}
                for scale in ('coarse','fine'):np.testing.assert_allclose(matrix(values,scale),pred['matrices_kcal_mol_A2'][scale],atol=1e-10,rtol=0)
        self.assertEqual(a['rows']['GGR_extended']['La']['matrices_kcal_mol_A2'],a['rows']['GGR_connected']['La']['matrices_kcal_mol_A2'])

    @unittest.skipUnless(PREP.exists(),'requires actual fixed native validation preparation')
    def test_fixed_validation_coordinates_method_and_charge(self):
        a,initial,p,r,gm,states,selected=inputs(ASSESS);top=read_json(PREP)
        self.assertEqual(validate(PREP.parent/'mace/manifest.json')['tasks'],2*len(top['states']))
        self.assertEqual(validate_quantum(PREP.parent/'quantum/manifest.json')['status'],'dry_run_pass')
        qm=read_json(PREP.parent/'quantum/manifest.json');self.assertEqual(len(qm['tasks']),len(top['states']))
        for t in qm['tasks']:
            source=states[t['case_id']][t['metal']]['core'];old=xyz(verify(source['xyz']));new=xyz(verify(t['xyz']))
            self.assertEqual(old[1:],new[1:]);self.assertEqual(t['charge'],source['charge'])
            np.testing.assert_allclose(np.array(new[0][1:])-old[0][1:],selected[t['task_id']]['prediction']['displacement_A'],atol=5e-10,rtol=0)
            text=verify(t['input']).read_text();self.assertIn('r2SCAN-3c NoAutostart DefGrid3 TightSCF EnGrad',text)
            self.assertNotIn('CPCM',text);self.assertNotIn('NumGrad',text)
        self.assertEqual(set(top['states'])|set(top['excluded']),set(selected))
        self.assertEqual(top['unavailable_predictions'],{})

    def test_corrupted_real_prediction_not_admitted(self):
        a=read_json(ASSESS);a['rows']['ALPHA_1F6S']['La']['prediction']['displacement_A'][0]+=.01
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted_actual_prediction.json';path.write_text(json.dumps(a))
            with self.assertRaises(InvalidArtifact):inputs(path)

    @unittest.skipUnless((BASE/'native_report_v2/result.json').exists(),'requires actual native DFT and MACE endpoints')
    def test_actual_ordering_and_failed_qualification_remain_separate(self):
        archived=read_json(BASE/'native_report_v2/result.json')
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'replayed_report'
            report(PREP,PREP.parent/'quantum/manifest.json',PREP.parent/'mace/manifest.json',out)
            actual=read_json(out/'result.json')
        for key in ('rows','scores','partition','contrasts'):
            self.assertEqual(actual[key],archived[key])
        self.assertTrue(actual['all_four_actual_direction_gate'])
        self.assertFalse(actual['all_four_direction_gate'])
        failures={key:[name for name,passed in row['checks'].items() if not passed]
                  for key,row in actual['rows'].items() if not row['pass']}
        self.assertEqual(failures,{'ALPHA_1F6S_Ca':['boundary_radial_sign'],
                                   'GGR_extended_Ca':['boundary_radial_sign']})
        self.assertTrue(actual['partition']['actual_response_threshold_pass'])
        self.assertFalse(actual['partition']['actual_final_threshold_pass'])
        self.assertAlmostEqual(actual['partition']['actual_final_difference_kcal'],2.023524103802629,places=8)
        for name,row in actual['scores'].items():
            ends=[actual['rows'][name+'_'+metal] for metal in ('Ca','La')]
            delta=ends[0]['actual_energy_change_kcal_mol']-ends[1]['actual_energy_change_kcal_mol']
            self.assertEqual(row['actual_delta_R_kcal_scale'],delta)
            self.assertIsNone(row['qualified_R_kcal_scale'])
            self.assertIsNone(row['calibrated_class'])
        for row in actual['contrasts']:
            self.assertGreater(row['actual_margin_kcal'],8.7)
            self.assertIsNone(row['qualified_margin_kcal'])

if __name__=='__main__':unittest.main()
