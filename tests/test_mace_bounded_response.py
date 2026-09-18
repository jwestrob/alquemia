"""Bounded algebra and execution invariants from real completed response grids."""
from pathlib import Path
import copy,json,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_bounded_response import ball,predict,sources,validate,endpoint_checks
BASE=ROOT/'workspaces/mace_metal_response_20260918'
ASSESS=BASE/'assessment_v1/result.json'
REUSE=BASE/'minimum_report_v2/result.json'
PREP=ROOT/'workspaces/mace_bounded_response_20260918/prepared_v3/preparation.json'

@unittest.skipUnless(ASSESS.exists() and REUSE.exists(),'requires actual completed response grids and native GGR results')
class BoundedResponseTests(unittest.TestCase):
    def test_real_model_sphere_equations_energy_and_exact_interior_reuse(self):
        a,p,reuse,rp=sources(ASSESS,REUSE)
        for name,row in a['rows'].items():
            for metal,state in row.items():
                b=predict(state);u=np.array(b['displacement_A']);g=np.array(b['gradient_kcal_mol_A'])
                k=sum(np.array(x['fine']) for x in state['prediction']['matrices_kcal_mol_A2'].values())
                self.assertEqual(b['status'],'eligible_for_native_check')
                self.assertLess(np.linalg.norm((k+b['lambda_kcal_mol_A2']*np.eye(3))@u+g),1e-8)
                self.assertAlmostEqual(g@u+.5*u@k@u,b['predicted_change_kcal_mol'],places=10)
                if name.startswith('ALPHA'):
                    self.assertTrue(b['boundary_active']);self.assertAlmostEqual(np.linalg.norm(u),.20,places=10)
                    self.assertGreater(b['lambda_kcal_mol_A2'],0)
                    # The unconstrained shorthand is incorrect at an active sphere.
                    self.assertGreater(abs(b['predicted_change_kcal_mol']+.5*g@np.linalg.solve(k,g)),.1)
                else:
                    self.assertFalse(b['boundary_active'])
                    self.assertEqual(b['displacement_A'],rp['states'][name+'_'+metal]['prediction']['proposed_displacement_A'])
                    old=reuse['rows'][name+'_'+metal]
                    check=endpoint_checks(b,p['centers'][name][metal],old['DFT_change_kcal_mol'],old['J_change_kcal_mol'],old['actual_gradient_kcal_mol_A'])
                    self.assertTrue(check['pass'])
                    self.assertEqual(check['actual_energy_change_kcal_mol'],old['actual_energy_change_kcal_mol'])

    def test_corrupted_real_matrix_and_prediction_are_rejected(self):
        state=read_json(ASSESS)['rows']['ALPHA_1F6S']['La']
        k=sum(np.array(x['fine']) for x in state['prediction']['matrices_kcal_mol_A2'].values())
        g=np.array(state['DFT_gradient_kcal_mol_A'])+state['J_gradient_kcal_mol_A']
        # Explicitly malformed real-derived matrix, not scientific data.
        damaged=k.copy();damaged[0,1]+=1
        with self.assertRaises(InvalidArtifact):ball(damaged,g)
        self.assertEqual(ball(-k,g)['status'],'unstable_curvature')
        damaged=copy.deepcopy(state);damaged['prediction']['proposed_displacement_A'][0]+=.01
        with self.assertRaises(InvalidArtifact):predict(damaged)

    @unittest.skipUnless(PREP.exists(),'requires the bounded preparation')
    def test_fixed_native_inputs_and_only_selected_metal_moves(self):
        a,p,reuse,rp=sources(ASSESS,REUSE);top=read_json(PREP)
        self.assertEqual(set(top['reused']),set(reuse['rows']));self.assertEqual(len(top['states']),4)
        self.assertEqual(top['excluded'],{})
        self.assertEqual(validate(PREP.parent/'short/manifest.json')['new_short_calls'],8)
        qm=read_json(PREP.parent/'quantum/manifest.json');self.assertEqual(len(qm['tasks']),4)
        for t in qm['tasks']:
            self.assertTrue(t['case_id'].startswith('ALPHA'))
            native=verify(t['input']).read_text()
            self.assertIn('r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF EnGrad',native)
            self.assertNotIn('NumGrad',native)
            old=xyz(verify(p['centers'][t['case_id']][t['metal']]['core_state']['xyz']))
            moved=xyz(verify(t['xyz']));self.assertEqual(old[1:],moved[1:])
            self.assertAlmostEqual(np.linalg.norm(np.array(moved[0][1:])-old[0][1:]),.20,places=8)
            self.assertEqual(t['charge'],p['centers'][t['case_id']][t['metal']]['core_state']['charge'])
        # Explicitly corrupted manifest keeps its original key: must not hit cache.
        manifest=read_json(PREP.parent/'short/manifest.json');manifest['tasks'][0]['charge']+=1
        with tempfile.TemporaryDirectory() as td:
            damaged=Path(td)/'corrupted_real_manifest.json';damaged.write_text(json.dumps(manifest))
            with self.assertRaises(InvalidArtifact):validate(damaged)

    def test_actual_native_result_signs_gradients_and_margins(self):
        path=ROOT/'workspaces/mace_bounded_response_20260918/report_v1/result.json'
        if not path.exists():self.skipTest('requires completed real bounded native validation')
        from affordable_common import HA_TO_KCAL,BOHR_TO_A
        from affordable_response import read_engrad
        from mace_hybrid import EV_TO_KCAL
        result=read_json(path);verify(result['implementation'])
        a,p,reuse,rp=sources(ASSESS,REUSE)
        dc=read_json(verify(result['sources']['DFT']));sc=read_json(verify(result['sources']['short']))
        self.assertEqual(len(result['rows']),8);self.assertTrue(result['partition_check']['pass'])
        self.assertTrue(all(row['pass'] for row in result['rows'].values()))
        for key,computed in dc['rows'].items():
            name,metal=key.rsplit('_',1);center=p['centers'][name][metal]
            dft=(computed['energy_hartree']-center['DFT_energy_Ha'])*HA_TO_KCAL
            full,core=sc['rows']['full__'+key],sc['rows']['core__'+key]
            j=((full['energy_eV']-center['full_short']['energy_eV'])-(core['energy_eV']-center['core_short']['energy_eV']))*EV_TO_KCAL
            self.assertEqual(result['rows'][key]['actual_energy_change_kcal_mol'],dft+j)
            raw=read_engrad(verify(dc['gradients'][key]['artifacts']['engrad']))
            gd=np.array(raw['gradient_Ha_per_bohr'][0])*HA_TO_KCAL/BOHR_TO_A
            gj=(-np.load(verify(full['forces']))[center['full_metal_index']]+np.load(verify(core['forces']))[0])*EV_TO_KCAL
            np.testing.assert_allclose(result['rows'][key]['actual_gradient_kcal_mol_A'],gd+gj,atol=1e-10,rtol=0)
        for name,row in result['scores'].items():
            expected=result['rows'][name+'_Ca']['actual_energy_change_kcal_mol']-result['rows'][name+'_La']['actual_energy_change_kcal_mol']
            self.assertEqual(row['qualified_delta_R_kcal_mol'],expected)
            self.assertIsNone(row['calibrated_class']);self.assertIsNone(row['entropy_correction_kcal_mol'])
        for name,row in result['comparisons'].items():
            self.assertFalse(row['qualified_direction_pass'])
            self.assertGreater(row['actual_margin_change_kcal_mol'],0)
            expected=result['scores'][name]['qualified_R_kcal_mol']-result['scores']['GGR_extended']['qualified_R_kcal_mol']
            self.assertEqual(row['qualified_R_kcal_mol'],expected)

if __name__=='__main__':unittest.main()
