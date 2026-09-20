"""Real-state/mapping and archived-gradient checks; no simulated solver success."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz,paired
from mace_site_kinematics import Kinematics
import accommodation_nonlinear as nl


class Nonlinear(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=read_json(ROOT/'workspaces/accommodation_torsion_20260920/prepared_v3/design.json')
        cls.primary=read_json(ROOT/'workspaces/accommodation_torsion_20260920/primary_result_v1.json')
        cls.response=read_json(ROOT/'workspaces/accommodation_response_20260920/result_v2.json')
        cls.ad=read_json(ROOT/'workspaces/accommodation_response_20260920/prepared_v2/design.json')

    def test_actual_eight_endpoint_coordinates_charge_and_moving_maps_match(self):
        self.assertEqual(tuple(c['case_id'] for c in self.source['cases']),nl.CASES)
        for c in self.source['cases']:
            ca,la=[c['origins'][z] for z in ('Ca','La')]
            paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
            self.assertEqual(read_json(verify(c['maps']['Ca'])),read_json(verify(c['maps']['La'])))
            kin=Kinematics(read_json(verify(c['maps']['Ca']))['context']);names=[m['id'] for m in kin.modes]
            active=[names.index(mid) for mid in c['modes'].values()]
            self.assertTrue(all(kin.modes[i]['unit']=='radian' for i in active))
            self.assertTrue(all(not names[i].startswith('metal_') for i in active))
            q=np.zeros(len(names));q[active]=.01
            self.assertTrue(kin.check(q)['pass'])
            # Source nuclei outside the active terminal groups remain fixed.
            moving={a for i in active for a in kin.modes[i]['moving_indices']}
            fixed=[i for i in range(len(kin.positions)) if i not in moving]
            # Existing rotation algebra reconstructs inactive bonds at angle0;
            # observed ~2e-16 A arithmetic roundoff is below its existing map
            # validation tolerance. No source coordinates are edited here.
            np.testing.assert_allclose(kin.evaluate(q)[0][fixed],kin.positions[fixed],atol=1e-12,rtol=0)

    def test_actual_analytic_components_sign_units_and_current_jacobian(self):
        collections=[]
        for stage,job in [('gate',1203465),('followup',1203499)]:
            collections+=read_json(ROOT/f'workspaces/accommodation_response_20260920/prepared_v2/{stage}/collection_{job}.json')['rows']
        expected=next(c for c in self.response['cases'] if c['case_id']=='1H4I')
        for z in ('Ca','La'):
            source=self.ad['sources']['1H4I__'+z];kin=Kinematics(read_json(verify(self.ad['maps']['1H4I__'+z]))['context'])
            forces=np.load(verify(source['forces']),allow_pickle=False)
            low={medium:next(r for r in collections if r['case_id']=='1H4I' and r['metal']==z and r['medium']==medium and r['kind']=='center') for medium in ('vacuum','alpb')}
            jac=kin.evaluate(np.zeros(len(kin.modes)))[3]
            ng,sg,cg=nl.project_components(jac,forces,low['vacuum']['gradient_kcal_mol_A'],low['alpb']['gradient_kcal_mol_A'])
            for actual,name in ((ng,'vacuum_MACE'),(sg,'solvent_transfer'),(cg,'composite')):
                np.testing.assert_allclose(actual,expected['endpoints'][z][name],atol=1e-10,rtol=0)
            i=next(i for i,m in enumerate(kin.modes) if m['id']=='A/177/chi3')
            q=np.zeros(len(kin.modes));q[i]=.1
            # Different J(q), with exactly the same actual Cartesian force array,
            # confirms that a stale origin projection is not baked into helper.
            fresh=nl.project_components(kin.evaluate(q)[3],forces,low['vacuum']['gradient_kcal_mol_A'],low['alpb']['gradient_kcal_mol_A'])[2]
            self.assertGreater(abs(fresh[i]-cg[i]),.01)

    def test_real_archived_energy_work_and_score_algebra(self):
        points={p['task_id']:p for p in self.primary['points']}
        for profile in self.primary['profiles']:
            for row in profile['rows']:
                for z,endpoint in row['endpoints'].items():
                    if endpoint.get('work_kcal_mol') is None:continue
                    point=points[endpoint['task_id']];origin=points[profile['case_id']+'__origin__'+z]
                    if point['status']!='complete':
                        self.assertIsNone(endpoint['work_kcal_mol'].get('composite'))
                        continue  # actual archived solver failure remains missing
                    work=nl.relative_components(point,origin);expected=endpoint['work_kcal_mol']
                    self.assertAlmostEqual(work['native_MACE_kcal_mol'],expected['MACE'],places=6)
                    self.assertAlmostEqual(work['solvent_transfer_kcal_mol'],expected['solvent_transfer'],places=6)
                    self.assertAlmostEqual(work['composite_kcal_mol'],expected['composite'],places=6)
        self.assertEqual(nl.relative_components(points['1H4I__origin__Ca'],points['1H4I__origin__Ca'])['composite_kcal_mol'],0.)

    def test_actual_nonstationary_gradient_and_unrefined_curvature_rejected(self):
        case=next(c for c in self.response['cases'] if c['case_id']=='1H4I')
        i=case['mode_ids'].index('A/177/chi3');g=case['endpoints']['Ca']['composite'][i]
        self.assertGreater(abs(g),.2)
        gate=nl.candidate_gate(True,[0.],[g],0.,{'pass':True})
        self.assertFalse(gate['eligible_for_curvature'])
        # Pure acceptance algebra using two real coarse energy-profile curvatures;
        # this is not a new analytic Hessian or a candidate qualification.
        actual=self.primary['profiles'][0]['finite_difference_diagnostics']['Ca']['central_curvature_kcal_mol_rad2']
        checked=nl.curvature_gate([[actual['0.2']]],[[actual['0.4']]])
        self.assertEqual(checked['status'],'fail')
        self.assertEqual(checked['eigenvalues_unclipped'],[actual['0.2']])
        self.assertIsNone(checked['entropy'])

    def test_native_analytic_recipe_matches_actual_converged_input(self):
        from accommodation_response import recipe
        mp=ROOT/'workspaces/accommodation_response_20260920/prepared_v2/gate/manifest.json'
        m=read_json(mp)
        for t in m['tasks']:
            self.assertEqual(recipe(t['charge'],t['medium'],True),verify(t['input']).read_text())
            self.assertIn('EnGrad',verify(t['input']).read_text())
            self.assertNotIn('NumGrad',verify(t['input']).read_text())
        self.assertEqual(nl.SETTINGS['maxiter'],24);self.assertEqual(nl.SETTINGS['bounds_radian'],[-.8,.8])

    def test_actual_nonlinear_candidates_retain_pass_and_refinement_failure(self):
        # Executed scientific fixtures from the frozen eight-endpoint pilot.
        # Replaying acceptance algebra must preserve a failed check even when
        # the candidate is stationary and its reported curvature is positive.
        directory=ROOT/'workspaces/accommodation_nonlinear_20260920/pilot_v2/endpoints'
        for metal,expected in [('Ca',True),('La',False)]:
            receipt=directory/f'1H4I__{metal}'/'result.json'
            if not receipt.exists():self.skipTest('executed nonlinear fixture unavailable')
            row=read_json(receipt);point=row['candidate']
            work=nl.relative_components(point['components'],row['origin']['components'])
            self.assertEqual(work,row['accommodation_work'])
            gate=nl.candidate_gate(row['optimizer']['success'],point['active_q_radian'],
                point['gradient_kcal_mol_rad'],work['composite_kcal_mol'],point['geometry_checks'])
            self.assertEqual(gate,row['stationarity'])
            self.assertTrue(gate['eligible_for_curvature'])
            h=row['curvature'];replayed=nl.curvature_gate(h['fine_H_kcal_mol_rad2'],h['coarse_H_kcal_mol_rad2'])
            self.assertEqual(replayed['status'],h['status'])
            self.assertEqual(replayed['status']=='pass',expected)
            self.assertEqual(row['stationary_minimum_qualified'],expected)
            components=nl.curvature_components(receipt.parent,row)
            self.assertEqual(components['status'],'complete')
            for name in ('fine_H_kcal_mol_rad2','coarse_H_kcal_mol_rad2'):
                native=np.asarray(components['components']['native_MACE'][name])
                solvent=np.asarray(components['components']['solvent_transfer'][name])
                np.testing.assert_allclose(native+solvent,h[name],atol=1e-10,rtol=0)
            self.assertIsNone(row['numerical_relaxation_free_energy'])
            self.assertIsNone(row['entropy'])

    def test_actual_projected_optimizer_stop_does_not_override_raw_gradient(self):
        receipt=ROOT/'workspaces/accommodation_nonlinear_20260920/pilot_v2/endpoints/PQQSEQ_83440678cbbd658047c9__La/result.json'
        if not receipt.exists():self.skipTest('executed nonlinear fixture unavailable')
        row=read_json(receipt);point=row['candidate']
        self.assertTrue(row['optimizer']['success'])
        self.assertGreater(max(abs(g) for g in point['gradient_kcal_mol_rad']),nl.SETTINGS['gtol'])
        gate=nl.candidate_gate(row['optimizer']['success'],point['active_q_radian'],point['gradient_kcal_mol_rad'],
            row['accommodation_work']['composite_kcal_mol'],point['geometry_checks'])
        self.assertEqual(gate,row['stationarity'])
        self.assertFalse(gate['eligible_for_curvature'])
        self.assertFalse(row['stationary_minimum_qualified'])

if __name__=='__main__':unittest.main()
