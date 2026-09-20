"""Real compact-context forces, mappings and frozen finite response manifests."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import accommodation_response as ar
from affordable_common import read_json,verify,xyz
from mace_site_kinematics import Kinematics
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'workspaces/accommodation_response_20260920/prepared_v2'


class AccommodationResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.d=read_json(WORK/'design.json')

    def test_complete_finite_denominator_and_analytic_inputs(self):
        self.assertEqual(len(self.d['all_tasks']),60)
        self.assertEqual(sum(t['kind']=='center' for t in self.d['all_tasks']),28)
        self.assertEqual(ar.validate(WORK/'gate/manifest.json')['tasks'],4)
        self.assertEqual(ar.validate(WORK/'followup/manifest.json')['tasks'],56)
        for t in self.d['all_tasks']:
            body=verify(t['input']).read_text()
            self.assertEqual('EnGrad' in body,t['kind']=='center')
            self.assertNotIn('NumGrad',body)
            self.assertEqual('ALPB(Water)' in body,t['medium']=='alpb')

    def test_exact_center_force_sources_and_projection(self):
        self.assertEqual(len(self.d['sources']),14)
        for key,source in self.d['sources'].items():
            r=read_json(verify(source['native_MACE_receipt']))
            f=np.load(verify(source['forces']),allow_pickle=False)
            kin=Kinematics(read_json(verify(self.d['maps'][key]))['context'])
            _,positions,_,jac=kin.evaluate(np.zeros(len(kin.modes)))
            np.testing.assert_allclose(positions,[x[1:] for x in xyz(verify(source['xyz']))],atol=1e-10,rtol=0)
            g=-ar.EV_TO_KCAL*np.einsum('mij,ij->m',jac,f)
            np.testing.assert_allclose(g,source['MACE_projected_gradient_kcal_per_unit'],atol=1e-12,rtol=0)
            self.assertEqual(r['status'],'computed')

    def test_physical_modes_and_caps_preserve_bonds(self):
        for key,pin in self.d['maps'].items():
            kin=Kinematics(read_json(verify(pin))['context'])
            self.assertTrue(kin.check(np.full(len(kin.modes),.001))['pass'])
        for t in self.d['all_tasks']:
            if t['kind']!='displacement':continue
            kin=Kinematics(read_json(verify(t['physical_mapping']))['context'])
            q=np.zeros(len(kin.modes));i=next(i for i,m in enumerate(kin.modes) if m['id']==t['mode_id']);q[i]=t['step_radian']
            self.assertLessEqual(kin.displacement(q),.05)
            self.assertTrue(kin.check(q)['pass'])
            self.assertEqual(t['mode_id'],ar.MODES[t['case_id']])

    def test_common_donor_measure_endpoint_water_H_kept_distinct(self):
        for case in ar.CASES:
            ca,la=[self.d['sources'][case+'__'+m] for m in ('Ca','La')]
            self.assertEqual(ca['mode_ids'],la['mode_ids']);self.assertEqual(ca['mode_units'],la['mode_units'])
            self.assertEqual(la['charge']-ca['charge'],1)
            a,b=[xyz(verify(x['xyz'])) for x in (ca,la)]
            self.assertEqual([x[0] for x in a[1:]],[x[0] for x in b[1:]])
            for xa,xb in zip(a,b):
                if xa[0]!='H':self.assertEqual(xa[1:],xb[1:])
            jac=[]
            for m in ('Ca','La'):
                kin=Kinematics(read_json(verify(self.d['maps'][case+'__'+m]))['context'])
                jac.append(kin.evaluate(np.zeros(len(kin.modes)))[3])
            np.testing.assert_allclose(jac[0],jac[1],atol=1e-10,rtol=0)

    def test_actual_native_gradients_include_matching_solvent_driver(self):
        for stage,job in [('gate',1203465),('followup',1203499)]:
            m=read_json(WORK/stage/'manifest.json');c=read_json(WORK/stage/f'collection_{job}.json')
            self.assertEqual(c['status'],'complete')
            for row in c['rows']:
                if row['kind']!='center':continue
                task=next(t for t in m['tasks'] if t['task_id']==row['task_id'])
                parsed=ar.parse_gradient(task,row)
                self.assertFalse(parsed['numerical_gradient_used'])
                np.testing.assert_allclose(parsed['projected_gradient_kcal_per_unit'],row['projected_gradient_kcal_per_unit'],atol=1e-10,rtol=0)
                self.assertEqual(parsed['gradient_kcal_mol_A'],row['gradient_kcal_mol_A'])

    def test_actual_finite_energy_validation_and_composite_sign(self):
        result=read_json(WORK.parent/'result_v2.json')
        self.assertEqual((result['complete_tasks'],result['task_denominator']),(60,60))
        self.assertEqual(len(result['derivative_checks']),12)
        self.assertEqual(result['derivative_check_status'],'pass')
        for r in result['derivative_checks']:
            error=abs(r['central_energy_derivative_kcal_mol_rad']['0.005']-r['analytic_kcal_mol_rad'])
            self.assertAlmostEqual(error,r['analytic_error_kcal_mol_rad'],places=12)
            self.assertLessEqual(error,r['tolerance_kcal_mol_rad'])
        for row in result['cases']:
            delta=row['Ca_minus_La_projected_derivatives']
            np.testing.assert_allclose(np.array(delta['vacuum_MACE'])+delta['solvent_transfer'],delta['composite'],atol=1e-10,rtol=0)
            self.assertIsNone(row['affinity_correction'])
        self.assertIsNone(result['entropy'])


if __name__=='__main__':unittest.main()
