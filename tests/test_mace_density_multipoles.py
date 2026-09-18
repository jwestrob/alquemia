"""Actual framework/density fixtures; algebra and executed integration separated."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,read_json,verify
from mace_density_multipoles import STEPS,points,hessian,monopole_hessian,validate,parse_moments,coupling

W=ROOT/'workspaces/mace_omol_20260917/density_multipole_coupling_v1'
F=ROOT/'workspaces/mace_omol_20260917/qm_electric_field_v1'


class ActualMonopoleAlgebra(unittest.TestCase):
    def test_second_derivative_sign_and_off_diagonal_order(self):
        if not (F/'manifest.json').exists():self.skipTest('actual quantum probe/charge fixture unavailable')
        t=read_json(F/'manifest.json')['tasks'][0];data=read_json(verify(t['probes']))
        # These are analytic potentials of real saved charges, not quantum output.
        centers=np.array(data['positions_bohr'])[:24];source=np.array(data['projected_positions_bohr']);q=np.array(data['projected_charge_e'])
        probes=points(centers);r=np.linalg.norm(probes[:,None,:]-source[None,:,:],axis=2)
        phi=np.sum(q[None,:]/r,axis=1);_,H=hessian(phi,len(centers));exact=monopole_hessian(q,source,centers)
        coarse=np.max(abs(H[0]-exact));fine=np.max(abs(H[1]-exact))
        # Far-field second differences can already be roundoff dominated.
        # This bound uses the actual sum of absolute Coulomb contributions.
        roundoff=64*np.finfo(float).eps*np.max(np.sum(abs(q)[None,:]/r,axis=1))/STEPS[1]**2
        self.assertLess(fine,max(coarse*.4,roundoff))
        self.assertLess(fine,1e-5)
        self.assertTrue(np.array_equal(H,H.transpose(0,1,3,2)))


class NativeAndDensity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('native moment/quantum preparation unavailable')
        cls.m=validate(W/'manifest.json')

    def test_three_actual_native_moment_exports(self):
        for pin in self.m['moment_exports'].values():
            r=read_json(verify(pin));text=verify(r['receipt']['log']).read_text();p=parse_moments(text,r['task'])
            self.assertTrue(p['pass_']);self.assertEqual(p['new_energy_calls'],0)
            self.assertGreater(p['max_absolute_native_quadrupole_trace_eA2'],0)
            damaged='\n'.join(s for s in text.splitlines() if s.split()[:2]!=['GLOBAL_POLE','1'])
            with self.assertRaisesRegex(InvalidArtifact,'missing native moment'):
                parse_moments(damaged,r['task'])

    def test_real_paired_moments_and_units(self):
        for case in {t['case_id'] for t in self.m['tasks']}:
            pair=[t for t in self.m['tasks'] if t['case_id']==case]
            with np.load(verify(pair[0]['multipole_inputs'])) as a,np.load(verify(pair[1]['multipole_inputs'])) as b:
                for key in ('charge','dipole','quadrupole'):np.testing.assert_array_equal(a[key],b[key])

    def test_native_moment_angstrom_formula_matches_atomic_unit_contraction(self):
        t=self.m['tasks'][0];d=read_json(verify(t['probes']))
        with np.load(verify(t['multipole_inputs'])) as a:
            q=a['charge'][:24];mu=a['dipole'][:24];Q=a['quadrupole'][:24]
        src=np.array(d['projected_positions_bohr']);obs=np.array(d['positions_bohr'])[:24];charge=np.array(d['projected_charge_e'])
        delta=obs[:,None,:]-src[None,:,:];r=np.linalg.norm(delta,axis=2)
        phi=np.sum(charge[None,:]/r,axis=1);E=np.sum(charge[None,:,None]*delta/r[:,:,None]**3,axis=1)
        H=monopole_hessian(charge,src,obs);parts=coupling(q,mu,Q,phi,E,H)
        # Independent Cartesian expansion in native eÅ/eÅ² units, including
        # the small actual trace of rounded parameters. No quantum data claimed.
        ra=delta*BOHR_TO_A;rr=r*BOHR_TO_A;ma=mu*BOHR_TO_A;qa=Q*BOHR_TO_A**2
        value=charge[None,:]*(q[:,None]/rr-np.einsum('ni,nki->nk',ma,ra)/rr**3+
            3*np.einsum('nki,nij,nkj->nk',ra,qa,ra)/rr**5-np.trace(qa,axis1=1,axis2=2)[:,None]/rr**3)
        native_units=float(np.sum(value))*HA_TO_KCAL*BOHR_TO_A
        self.assertAlmostEqual(sum(parts[k] for k in ('charge','dipole','quadrupole')),native_units,places=10)

    def test_actual_quantum_utility_and_component_accounting(self):
        reports=sorted(W.glob('report_job_*/result.json'))
        if not reports:self.skipTest('native saved-density integrations unrun')
        r=read_json(reports[-1]);self.assertTrue(r['complete']);self.assertIsNone(r['numerical_score'])
        for row in r['rows'].values():
            receipt=read_json(verify(row['execution_receipt']));self.assertEqual(receipt['status'],'complete');verify(receipt['potential'])
            c=row['components_fine_kcal_mol']
            self.assertEqual(row['total_component_diagnostic_kcal_mol'],sum(c[k] for k in ('charge','dipole','quadrupole')))
            self.assertEqual(row['quadrupole_refinement_pass'],abs(row['quadrupole_refinement_kcal_mol'])<=.02)
        for case,p in r['pairs'].items():
            c=r['rows'][case+'_Ca']['components_fine_kcal_mol'];l=r['rows'][case+'_La']['components_fine_kcal_mol']
            self.assertEqual(p['Ca_minus_La_component_diagnostic_kcal_mol'],{k:c[k]-l[k] for k in c})


if __name__=='__main__':unittest.main()
