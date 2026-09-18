"""Real-fixture coupled kinematics and finite-manifest regressions."""
from pathlib import Path
import copy,json,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,InvalidArtifact
from mace_site_kinematics import Kinematics
from mace_site_path import validate
BASE=ROOT/'workspaces/mace_site_response_20260918'
P=BASE/'coupled_prepared_v1/preparation.json'

@unittest.skipUnless(P.exists(),'requires actual coupled paths prepared from archived gradients')
class CoupledPathTests(unittest.TestCase):
    def test_analytic_composed_jacobians_and_actual_center_gradient(self):
        p=read_json(P)
        for key in ('GGR_2FW0_connected_Ca','GGR_2FVY_extended_La','ALPHA_6IP9_La'):
            s=p['states'][key];kin=Kinematics(read_json(verify(s['kinematics'])));row=read_json(verify(s['source_projection']));zero=np.zeros(len(kin.modes))
            full,core,j,cj=kin.evaluate(zero)
            np.testing.assert_allclose(j,np.load(verify(row['full_jacobian'])),atol=1e-12,rtol=0)
            np.testing.assert_allclose(cj,np.load(verify(row['core_jacobian'])),atol=1e-12,rtol=0)
            q=np.array(s['direction']['q']);self.assertLess(np.array(row['hybrid_projected_gradient'])@q,0)
            for fraction in (.5,1.):
                check=kin.check(q*fraction);self.assertTrue(check['pass'],check);self.assertLessEqual(kin.displacement(q*fraction),.20+1e-12)
                np.testing.assert_array_equal(kin.positions_only(q*fraction),kin.evaluate(q*fraction)[0])
            # Archived full gradient independently projects through the map at 0.
            g=np.load(verify(row['hybrid_physical_gradient']))
            np.testing.assert_allclose(np.einsum('mij,ij->m',j,g),row['hybrid_projected_gradient'],atol=1e-9,rtol=0)

    @unittest.skipUnless((BASE/'coupled_initial_v1/initial.json').exists(),'finite initial manifest not yet prepared')
    def test_real_manifest_and_explicitly_corrupted_charge(self):
        initial=read_json(BASE/'coupled_initial_v1/initial.json');path=verify(initial['jobs'][0]);m=read_json(path)
        self.assertEqual(validate(path)['status'],'pass')
        broken=copy.deepcopy(m);broken['tasks'][0]['charge']+=1
        with tempfile.TemporaryDirectory() as d:
            q=Path(d)/'corrupted_real_manifest.json';q.write_text(json.dumps(broken))
            with self.assertRaises(InvalidArtifact):validate(q)

    @unittest.skipUnless((BASE/'coupled_selection_v1/result.json').exists(),'actual finite energy grid not yet collected')
    def test_actual_selection_and_corrupted_selected_point(self):
        from mace_site_path_native import check_selection
        a=read_json(BASE/'coupled_selection_v1/result.json');p=read_json(P)
        check_selection(a,p)
        broken=copy.deepcopy(a);key=next(iter(broken['selected']))
        chosen=broken['selected'][key]['point']
        broken['selected'][key]['point']='p025' if chosen!='p025' else 'p100'
        with self.assertRaises(InvalidArtifact):check_selection(broken,p)

    @unittest.skipUnless((BASE/'coupled_report_v1/result.json').exists(),'actual native calculations not yet reported')
    def test_actual_native_energy_algebra_and_report_replay(self):
        from affordable_common import HA_TO_KCAL
        from mace_hybrid import EV_TO_KCAL
        from mace_site_path_native import report
        old=read_json(BASE/'coupled_report_v1/result.json');selection=read_json(verify(old['selection']))
        for name,s in old['scores'].items():
            energies={}
            for metal in ('Ca','La'):
                key=name+'_'+metal;point=selection['selected'][key]['point'];q=old['quantum']['rows'][key]
                f,c=[old['MACE']['rows'][kind+'__'+key+'__'+point]['energy_eV'] for kind in ('full','core')]
                energies[metal]=(q['energy_hartree'],f,c)
            ca,la=energies['Ca'],energies['La']
            direct=(ca[0]-la[0])*HA_TO_KCAL+((ca[1]-la[1])-(ca[2]-la[2]))*EV_TO_KCAL
            self.assertAlmostEqual(direct,s['actual_R_kcal'],places=6)
        self.assertEqual(len(old['contrasts']),12)
        with tempfile.TemporaryDirectory() as d:
            report(BASE/'coupled_native_v1/preparation.json',Path(d)/'replay');new=read_json(Path(d)/'replay/result.json')
            for key in ('rows','scores','partition','contrasts','raw_pass_count','qualified_pass_count'):
                self.assertEqual(old[key],new[key])

if __name__=='__main__':unittest.main()
