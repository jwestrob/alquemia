"""Actual archived La geometry/gradient algebra, not molecular qualification."""
import sys,unittest,copy
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import read_json,verify,InvalidArtifact
from affordable_response import read_engrad
from metal_environment_reference import read_pcgrad
import metal_environment_force_assembly as f
ROOT=Path(__file__).resolve().parents[1]
INPUTS=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2/INPUTS.json'
COLLECTION=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_hans_scout_v1/FINAL_COLLECTION.json'
class ActualLanMAssembly(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mp=f.prepare_mapping(INPUTS,'A','La');r=read_json(COLLECTION)['rows']['La_A'];cls.q=read_engrad(verify(r['engrad']))['gradient_Ha_per_bohr'];cls.pc=read_pcgrad(verify(r['pointcharge_gradient']),len(cls.mp['field']))
    def test_both_real_archived_configurations_and_no_silent_classical_zero(self):
        for state in ('A','B'):
            r=f.assemble(INPUTS,state,'La',COLLECTION,'La_'+state);self.assertEqual(len(r['physical_atoms']),1891);self.assertIsNone(r['total_gradient_kcal_mol_A']);self.assertEqual(r['classical_status'],'unavailable');self.assertLess(np.max(np.abs(r['electronic']['mapping_translation_conservation_residual_kcal_mol_A'])),1e-10)
    def test_real_caps_both_anchor_chain_rule(self):
        q=np.array([a['xyz_A'] for a in self.mp['physical_atoms']]);core=self.mp['core'];delta=np.array([.173,.117,.231]);delta/=np.linalg.norm(delta);h=1e-5
        for cap in self.mp['caps']:
            g=self.q[cap['qm_index']];length=core[cap['qm_index']]['length_A'];ri=cap['retained_index'];oi=cap['omitted_index']
            def energy(r,o):
                # Linear contraction with an ACTUAL archived cap gradient;
                # algebraic chain-rule functional, not a molecular energy.
                return float(g@(r+length*(o-r)/np.linalg.norm(o-r)))
            for anchor,j in [('retained','jacobian_retained'),('omitted','jacobian_omitted')]:
                r=q[ri];o=q[oi]
                if anchor=='retained':fd=(energy(r+h*delta,o)-energy(r-h*delta,o))/(2*h)
                else:fd=(energy(r,o+h*delta)-energy(r,o-h*delta))/(2*h)
                self.assertAlmostEqual(fd,float(g@np.array(cap[j])@delta),places=8)
    def test_actual_gradient_rigid_covariance_algebra(self):
        theta=.37;rot=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
        mm=copy.deepcopy(self.mp)
        for cap in mm['caps']:
            for k in ('jacobian_retained','jacobian_omitted'):cap[k]=(rot@np.array(cap[k])@rot.T).tolist()
        a=f.scatter(self.mp,self.q,self.pc);b=f.scatter(mm,self.q@rot.T,self.pc@rot.T);np.testing.assert_allclose(b['gradient_kcal_mol_A'],np.array(a['gradient_kcal_mol_A'])@rot.T,atol=1e-11,rtol=0)
    def test_wrong_source_state_and_missing_electronic_rejected(self):
        with self.assertRaisesRegex(InvalidArtifact,'geometry/field'):f.assemble(INPUTS,'B','La',COLLECTION,'La_A')
        with self.assertRaisesRegex(InvalidArtifact,'unavailable'):f.assemble(INPUTS,'A','Dy',COLLECTION,'Dy_A')
    def test_classical_undeclared_coulomb_rejected(self):
        r=f.assemble(INPUTS,'A','La',COLLECTION,'La_A')
        with self.assertRaisesRegex(InvalidArtifact,'exclusions'):f.add_classical(r,dict(quantity='gradient',units='kcal/mol/angstrom',QM_MM_Coulomb_included=True,C4_induction_included=False))
if __name__=='__main__':unittest.main()
