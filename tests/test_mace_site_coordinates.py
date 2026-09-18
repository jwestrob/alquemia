"""Real matched-hybrid donor kinematics, caps and gradient projection."""
from pathlib import Path
import copy,sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,InvalidArtifact
from mace_site_coordinates import SiteCoordinates
from mace_hybrid import EV_TO_KCAL
P=ROOT/'workspaces/mace_site_response_20260918/prepared_v1/preparation.json'

@unittest.skipUnless(P.exists(),'requires real archived matched-hybrid gradients and source preparations')
class SiteCoordinateTests(unittest.TestCase):
    def test_real_bonds_and_caps_under_all_physical_coordinates(self):
        p=read_json(P)
        for key,ref in p['rows'].items():
            if not key.endswith('_Ca'):continue
            row=read_json(verify(ref));s=row['source'];g=SiteCoordinates(read_json(verify(s['repair'])),read_json(verify(s['whole_preparation'])),verify(s['core_state']['xyz']))
            self.assertEqual(g.modes,row['modes']);self.assertTrue(g.checks()['pass'])
            np.testing.assert_array_equal(g.full_jacobian,np.load(verify(row['full_jacobian'])))
            np.testing.assert_array_equal(g.core_jacobian,np.load(verify(row['core_jacobian'])))
            # A simultaneous small rotation tests coupled bond preservation,
            # without inventing a protein, force, energy or curvature model.
            q=np.array([m['linear_probe_amplitude']/4 for m in g.modes]);moved=g.full_positions(q)
            for i,j in g.bonds:self.assertAlmostEqual(np.linalg.norm(moved[i]-moved[j]),np.linalg.norm(g.positions[i]-g.positions[j]),places=9)
            for i,a in enumerate(g.atoms):
                if a['kind']=='retained_site_water':np.testing.assert_array_equal(moved[i],g.positions[i])

    def test_actual_gradient_chain_rule_and_sign(self):
        p=read_json(P)
        for key,ref in p['rows'].items():
            row=read_json(verify(ref));s=row['source'];g=SiteCoordinates(read_json(verify(s['repair'])),read_json(verify(s['whole_preparation'])),verify(s['core_state']['xyz']))
            gd=np.array(s['DFT']['gradient_kcal_mol_per_A']);gc=np.load(verify(s['MACE_core']['gradient']))*EV_TO_KCAL;gf=np.load(verify(s['MACE_full']['gradient']))*EV_TO_KCAL
            expected=np.einsum('mij,ij->m',g.core_jacobian,gd-gc)+np.einsum('mij,ij->m',g.full_jacobian,gf)
            np.testing.assert_allclose(expected,row['hybrid_projected_gradient'],atol=1e-10,rtol=0)
            np.testing.assert_allclose(expected,np.einsum('mij,ij->m',g.full_jacobian,g.source_gradient(gd-gc)+gf),atol=1e-9,rtol=0)
            self.assertIsNone(row['score_correction'])
        r=read_json(P.parent.parent/'report_v1/result.json')
        for name,pair in r['pairs'].items():
            ca,la=[read_json(verify(p['rows'][name+'_'+m])) for m in ('Ca','La')]
            np.testing.assert_allclose(pair['projected_gradient_R'],np.array(ca['hybrid_projected_gradient'])-la['hybrid_projected_gradient'],atol=1e-10,rtol=0)
            self.assertFalse(pair['is_energy_evaluation'])

    def test_corrupted_real_connectivity_is_rejected(self):
        p=read_json(P);row=read_json(verify(p['rows']['GGR_1GLG_extended_Ca']));s=row['source'];w=read_json(verify(s['whole_preparation']));repair=read_json(verify(s['repair']))
        broken=copy.deepcopy(w);broken['preparation_details']['bonds']=[b for b in broken['preparation_details']['bonds'] if {b['atom_a_id'],b['atom_b_id']}!={'A/134//CA','A/134//CB'}]
        with self.assertRaises(InvalidArtifact):SiteCoordinates(repair,broken,verify(s['core_state']['xyz']))

if __name__=='__main__':unittest.main()
