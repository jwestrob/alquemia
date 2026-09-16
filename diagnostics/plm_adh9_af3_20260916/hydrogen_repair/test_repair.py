"""Actual-artifact regressions: failed H preparation must never be admitted."""
import copy
import json
import unittest
import gemmi
import prepare_repaired as prep
from validate_geometry import validate

class RepairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        wrapper=prep.load_wrapper();fixed=wrapper.fixed
        cls.old=prep.read(prep.helper.OLD/(prep.helper.STEM+'_carve_manifest.json'))['qm_fragments']
        cls.min=prep.read(prep.BASE/'recovered_exact_objective/minimization.json')
        structure=gemmi.read_structure(str(prep.BASE/'recovered_exact_objective/protonated_recovered.pdb'))
        target=prep.read(prep.ROOT/'workspaces/plm_adh9_af3_20260916/hydrogen_retry_01/candidate_manifest.json')['targets'][0]
        keys,residues,partner=wrapper.role_state(structure[0],target)
        site,*_=wrapper.site_state(structure[0],target);cls.new=[copy.deepcopy(cls.old[0])]
        for role in fixed.ROLE_ORDER:
            if role=='catalytic_asp_cationic_partner':_,f=fixed.cationic_sidechain_fragment(residues[role],keys[role],site.atom.pos)
            else:_,f=fixed.canonical_sidechain_fragment(residues[role],keys[role],role,site.atom.pos)
            cls.new.append(f)
    def test_original_overlaps_rejected(self):
        with self.assertRaisesRegex(ValueError,'Overlapping'):validate(self.old,self.old,self.min)
    def test_repaired_actual_core_passes(self):
        d=validate(self.old,self.new,self.min);self.assertEqual(d['changed_source_hydrogens'],23);self.assertGreater(d['minimum_core_HH_distance_A'],.5)
    def mutate_coordinate(self,origin):
        new=copy.deepcopy(self.new)
        for f in new:
            for a in f['atom_records']:
                if a['origin']==origin:a['xyz_A'][0]+=.001;return new
        self.fail('Required atom missing')
    def test_heavy_move_rejected(self):
        with self.assertRaisesRegex(ValueError,'Frozen'):validate(self.old,self.mutate_coordinate('source_heavy_atom'),self.min)
    def test_PQQ_hydrogen_move_rejected(self):
        with self.assertRaisesRegex(ValueError,'Frozen'):validate(self.old,self.mutate_coordinate('generated_opposite_bisector'),self.min)
    def test_cap_move_rejected(self):
        with self.assertRaisesRegex(ValueError,'Frozen'):validate(self.old,self.mutate_coordinate('generated_Cbeta_link_cap'),self.min)
    def test_identity_move_rejected(self):
        new=copy.deepcopy(self.new);new[1]['atom_records'][0]['name']='WRONG'
        with self.assertRaisesRegex(ValueError,'identity/order'):validate(self.old,new,self.min)
    def test_charge_change_rejected(self):
        new=copy.deepcopy(self.new);new[1]['formal_charge']=0
        with self.assertRaisesRegex(ValueError,'charge'):validate(self.old,new,self.min)
    def test_unconverged_rejected(self):
        m=copy.deepcopy(self.min);m['final']['hydrogen_force_rms_kJ_mol_nm']=2
        with self.assertRaisesRegex(ValueError,'converged'):validate(self.old,self.new,m)
    def test_nonfinite_rejected(self):
        new=copy.deepcopy(self.new);new[1]['atom_records'][5]['xyz_A'][0]=float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):validate(self.old,new,self.min)

if __name__=='__main__':unittest.main(verbosity=2)
