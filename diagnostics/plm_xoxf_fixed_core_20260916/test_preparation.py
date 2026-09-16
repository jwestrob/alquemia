"""Read-only regression checks for selected models and actual prepared cores."""
import copy
import unittest
import select_inputs as inputs
from geometry_checks import check,unchanged_non_H

class PreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pp=inputs.read(inputs.OUT/'prepared/prepared_pairs.json')
        cls.carves=[inputs.read(c['carve_manifest']['path']) for c in cls.pp['cases']]
    def test_two_exact_targets(self):
        self.assertEqual({c['target_id'] for c in self.pp['cases']},set(inputs.EXPECTED))
        self.assertTrue(all(c['rank']==0 for c in self.pp['cases']))
    def test_actual_cores_pass(self):
        for m in self.carves:
            d=check(m['qm_fragments']);self.assertEqual(d['status'],'PASS');self.assertEqual(d['total_core_atom_count'],80)
            self.assertEqual(d['formula_excluding_metal'],{'C':27,'H':31,'N':6,'O':15})
    def test_old_failed_hydrogens_rejected(self):
        p=inputs.A/'workspaces/plm_adh9_af3_20260916/preparation/pairs/PQQSEQ_13d74836d4b7a3e02140_AF3_sample1/PQQSEQ_13d74836d4b7a3e02140_AF3_sample1_carve_manifest.json'
        self.assertEqual(check(inputs.read(p)['qm_fragments'])['status'],'FAIL')
    def test_core_overlap_rejected(self):
        f=copy.deepcopy(self.carves[0]['qm_fragments']);a=[a for g in f for a in g['atom_records'] if a['element']=='H']
        a[1]['xyz_A']=list(a[0]['xyz_A']);self.assertEqual(check(f)['status'],'FAIL')
    def test_nonfinite_rejected(self):
        f=copy.deepcopy(self.carves[0]['qm_fragments']);f[0]['atom_records'][0]['xyz_A'][0]=float('nan')
        self.assertEqual(check(f)['status'],'FAIL')
    def test_heavy_movement_rejected(self):
        f=self.carves[0]['qm_fragments'];new=copy.deepcopy(f);new[0]['atom_records'][0]['xyz_A'][0]+=.001
        with self.assertRaisesRegex(ValueError,'Frozen'):unchanged_non_H(f,new)
    def test_selection_tie_break(self):
        rows=[{'CN':7,'protein_La_iptm':.98,'protein_CN':4,'sample':i} for i in [2,1,0]]
        self.assertEqual(inputs.select(rows)['sample'],0)
    def test_CN_has_priority_without_admission_rescue(self):
        rows=[{'CN':7,'protein_La_iptm':.99,'protein_CN':4,'sample':0},{'CN':8,'protein_La_iptm':.8,'protein_CN':5,'sample':1}]
        self.assertEqual(inputs.select(rows)['sample'],1)
    def test_pair_coordinates_and_charges(self):
        for m in self.carves:
            la=inputs.Path(m['outputs']['La_xyz']['path']).read_text().splitlines();ca=inputs.Path(m['outputs']['Ca_xyz']['path']).read_text().splitlines()
            self.assertEqual(la[3:],ca[3:]);self.assertEqual(la[2].split()[1:],ca[2].split()[1:])
            self.assertEqual((m['charge_ledger']['La_total'],m['charge_ledger']['Ca_total']),(-2,-3))

if __name__=='__main__':unittest.main(verbosity=2)
