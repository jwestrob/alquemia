"""Real source states and actual native boundary initializations."""
from pathlib import Path
import copy
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_density_gk_boundary import boundary,validate,verify_native,comparable

W=ROOT/'workspaces/mace_omol_20260917/density_gk_boundary_v1'


class Boundary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('real native boundary preparation absent')
        cls.m=validate(W/'manifest.json')

    def test_real_graph_charge_ledger_and_corrupted_formal_state(self):
        for pin in self.m['cases'].values():
            meta=read_json(verify(pin));s=read_json(verify(meta['state']));mapping=read_json(verify(meta['base_mapping']));p=read_json(verify(meta['base_parameters']))
            b=boundary(s,mapping,p);self.assertEqual(b['ledger'],meta['ledger'])
            self.assertAlmostEqual(b['actual_environment_charge_e'],b['environment_formal_charge_e'],places=9)
            self.assertTrue(all(b['environment_charge_by_id'][i]==0 for i in b['source_support_ids']))
            bad=copy.deepcopy(s);rid=next(r['residue'] for r in s['boundary_ledger'] if r['residue'] in s['ligand_ledger'])
            bad['ligand_ledger'][rid]+=1
            with self.assertRaisesRegex(InvalidArtifact,'chemical formal-charge ledger'):
                boundary(bad,mapping,p)

    def test_sixteen_real_native_preparations_and_source_atom_retention(self):
        for t in self.m['tasks']:
            r=read_json(W/'tasks'/t['task_id']/'result.json');meta=read_json(verify(t['boundary']))
            p=verify_native(verify(r['receipt']['log']).read_text(),t,meta,read_json(verify(meta['base_parameters'])))
            self.assertEqual(p['status'],'native_density_boundary_prepared');self.assertTrue(p['pass_'])
            n=len(meta['physical_ids']);self.assertEqual(p['native_multipole_map'][n],[n,n,13])
            self.assertEqual(p['parameters']['atoms'][-1]['polarizability_A3'],0)
            self.assertIsNone(p['numerical_score'])

    def test_paired_actual_environment_and_cavity_identity(self):
        for case in self.m['cases']:
            ca=read_json(W/'tasks'/(case+'_Ca_zero_source')/'result.json')
            la=read_json(W/'tasks'/(case+'_La_zero_source')/'result.json')
            self.assertEqual(comparable(ca),comparable(la))
            damaged=copy.deepcopy(la);damaged['parameters']['atoms'][-1]['radius_A']+=.01
            self.assertNotEqual(comparable(ca),comparable(damaged))


if __name__=='__main__':unittest.main()
