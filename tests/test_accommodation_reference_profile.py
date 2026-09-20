"""One actual prepared reference and its finite computed grid; no invented energy."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import read_json,verify,xyz
from mace_hybrid import EV_TO_KCAL
from compact_solvation import completed
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'workspaces/accommodation_reference_profile_20260920'

class ReferenceProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.d=read_json(WORK/'prepared_v1/design.json')

    def test_fixed_source_pair_and_complete_grid(self):
        self.assertEqual(self.d['source_conditioning_metal'],'Ca')
        self.assertEqual(self.d['new_calls'],{'MACE':4,'GFN2':8,'DFT':0})
        self.assertEqual({(p['metal'],p['angle_radian']) for p in self.d['points']},{(m,q) for m in ('Ca','La') for q in (-.2,.2)})
        for q in (-.2,.2):
            pair={p['metal']:p for p in self.d['points'] if p['angle_radian']==q}
            ca,la=[xyz(verify(pair[z]['xyz'])) for z in ('Ca','La')]
            self.assertEqual(ca[1:],la[1:]);self.assertEqual(ca[0][1:],la[0][1:])
            self.assertEqual(pair['La']['charge']-pair['Ca']['charge'],1)
            self.assertTrue(all(p['geometry_checks']['pass'] for p in pair.values()))

    def test_exact_archived_zero_receipts(self):
        for z,c in self.d['centers'].items():
            native=read_json(verify(c['MACE_receipt']))
            self.assertEqual(native['energy_eV'],c['MACE_energy_eV'])
            for medium,pin in c['GFN2'].items():
                current=completed(verify(pin['manifest']),pin['task_id'])
                self.assertIsNotNone(current);self.assertEqual(current['energy_hartree'],pin['energy_hartree'])

    def test_actual_grid_lowering_and_unit_sign(self):
        r=read_json(WORK/'result_v1.json');self.assertEqual(len(r['points']),6)
        for model,result in r['models'].items():
            self.assertEqual(result['status'],'available')
            field=model+'_kcal_mol';group={z:[p for p in r['points'] if p['metal']==z] for z in ('Ca','La')}
            actual=min(p[field] for p in group['Ca'])-min(p[field] for p in group['La'])
            self.assertAlmostEqual(actual,result['R_grid_kcal_mol'],places=9)
            self.assertAlmostEqual(result['endpoint_lowering_kcal_mol']['Ca']-result['endpoint_lowering_kcal_mol']['La'],result['delta_R_grid_kcal_mol'],places=8)
            self.assertTrue(result['old_band_transfer_only']);self.assertIsNone(result['free_energy'])
        self.assertFalse(r['threshold_refitted']);self.assertFalse(r['baseline_changed'])

if __name__=='__main__':unittest.main()
