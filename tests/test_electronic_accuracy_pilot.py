"""Integrity checks on the real archived six-endpoint panel; no computed fixtures."""
import sys
from pathlib import Path
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from electronic_accuracy_pilot import sources,electronic_input
from affordable_common import xyz,verify,energy
from electronic_accuracy_collect import compare,cc_output
from electronic_accuracy_autocollect import report,finalize

class RealElectronicPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rows,cls.parent=sources(ROOT)

    def test_real_source_energies_and_charge_pairs(self):
        self.assertEqual(len(self.rows),6)
        for r in self.rows:
            self.assertAlmostEqual(energy(verify(r['baseline']['output'])),r['baseline']['energy_hartree'],places=9)
        for case in {r['case'] for r in self.rows}:
            rows={r['metal']:r for r in self.rows if r['case']==case}
            self.assertEqual(rows['La']['charge']-rows['Ca']['charge'],1)
            ca,la=[xyz(verify(rows[m]['xyz'])) for m in ('Ca','La')]
            self.assertEqual([a[0] for a in ca[1:]],[a[0] for a in la[1:]])
            self.assertTrue(np.array_equal(np.array([a[1:] for a in ca if a[0]!='H']),np.array([a[1:] for a in la if a[0]!='H'])))

    def test_context_preparation_only_moves_water_h(self):
        allowed={35,36,38,39}
        for metal in ('Ca','La'):
            rows={r['case']:r for r in self.rows if r['metal']==metal}
            a=xyz(verify(rows['ALPHA_1F6S_original']['xyz']))
            b=xyz(verify(rows['ALPHA_1F6S_water_prepared']['xyz']))
            self.assertEqual(len(a),40)
            for i,(x,y) in enumerate(zip(a,b)):
                self.assertEqual(x[0],y[0])
                if i not in allowed:self.assertEqual(x,y)
                else:self.assertEqual(x[0],'H')

    def test_same_explicit_correlated_method_no_composite_corrections(self):
        for r in self.rows:
            s=electronic_input(r['charge'],r['metal'])
            self.assertTrue(s.endswith('\n'))
            self.assertIn('DLPNO-CCSD(T1) TightPNO',s)
            self.assertIn('CPCMccm 2',s)
            self.assertIn('NewNCore Ca 10',s)
            self.assertIn('NewNCore La 46',s)
            self.assertNotIn('D3',s)
            self.assertNotIn('gCP',s)

    def test_archived_preparation_effect_and_sign(self):
        es={}
        for r in self.rows:es.setdefault(r['case'],{})[r['metal']]=energy(verify(r['baseline']['output']))
        result=compare(es)
        self.assertAlmostEqual(result['alpha_minus_GGR_kcal_mol']['ALPHA_1F6S_original'],-14.5592454084387,places=7)
        self.assertAlmostEqual(result['alpha_minus_GGR_kcal_mol']['ALPHA_1F6S_water_prepared'],10.923654263172676,places=7)
        self.assertAlmostEqual(result['water_preparation_delta_R_kcal_mol'],25.482899671611374,places=7)
        self.assertIsNone(result['calibrated_classification'])

    def test_baseline_cannot_replace_cc_output(self):
        with self.assertRaisesRegex(ValueError,'missing correlated'):
            cc_output(verify(self.rows[0]['baseline']['output']))

    def test_live_incomplete_cc_has_no_derived_energy(self):
        from electronic_accuracy_collect import collect
        from electronic_accuracy_autocollect import scheduler
        m=ROOT/'workspaces/electronic_accuracy_20260919/direct_v1/manifest.json'
        result=collect(m);actual=scheduler('1202429')
        body=report(result,actual,m)
        self.assertIn(f"{result['complete_endpoints']}/6",body)
        for r in result['rows']:
            if r['status']!='complete':self.assertIsNone(r['result'])
        if result['complete_endpoints']<6:self.assertIn('full electronic-method comparison is unavailable',body)
        # Actual current running state: finalizer must refuse before it even
        # reads a destination, sends an email or attempts cleanup.
        if not actual['terminal']:
            with self.assertRaisesRegex(RuntimeError,'terminal job'):
                finalize({},actual)

if __name__=='__main__':unittest.main()
