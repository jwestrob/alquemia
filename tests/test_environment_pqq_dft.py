"""Finite native check using pinned actual PQQ preparations and CPCM outputs."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL, read_json, verify
from environment_pqq_dft import CASES, compare_state, validate


class PQQNativeContext(unittest.TestCase):
    def test_frozen_eight_actual_context_inputs(self):
        path=ROOT/'workspaces/environment_pqq_dft_20260919/prepared_v2/manifest.json'
        self.assertEqual(validate(path)['tasks'],8)
        m=read_json(path)
        self.assertEqual({t['case'] for t in m['tasks']},set(CASES))
        for t in m['tasks']:
            self.assertEqual(t['xyz']['sha256'],t['source_OMOL_task']['xyz']['sha256'])
            self.assertEqual(t['charge'],t['source_OMOL_task']['charge'])
        self.assertIsNone(m['reference']); self.assertIsNone(m['calibrated_decision'])

    def test_actual_neutral_cpcm_contrast_and_diagnostic_algebra(self):
        old=read_json(ROOT/'workspaces/second_shell_20260919/prepared_v2/collection_1202082.json')
        cfg=read_json(ROOT/'diagnostics/second_shell_20260919/CONFIG.json')
        mace=read_json(ROOT/'workspaces/environment_pqq_20260919/prepared_v1/result_1202474.json')
        for case in ('1H4I','4MAE'):
            c=next(s['endpoints'] for s in cfg['cases'] if s['case']==case)
            prior=next(r for r in old['rows'] if r['case']==case)
            row=compare_state(case,c,prior['endpoints'],next(r for r in mace['rows'] if r['case']==case))
            self.assertAlmostEqual(row['delta_R_kcal_mol'],prior['delta_R_kcal_mol'],places=8)
            terms=row['CPCM_dielectric_diagnostic']
            self.assertAlmostEqual(row['CPCM_dielectric_delta_contrast_kcal_mol'],(terms['Ca']['delta_hartree']-terms['La']['delta_hartree'])*HA_TO_KCAL,places=10)
            self.assertTrue(row['CPCM_term_already_in_total_energy'])
            self.assertIsNone(row['unique_physical_component_attribution'])


if __name__=='__main__':unittest.main()
