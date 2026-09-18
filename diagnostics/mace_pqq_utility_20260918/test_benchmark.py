"""Parser/algebra regression using the real frozen PQQ results; no inference."""
import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('benchmark',Path(__file__).with_name('benchmark.py'))
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)


class RealReferenceTests(unittest.TestCase):
    def test_dft_algebra_and_released_primary_bands(self):
        ref=b.read(ROOT/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json')
        self.assertEqual(len(ref['scores']),25)
        for row in ref['scores']:
            score=b.dft_score(row['energies_hartree'],ref)
            self.assertEqual(score['score'],row['S_aquo_gauge_kcal_mol'])
            self.assertEqual(score['raw_R_kcal_mol'],row['R_kcal_mol'])
            self.assertEqual(score['class'],row['class'])

    def test_masked_two_call_algebra_bands_and_missing_case(self):
        ref=b.read(ROOT/'workspaces/mace_omol_20260917/masked_calibration_v1/reference.json')
        bands=ref['bands'];missing=[]
        for r in ref['scores']:
            if r['status']!='computed':
                missing.append(r['case_id']);self.assertIsNone(r['two_call_model_kcal']);continue
            bound=r['bound_model_eV'];free=r['disconnected_atom_model_eV']
            score=(bound['Ca']-bound['La']-(free['Ca']-free['La']))*23.06054783061903
            self.assertEqual(score,r['two_call_model_kcal'])
            self.assertEqual(b.classify(score,bands['Ca_max_inclusive_model_kcal'],
                                       bands['La_min_inclusive_model_kcal']),r['expected_class'])
        self.assertEqual(missing,['1KB0'])


if __name__=='__main__':unittest.main()
