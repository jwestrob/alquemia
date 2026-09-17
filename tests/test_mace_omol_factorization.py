"""Exact model factorization checked against actual masked molecular outputs."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,write_new
from mace_file_checks import cached_file_checks
from mace_omol_factorization import verified,TOLERANCE,FACTORIZATION

WORK=ROOT/'workspaces/mace_omol_20260917'
REFERENCE=WORK/'factorization_report_v2/result.json'


@unittest.skipUnless(REFERENCE.exists(),'requires the actual saved-readout factorization audit')
class FactorizationTests(unittest.TestCase):
    @cached_file_checks
    def test_actual_five_case_scores_and_all_detached_variants_agree(self):
        r=verified(REFERENCE)
        self.assertEqual(r['status'],'pass');self.assertEqual(r['new_model_forwards'],0)
        self.assertEqual(len(r['scores']),5);self.assertEqual(len(r['checks']),60)
        self.assertTrue(all(abs(c['error_model_kcal'])<=TOLERANCE for c in r['checks']))
        self.assertTrue(all(abs(s['difference_model_kcal'])<=TOLERANCE for s in r['scores'].values()))

    @cached_file_checks
    def test_corrupted_model_reference_rejected(self):
        r=copy.deepcopy(read_json(REFERENCE));r['Ca_minus_La_disconnected_atom_model_eV']+=1.
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_actual_reference.json';write_new(p,r)
            with self.assertRaisesRegex(InvalidArtifact,'actual qualified readouts'):verified(p)

    @cached_file_checks
    def test_actual_two_bound_reuse_retains_score_without_inventing_detached_components(self):
        try:import openmm
        except ImportError:self.skipTest('prepared-input replay needs the recorded OpenMM driver')
        from mace_omol_prepared import validate,collect
        p=WORK/'prepared_interface_pqq_two_call_v2/manifest.json'
        if not p.exists():self.skipTest('requires the real two-call prepared interface')
        checked=validate(p);self.assertEqual(checked['tasks'],0);self.assertEqual(checked['reused_endpoints'],2)
        r=collect(p);old=read_json(REFERENCE)['scores']['PQQ_1H4I']['four_forward_model_kcal']
        self.assertEqual(r['score_evaluation'],FACTORIZATION)
        self.assertLessEqual(abs(r['R_mask_model_kcal']-old),TOLERANCE)
        self.assertIsNone(r['bound_minus_detached_model_eV'])
        self.assertIsNone(r['calibrated_class']);self.assertEqual(len(r['reused_rows']),2)


if __name__=='__main__':unittest.main()
