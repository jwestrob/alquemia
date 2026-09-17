"""Actual reference derivation and compatibility checks, never synthetic energies."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,write_new
from mace_file_checks import cached_file_checks
from mace_omol_mask_calibration import calculate,verified,decision

WORK=ROOT/'workspaces/mace_omol_20260917'
REFERENCE=WORK/'masked_calibration_v1/reference.json'
FACTOR=WORK/'factorization_report_v2/result.json'


class CalibrationGuardTests(unittest.TestCase):
    def test_native_canonical_report_cannot_supply_masked_reference(self):
        p=WORK/'intact_panel_report_v1/result.json'
        if not p.exists():self.skipTest('requires actual failed native canonical report')
        with self.assertRaisesRegex(InvalidArtifact,'masked-descriptor report'):calculate(p,FACTOR)


@unittest.skipUnless(REFERENCE.exists(),'requires actual completed canonical descriptor and factorization reference')
class ActualCalibrationTests(unittest.TestCase):
    @cached_file_checks
    def test_all_cases_and_unavailable_transfer_retained(self):
        r=verified(REFERENCE)
        self.assertEqual(len(r['scores']),28);self.assertEqual(r['calibration_total_count'],25)
        self.assertEqual(r['transfer_total_count'],3)
        broken=next(s for s in r['scores'] if s['case_id']=='1KB0')
        self.assertEqual(broken['status'],'unavailable');self.assertIsNone(broken['two_call_model_kcal'])
        self.assertIsNone(broken['calibrated_class']);self.assertEqual(r['new_model_forwards'],0)

    @cached_file_checks
    def test_corrupted_actual_reference_policy_rejected(self):
        r=copy.deepcopy(read_json(REFERENCE));r['decision_policy']['minimum_gap_kcal_mol']+=1.
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_real_reference.json';write_new(p,r)
            with self.assertRaisesRegex(InvalidArtifact,'actual source calculations'):verified(p)

    @cached_file_checks
    def test_pqq_replay_uses_its_recorded_transfer_decision(self):
        try:import openmm
        except ImportError:self.skipTest('prepared scoring needs recorded OpenMM driver')
        from mace_omol_prepared import collect
        ref=read_json(REFERENCE)
        if ref['status']!='pass':self.skipTest('actual calibration unavailable; no successful decision is fabricated')
        p=WORK/'prepared_interface_pqq_two_call_v2/manifest.json';m=read_json(p);score=collect(p)
        result=decision(score,m,REFERENCE)
        source=next(s for s in ref['scores'] if s['case_id']=='1H4I')
        self.assertEqual(result['calibrated_class'],source['calibrated_class'])
        self.assertEqual(result['decision_status'],'research_calibrated_PQQ_functional_association')
        wrong=copy.deepcopy(score);wrong['score_evaluation']='four_full_bound_detached_states_v1'
        self.assertEqual(decision(wrong,m,REFERENCE)['decision_status'],'incompatible_model_or_numeric_evaluation')

    @cached_file_checks
    def test_real_ggr_never_inherits_pqq_band(self):
        try:import openmm
        except ImportError:self.skipTest('prepared scoring needs recorded OpenMM driver')
        from mace_omol_prepared import collect
        ref=read_json(REFERENCE)
        if ref['status']!='pass':self.skipTest('actual calibration unavailable')
        p=WORK/'prepared_interface_ggr_two_call_v1/manifest.json'
        if not p.exists():self.skipTest('requires actual GGR reuse manifest')
        result=decision(collect(p),read_json(p),REFERENCE)
        self.assertIsNone(result['calibrated_class'])
        self.assertEqual(result['decision_status'],'outside_canonical_PQQ_calibration_scope')


if __name__=='__main__':unittest.main()
