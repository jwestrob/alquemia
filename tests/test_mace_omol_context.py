"""Actual partition cancellation and rejection of altered archived results."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,write_new
from mace_hybrid import EV_TO_KCAL
from mace_omol_context import partition

class ContextPartition(unittest.TestCase):
    def setUp(self):
        self.work=ROOT/'workspaces/mace_omol_20260917'
        if not (self.work/'masked_context_partition_v1/result.json').exists():self.skipTest('requires actual cached partition result')

    def test_actual_endpoint_algebra_and_failed_gate(self):
        r=read_json(self.work/'masked_context_partition_v1/result.json');terms={}
        for name in ('GGR_extended','GGR_connected'):
            e=r['rows'][name]['endpoints']
            terms[name]=(e['Ca']['DFT_energy_hartree']-e['La']['DFT_energy_hartree'])*HA_TO_KCAL-(e['Ca']['masked_energy_eV']-e['La']['masked_energy_eV'])*EV_TO_KCAL
        self.assertAlmostEqual(terms['GGR_connected']-terms['GGR_extended'],r['hybrid_partition_shift_kcal_scale'],places=7)
        self.assertFalse(r['partition_gate_pass']);self.assertFalse(r['whole_inference_eligible'])
        self.assertIsNone(r['whole_context_scores']);self.assertEqual(r['new_model_calls'],0)
        self.assertTrue(all(v['maximum_coordinate_error_A']<=1e-6 for v in r['mapping_checks'].values()))

    def test_corrupted_copy_of_real_response_cannot_change_energy(self):
        r=copy.deepcopy(read_json(self.work/'masked_response_report_v1/result.json'))
        r['rows']['GGR_extended_center_Ca']['energy_eV']+=1.
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'corrupted_real_response.json';write_new(p,r)
            with self.assertRaisesRegex(InvalidArtifact,'actual numerically qualified'):
                partition(p,ROOT/'diagnostics/mace_omol_20260917/MASKED_SUBTRACTIVE_CONTEXT_PLAN.md',Path(tmp)/'unexecuted')

if __name__=='__main__':unittest.main()
