"""Actual completed GGR robustness results, including the failed direction gate."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from mace_omol_ggr_robustness import report


class RealGGRRobustness(unittest.TestCase):
    def test_actual_failure_is_retained_despite_numerical_success(self):
        w = ROOT/'workspaces/mace_omol_20260917'
        if not (w/'ggr_structure_report_v1/result.json').exists():
            self.skipTest('requires actual completed four-forward GGR experiment')
        with tempfile.TemporaryDirectory() as tmp:
            r = report(w/'charge_ablation_report_v1/result.json',
                       w/'ggr_masked_2fw0_report_v1/result.json',
                       w/'ggr_masked_2fvy_report_v1/result.json',
                       ROOT/'diagnostics/mace_omol_20260917/GGR_STRUCTURE_ROBUSTNESS_PLAN.md',
                       ROOT/'diagnostics/mace_omol_20260917/GGR_SAVED_READOUT_PLAN.md',
                       w/'ggr_structure_sacct_v1.tsv', Path(tmp)/'report')
        self.assertTrue(r['numerical_gate_pass'])
        self.assertFalse(r['robustness_gate_pass'])
        self.assertEqual((r['passed_margins'], r['total_margins']), (2, 6))
        self.assertEqual(r['biological_groups'], 2)
        self.assertEqual(r['cost']['new_successful_forwards'], 4)
        self.assertAlmostEqual(r['GGR_range_model_kcal'], 21.988078952339475, places=8)
        self.assertFalse(r['broad_affinity_validated'])


if __name__ == '__main__': unittest.main()
