"""Health parsing on retained actual failed-Dy/successful-La ORCA outputs only."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from slurm_native_health_watch import inspect_text,snapshot,POLICY
W=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_hans_scout_v1'

class RealSCFHealth(unittest.TestCase):
    def test_actual_dy_pathologies(self):
        for name in ('Dy_A','Dy_B'):
            with self.subTest(endpoint=name):
                result=inspect_text((W/name/'endpoint.out').read_text())
                self.assertIn('repeated_extreme_trah_micro_steps',result['alerts'])
                self.assertFalse(result['normal_termination'])
                self.assertGreaterEqual(len(result['extreme_micro_distinct_macros']),3)
        result=inspect_text((W/'Dy_B/endpoint.out').read_text())
        self.assertIn('persistent_large_trah_error_without_twofold_improvement',result['alerts'])

    def test_actual_completed_la_not_mislabeled(self):
        for name in ('La_A','La_B'):
            result=inspect_text((W/name/'endpoint.out').read_text())
            self.assertTrue(result['normal_termination'])
            self.assertTrue(result['scf_converged'])
            self.assertEqual(result['alerts'],[])

    def test_real_partial_execution_resource_accounting(self):
        result=snapshot(W/'manifest.json',224)
        self.assertEqual(result['terminal_endpoints'],2)
        self.assertEqual(result['maximum_remaining_endpoint_rank_demand'],112)
        self.assertEqual(result['allocated_cpus_without_remaining_endpoint_work'],112)
        self.assertIn('not measured',result['resource_note'])

    def test_live_prefix_can_alert_before_terminal(self):
        # Prefix of a real retained log: exactly what was available at that point.
        lines=(W/'Dy_B/endpoint.out').read_text().splitlines()
        stop=next(i for i,line in enumerate(lines) if '(TRAH MAcro)' in line and line.split()[0]=='20')
        result=inspect_text('\n'.join(lines[:stop]))
        self.assertFalse(result['normal_termination'])
        self.assertTrue(result['alerts'])
        self.assertLess(result['trah_macro_count'],21)
        self.assertFalse(POLICY['automatic_cancellation'])

if __name__=='__main__':unittest.main()
