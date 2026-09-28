"""Health parsing on retained actual failed-Dy/successful-La ORCA outputs only."""
import sys
import re
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

    def test_gross_diis_alert_before_trah_on_actual_prefixes(self):
        for name in ('Dy_A','Dy_B','La_A','La_B'):
            with self.subTest(endpoint=name):
                text=(W/name/'endpoint.out').read_text()
                start=text.index('Iteration    Energy')
                next_iteration=re.search(r'^\s*17\s+-',text[start:],re.M)
                self.assertIsNotNone(next_iteration)
                prefix=text[:start+next_iteration.start()]
                result=inspect_text(prefix)
                self.assertEqual(result['last_diis_iteration'],16)
                self.assertTrue(result['diis_current_phase'])
                self.assertEqual('gross_diis_nonprogress' in result['alerts'],name.startswith('Dy'))

    def test_previous_diis_rows_do_not_trigger_after_trah_transition(self):
        result=inspect_text((W/'Dy_A/endpoint.out').read_text())
        self.assertFalse(result['diis_current_phase'])
        self.assertNotIn('gross_diis_nonprogress',result['alerts'])

    def test_actual_small_dy_near_root_stagnation(self):
        small=ROOT/'workspaces/metal_environment_response_20260926/dy_small_guess_v3'
        for name in ('PModel','HCore'):
            text=(small/name/'endpoint.out').read_text()
            # Actual printed molecular log before native error/footer, preserving
            # all real convergence rows rather than inventing residual values.
            prefix=text.split('ORCA finished by error termination')[0]
            result=inspect_text(prefix)
            self.assertEqual(result['printed_trah_gradient_tolerance'],1e-5)
            self.assertIn('near_root_stagnation_review_only',result['alerts'])
            self.assertLess(result['near_root_energy_span_hartree'],1e-6)
            self.assertEqual(len(result['near_root_macro_window']),12)
            self.assertTrue(all(row['error_norm']>5e-5 for row in result['near_root_macro_window']))
        hcore=inspect_text((small/'HCore/endpoint.out').read_text())
        self.assertEqual([r['iteration'] for r in hcore['near_root_macro_window']],list(range(81,93)))
        self.assertTrue(any(r['solver']=='NR' for r in hcore['near_root_macro_window']))

    def test_near_root_rule_uses_actual_tolerance_and_never_invents_one(self):
        text=(ROOT/'workspaces/metal_environment_response_20260926/dy_small_guess_v3/HCore/endpoint.out').read_text()
        # Explicitly corrupted copy removes the printed threshold: unavailable
        # evidence must not acquire a guessed tolerance or trigger this rule.
        corrupt=re.sub(r'^.*Converg\. threshold.*$', '', text, flags=re.M)
        self.assertNotEqual(corrupt,text)
        result=inspect_text(corrupt)
        self.assertIsNone(result['printed_trah_gradient_tolerance'])
        self.assertNotIn('near_root_stagnation_review_only',result['alerts'])

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
