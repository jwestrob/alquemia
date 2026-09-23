"""Pinned real-context numerical-policy replay; no fabricated chemistry."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import slsqp_precision as s
import adaptive_completion as original
import union_adaptive
import nikasha_pool
M=ROOT/'workspaces/slsqp_precision_20260923/proposals_v1/manifest.json'
class Precision(unittest.TestCase):
 def test_exact_four_sources_eight_origins_and_selector(self):
  v=s.validate(M);self.assertEqual((v['cases'],v['optimizer_starts'],v['new_q0_calls']),(4,8,0))
 def test_optimizer_is_private_and_only_ftol_changes(self):
  before=dict(original.SETTINGS);uv=union_adaptive.validate;pv=nikasha_pool.validate
  u,p=s.engine()
  self.assertIsNot(u,union_adaptive);self.assertIsNot(u.completion,original);self.assertIsNot(p,nikasha_pool)
  self.assertEqual(original.SETTINGS,before);self.assertIs(union_adaptive.validate,uv);self.assertIs(nikasha_pool.validate,pv)
  self.assertEqual({k for k in before if before[k]!=u.completion.SETTINGS[k]},{'optimizer_ftol'})
  self.assertEqual(u.completion.SETTINGS['optimizer_ftol'],1e-8)
 def test_real_pathological_source_retained(self):
  m=s.read_json(M);src=m['sources'][0];self.assertEqual(src['case_id'],s.CASES[0])
  r=s.read_json(s.verify(src['parent']).parent/'proposals'/(s.CASES[0]+'__La')/'result.json')
  self.assertEqual(r['status'],'proposal_available');self.assertEqual(r['optimizer']['function_evaluations'],1304)
  self.assertLessEqual(r['final_geometry']['maximum_heavy_displacement_A'],.8+1e-7)
 def test_actual_comparison_complete_when_run(self):
  p=M.parent.parent/'COMPARISON_v1.json'
  if not p.exists():self.skipTest('approved molecular precision comparison not yet executed')
  r=s.read_json(p);self.assertEqual([v['case_id'] for v in r['rows']],list(s.CASES));self.assertFalse(r['calibration_changed'])
  for c in r['rows']:
   if c['new_R'] is not None:self.assertAlmostEqual(c['delta_R'],c['new_R']-c['old_R'],places=10)
   for z,v in c['endpoints'].items():
    self.assertEqual(v['native_energy_pass'],abs(v['native_proposal_delta_kcal_mol'])<=.001)
if __name__=='__main__':unittest.main()
