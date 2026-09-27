"""Matched real saved components; no fabricated energies or molecular calls."""
from pathlib import Path
import sys,unittest,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,InvalidArtifact,HA_TO_KCAL,BOHR_TO_A
from metal_environment_reference import read_pcgrad
from metal_environment_combine_saved import matches,physical_gradient,combine
B=ROOT/'workspaces/metal_environment_response_20260926'
class CombinedTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  paths=[B/'component_checks_v1/FINAL_COLLECTION_1219497.json',B/'reference_scout_v1/FINAL_COLLECTION_1219207.json',B/'force_checks_v1/COLLECTION_1219319.json']
  if not all(p.exists() for p in paths):raise unittest.SkipTest('Real completed sources unavailable')
  cls.paths=paths;cls.cm=read_json(B/'component_checks_v1/manifest.json');l=read_json(verify(cls.cm['ledger']));cls.particles=read_json(verify(l['artifacts']['particles.json']));cls.caps=read_json(verify(l['artifacts']['caps.json']));c=read_json(verify(l['inputs']['source_manifest']));cls.env=read_json(verify(c['environments']['A']['atoms']));cls.source=read_json(paths[1]);sm=read_json(verify(cls.source['manifest']));cls.qt=next(t for t in sm['tasks'] if t['task_id']=='Ca_A');cls.x=np.load(verify(next(t for t in cls.cm['tasks'] if t['task_id']=='Ca_A')['coordinates']))
 def test_actual_mapping_and_mismatched_B_rejected(self):
  self.assertLess(matches('A',self.x,self.qt,self.particles,self.caps,self.env),1e-11)
  wrong=np.load(verify(next(t for t in self.cm['tasks'] if t['task_id']=='Ca_B')['coordinates']))
  with self.assertRaisesRegex(InvalidArtifact,'unmatched'):matches('A',wrong,self.qt,self.particles,self.caps,self.env)
 def test_cap_chain_rule_preserves_total_gradient(self):
  r=self.source['rows']['Ca_A'];q=np.array(r['gradient_kcal_mol_per_A']);p=read_pcgrad(verify(r['pointcharge_gradient']),9087)*HA_TO_KCAL/BOHR_TO_A;g=physical_gradient(self.x,q,p,self.particles,self.caps,self.env);np.testing.assert_allclose(g.sum(axis=0),q.sum(axis=0)+p.sum(axis=0),atol=1e-10,rtol=0)
 def test_saved_combination_keeps_failures_and_missing_cells(self):
  r=combine(*self.paths);self.assertEqual(r['matched_configurations'],24);self.assertEqual(r['native_qualification_status'],'not_qualified');self.assertFalse(r['full_hybrid_qualified']);self.assertFalse(r['solvent_present']);self.assertEqual(r['La_minus_Ca_response_kcal_mol']['combined'],r['La_minus_Ca_response_kcal_mol']['electronic']);self.assertTrue(all(v['combined_finite_difference'] is None for v in r['unavailable_comparisons']));self.assertEqual(len(r['directional_FD_residuals']),8)
if __name__=='__main__':unittest.main()
