"""Real pinned state/manifest checks; no fabricated ORCA success output."""
import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from metal_environment_frozen_f import state,validate,input_text
from metal_environment_electronic_state import describe
class FrozenFTests(unittest.TestCase):
 def test_real_state_distinction(self):
  p=ROOT/'workspaces/lanm_series_followup_20260923/prepared_v1/endpoints/Hans_EF3__Dy/core.xyz'
  old=describe(p,'Dy',-1,6);new=state(p)
  self.assertEqual(old['explicit_electrons'],235);self.assertEqual(new['explicit_electrons'],208)
  self.assertEqual(new['physical_multiplicity'],6);self.assertEqual(new['effective_multiplicity'],1)
  self.assertEqual(old['ecp_core_electrons'],28)
 def test_actual_manifest(self):
  self.assertEqual(validate(ROOT/'workspaces/metal_environment_response_20260926/dy_frozen_f_scout_v1/manifest.json')['status'],'dry_run_pass')
 def test_real_displacement_axis(self):
  import numpy as np
  from metal_environment_frozen_f_response import direction,STEPS
  from affordable_common import xyz
  p=ROOT/'workspaces/lanm_series_followup_20260923/prepared_v1/endpoints/Hans_EF3__Dy/core.xyz'
  d=direction(p);a=xyz(p);q=np.array([x[1:] for x in a]);v=np.array(d['unit_vector'])
  self.assertAlmostEqual(np.linalg.norm(v),1)
  self.assertEqual(a[d['metal_index']][0],'Dy');self.assertEqual(a[d['donor_index']][0],'O')
  for step in STEPS.values():
   shifted=q.copy();shifted[d['metal_index']]+=step*v
   self.assertAlmostEqual(np.linalg.norm(shifted[d['donor_index']]-shifted[d['metal_index']]),np.linalg.norm(q[d['donor_index']]-q[d['metal_index']])-step)
 def test_actual_custom_ecp_header(self):
  import re
  p=ROOT/'workspaces/metal_environment_response_20260926/dy_frozen_f_scout_v1/Dy_origin/endpoint.out'
  hits=re.findall(r'Type\s+(\w+)\s+ECP(?:\s+(\S+))?\s+\(replacing\s+(\d+)\s+core electrons',p.read_text())
  self.assertEqual(hits,[('Dy','','55')])
 def test_executed_origin_gradient(self):
  from metal_environment_frozen_f import collect
  r=collect(ROOT/'workspaces/metal_environment_response_20260926/dy_frozen_f_scout_v1/manifest.json')['rows']['Dy_origin']
  self.assertEqual(r['status'],'complete');self.assertEqual(len(r['gradient_kcal_mol_per_A']),50)
  self.assertFalse(r['force_consistency_qualified'])
if __name__=='__main__':unittest.main()
