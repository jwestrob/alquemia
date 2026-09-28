"""Actual saved input/gradient consistency, without molecular calls."""
from pathlib import Path
import sys,json
import numpy as np
import RUN
from affordable_common import read_json,verify

def main():
 c,ledger,artifacts,p,xs,target,checks,lp=RUN.prepare();out=RUN.W/'glu91_classical_20260928';r=read_json(out/'RESULT.json');assert r['component_queries']==16 and r['QM_calls']==0
 for key,pin in r['endpoints'].items():
  d=read_json(verify(pin));metal,label=key.split('_');assert d['metal']==metal and d['configuration']==label
  assert np.array_equal(d['coordinates_A'],xs[label]);assert d['target_source_id']==target=='A/203//LA'
  assert len(d['source_ids'])==len(set(d['source_ids']))==1891
  g=np.array(d['gradient_kcal_mol_A']);assert g.shape==(1891,3) and np.isfinite(g).all()
  assert np.max(abs(g-sum(np.array(v) for v in d['component_gradients_kcal_mol_A'].values())))<1e-10
  assert abs(sum(d['components_kcal_mol'].values())-d['total_kcal_mol'])<1e-9
  assert d['quantity']=='gradient' and d['units']=='kcal/mol/angstrom'
  assert d['ledger']==r['ledger'] and d['inputs']==r['inputs']
  assert d['endpoint_state']==c['configurations'][label]['endpoints'][metal]
 for label in ['A','B']:
  la=read_json(out/f'La_{label}.json');dy=read_json(out/f'Dy_{label}.json')
  for name in ['retained_bonded','MM_LJ','MM_Coulomb']:
   assert la['components_kcal_mol'][name]==dy['components_kcal_mol'][name]
   assert np.array_equal(la['component_gradients_kcal_mol_A'][name],dy['component_gradients_kcal_mol_A'][name])
 print(json.dumps({'status':'passed','actual_configurations':4,'physical_atoms_each':1891,'source_cap_state_field_checks':checks,'molecular_calls':0},indent=2))
if __name__=='__main__':main()
