"""Executed-artifact tests; no molecular evaluations or invented energies."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
W=ROOT/'workspaces/metal_environment_response_20260926/scaffold_origin_transfer_v1'
def read(p):return json.loads(Path(p).read_text())
def test():
 result=read(W/'RESULT.json');assert result['component_queries']==16
 counts=[]
 for source,expected_spectator,nativeN,physicalN,fieldQ in [('Hans8FNR','Dy',2202,2206,9),('Mex8FNS','Nd',2045,2049,2)]:
  d=W/source;ledger=read(d/'LEDGER.json');p=read(d/'particles.json');c=read(ledger['inputs']['path']);native=read(read(c['export']['path'])['atoms']['path']);real=read(c['physical_source']['path']);env=read(c['environment_atoms']['path']);selected=set(read(c['boundary_mapping']['path'])['all_selected_QM_source_ids'])
  assert len(native)==nativeN and len(p)==physicalN
  assert [v['id'] for v in native]==[v['id'] for v in p[:nativeN]]
  assert np.array_equal([v['xyz_A'] for v in real],[v['xyz_A'] for v in p[:nativeN]])
  assert {v['id'] for v in p if v['region']=='QM'}==selected|{'target_EF3'}
  assert [v['element'] for v in p[nativeN:-1]]==[expected_spectator]*3
  assert abs(sum(v['charge_e'] for v in p)-fieldQ)<1e-9
  assert abs(sum(v['charge_e'] for v in env)-fieldQ)<1e-9
  assert ledger['counts']['HarmonicBondForce:cross_radial_stiffness']==4
  la,dy=[read(d/f'{m}_A.json') for m in ['La','Dy']]
  assert la['source_ids']==dy['source_ids'] and len(set(la['source_ids']))==physicalN
  assert la['target_source_id']==dy['target_source_id'] and la['target_source_id']!='target_EF3'
  assert np.array_equal(la['coordinates_A'],dy['coordinates_A'])
  for a in [la,dy]:
   assert not a['QM_MM_Coulomb_included'] and not a['C4_induction_included']
   assert a['quantity']=='gradient' and a['units']=='kcal/mol/angstrom'
   grad=np.array(a['gradient_kcal_mol_A']);assert grad.shape==(physicalN,3) and np.isfinite(grad).all()
   assert np.max(abs(grad-sum(np.array(g) for g in a['component_gradients_kcal_mol_A'].values())))<1e-10 # JSON sorts group keys; sum ordering roundoff
   assert abs(a['total_kcal_mol']-sum(a['components_kcal_mol'].values()))<1e-9
  for name in ['retained_bonded','MM_LJ','MM_Coulomb']:
   assert la['components_kcal_mol'][name]==dy['components_kcal_mol'][name]
   assert np.array_equal(la['component_gradients_kcal_mol_A'][name],dy['component_gradients_kcal_mol_A'][name])
  counts.append({'source':source,'physical_atoms':physicalN,'passed':True})
 print(json.dumps({'status':'passed','real_artifact_tests':counts,'molecular_calls':0},indent=2))
if __name__=='__main__':test()
