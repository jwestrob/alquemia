"""Compare actual repaired origin outputs; no molecular evaluations."""
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
def read(p):return json.loads(Path(p).read_text())
def main():
 old=ROOT/'workspaces/metal_environment_response_20260926/coupled_scaffold_classical_v1'
 new=ROOT/'workspaces/metal_environment_response_20260926/scaffold_origin_DQ2_v1'
 r=read(new/'RESULT.json');oldr=read(old/'RESULT.json');d=new/'Hans8DQ2'
 assert r['component_queries']==8
 # Exact frozen physical representation, FF parameters and retained terms.
 for name in ['exceptions.json','bonded_terms.json','cross_term_decisions.json']:
  assert read(d/name)==read(old/name)
 a,b=read(d/'particles.json'),read(old/'particles.json')
 for x,y in zip(a,b):
  assert {k:v for k,v in x.items() if k!='representation'}=={k:v for k,v in y.items() if k!='representation'}
 assert len(a)==len(b)==1891
 assert read(d/'LEDGER.json')['native_system']==read(old/'LEDGER.json')['native_system']
 residuals={}
 for m in ['La','Dy']:
  n=read(d/f'{m}_A.json');o=oldr['results'][m]['configurations']['origin']
  residuals[m]={k:n['components_kcal_mol'][k]-v for k,v in o['components_kcal_mol'].items()}
  assert max(map(abs,residuals[m].values()))<1e-9
  assert n['target_source_id']=='A/203//LA'
  assert n['source_ids'][-1]=='A/203//LA'
  assert np.array_equal(n['coordinates_A'],[v['xyz_A'] for v in b])
  g=np.array(n['gradient_kcal_mol_A']);assert g.shape==(1891,3) and np.isfinite(g).all()
  assert np.max(abs(g-sum(np.array(v) for v in n['component_gradients_kcal_mol_A'].values())))<1e-10
 print(json.dumps({'status':'passed','component_energy_residuals_kcal_mol':residuals,'native_and_physical_ledger_identical':True,'actual_gradient_atoms':1891,'quantum_calls':0},indent=2))
if __name__=='__main__':main()
