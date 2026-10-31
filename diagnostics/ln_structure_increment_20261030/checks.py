"""Saved real-record checks; no synthetic data and no refitting."""
import json,os,sys,hashlib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from curate import HERE,ROOT,W,KC,read
from structure import organic_species,potential_donors
assert os.environ.get('SLURM_JOB_ID')
checks=[]
def check(name,passed,details):
 checks.append(dict(name=name,passed=bool(passed),details=details))
rows=read(W/'primary.tsv'); byid={r['lig_id']:r for r in rows}
f=pd.read_csv(W/'features.tsv',sep='\t',dtype={'lig_id':str}).set_index('lig_id')
check('EDTA_target_sign_units',abs(float(byid['576']['y_kcal'])-(-KC*(18.30-15.36)))<1e-12,byid['576']['y_kcal'])
pd0=pd.read_csv(W/'historical/features.tsv',sep='\t',dtype={'lig_id':str})
old=pd.read_csv(ROOT/'workspaces/ln_selectivity_calibration_20261027/features.tsv',sep='\t',dtype={'lig_id':str})
check('historical_real_records_exact_replay',pd0.equals(old),len(pd0))
# Check every stored candidate's organic energy, energy closure and target-unit difference.
max_replay=0.;max_closure=0.;count=0;radial=[];stereo=[]
for row in rows:
 lid=row['lig_id']; trace=json.load(open(W/'geometry'/f'{lid}.json'))
 if trace.get('failure'): continue
 org=organic_species(row,Chem.MolFromSmiles(row['smiles']))
 AllChem.EmbedMolecule(org,randomSeed=20261030)
 prop=AllChem.MMFFGetMoleculeProperties(org,mmffVariant='MMFF94')
 ff=AllChem.MMFFGetMoleculeForceField(org,prop)
 for c in trace['candidates']:
  xyz=np.array(c['coordinates']); replay=ff.CalcEnergy(xyz.ravel().tolist())
  max_replay=max(max_replay,abs(replay-c['organic_energy_kcal']))
  max_closure=max(max_closure,abs(c['objective_kcal']-c['organic_energy_kcal']-c['restraint_energy_kcal']))
  count+=1
 for setting,sel in trace['selected'].items():
  if 'status' in sel: continue
  s=sel['Dy']['organic_energy_kcal']-sel['La']['organic_energy_kcal']
  assert abs(s-f.loc[lid,'s_'+setting])<1e-9
  for met,c in sel.items(): radial.append(dict(lig_id=lid,setting=setting,metal=met,rms_A=c['distance_rms_error_A']))
 # Real cyclic tetraacetate relative stereochemistry, coordinate triple products.
 if lid in ('597','598','599','600','601','602','603','604'):
  neutral=Chem.MolFromSmiles(row['smiles'])
  ring=next(r for r in neutral.GetRingInfo().AtomRings() if len(r) in (5,6))
  bound=[]
  for j,i in enumerate(ring):
   ns=[a.GetIdx() for a in neutral.GetAtomWithIdx(i).GetNeighbors() if a.GetSymbol()=='N']
   if ns: bound.append((j,i,ns[0]))
  expected=1 if row['name'].startswith('cis-') else -1
  vals=[]
  for c in trace['candidates']:
   xyz=np.array(c['coordinates']); signs=[]
   for j,i,n in bound:
    a=xyz[ring[(j-1)%len(ring)]]-xyz[i]
    b=xyz[ring[(j+1)%len(ring)]]-xyz[i]
    signs.append(np.dot(np.cross(a,b),xyz[n]-xyz[i]))
   vals.append(int(np.sign(np.prod(signs))))
  stereo.append(dict(lig_id=lid,name=row['name'],expected=expected,observed=sorted(set(vals)),n=len(vals),passed=all(x==expected for x in vals)))
check('all_actual_candidate_organic_energy_replay',max_replay<1e-7,dict(n=count,max_error_kcal=max_replay))
check('restraint_energy_excluded_and_objective_closes',max_closure<1e-7,dict(max_error_kcal=max_closure))
check('named_cyclic_tetraacetate_relative_stereo_preserved',all(r['passed'] for r in stereo),stereo)
# Actual organic+dummy analytic gradients versus finite differences on two real ligands.
for lid in ('3049','576'):
 row=byid[lid]; tr=json.load(open(W/'geometry'/f'{lid}.json'))
 c=next(c for c in tr['candidates'] if c['k']==20 and c['metal']=='La')
 org=organic_species(row,Chem.MolFromSmiles(row['smiles']))
 AllChem.EmbedMolecule(org,randomSeed=20261030)
 prop=AllChem.MMFFGetMoleculeProperties(org,mmffVariant='MMFF94');ff=AllChem.MMFFGetMoleculeForceField(org,prop)
 ds=np.array(tr['assignments'][c['assignment']]);target=np.array([1.216+(1.55 if org.GetAtomWithIdx(int(i)).GetSymbol()=='N' else 1.4) for i in ds])
 v=np.r_[np.array(c['coordinates']).ravel(),c['center']]
 def fun(v):
  xyz=v[:-3].reshape(-1,3);dx=xyz[ds]-v[-3:];dist=np.linalg.norm(dx,axis=1);delta=dist-target
  er=.5*20*(delta@delta);gr=20*delta[:,None]*dx/dist[:,None]
  grad=np.array(ff.CalcGrad(v[:-3].tolist())).reshape(-1,3);np.add.at(grad,ds,gr)
  return ff.CalcEnergy(v[:-3].tolist())+er,np.r_[grad.ravel(),-gr.sum(axis=0)]
 value,g=fun(v);errors=[]
 for i in [0,3*int(ds[0]),len(v)-3,len(v)-1]:
  plus=v.copy();minus=v.copy();plus[i]+=1e-5;minus[i]-=1e-5
  errors.append(abs((fun(plus)[0]-fun(minus)[0])/2e-5-g[i]))
 check('real_ligand_objective_gradient_'+lid,max(errors)<1e-4,dict(max_error=max(errors),n=len(errors)))
# Assert actual outer and inner normalization separation is also checked by evaluator.
pred=pd.read_csv(W/'predictions.tsv',sep='\t',dtype={'lig_id':str})
check('all_ligand_observations_in_one_family',all(len(set(g.family))==1 for _,g in pred.groupby('lig_id')),len(pred))
check('no_missing_descriptor_zero_substitution',set(pred[(pred.subset=='size_fit')&(pred.model=='M3')].lig_id)==set(f[f.descriptor_status=='ok'].index),int((f.descriptor_status=='ok').sum()))
pd.DataFrame(radial).to_csv(W/'selected_distance_residuals.tsv',sep='\t',index=False)
(W/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
