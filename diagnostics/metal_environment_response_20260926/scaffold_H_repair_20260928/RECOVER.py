"""Recover actual saved geometry/force diagnostics after numpy-bool serialization failure."""
from pathlib import Path
import json,sys
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
base=ROOT/'workspaces/metal_environment_response_20260926/scaffold_H_repair_v1'
rows=[]
for source in ['Hans8DQ2','Hans8FNR','Mex8FNS']:
 out=base/source;p=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1'/source;c=read_json(p/'INPUTS.json');ex=read_json(verify(c['export']));a=read_json(verify(ex['atoms']));ids={v['id']:i for i,v in enumerate(a)};adj={i:set() for i in range(len(a))}
 for u,v in read_json(verify(ex['bonds'])):i,j=ids[u],ids[v];adj[i].add(j);adj[j].add(i)
 free={i for i,v in enumerate(a) if v['element']=='H' and len(adj[i])==1 and a[next(iter(adj[i]))]['element']=='C'};fixed=sorted(set(range(len(a)))-free);xx=read_json(out/'coordinates_A.json');assert xx['source_ids']==[v['id'] for v in a];x=np.array(xx['xyz_A']);x0=np.array([v['xyz_A'] for v in a]);forces=read_json(out/'forces.json');ga=-np.array(forces['after_kcal_mol_A']);gb=-np.array(forces['before_kcal_mol_A']);system=mm.XmlSerializer.deserialize((out/'system_objective.xml').read_text());bf=next(v for v in system.getForces() if isinstance(v,mm.HarmonicBondForce));bondmax=0.
 for i in range(bf.getNumBonds()):
  u,v,r,k=bf.getBondParameters(i);u,v=int(u),int(v)
  if u in free or v in free:bondmax=max(bondmax,float(abs(np.linalg.norm(x[u]-x[v])-r.value_in_unit(unit.angstrom))))
 clashes=[];mind=float('inf')
 for i in sorted(free):
  exc={i}|adj[i]|set().union(*(adj[j] for j in adj[i]));other=np.array(sorted(set(range(len(a)))-exc));d=np.linalg.norm(x[other]-x[i],axis=1);mind=min(mind,float(d.min()))
  for j in other[d<.65]:
   if np.linalg.norm(x0[j]-x0[i])>=.65:clashes.append([a[i]['id'],a[int(j)]['id']])
 inversions=[]
 for i,nei in adj.items():
  if len(nei)!=4 or a[i]['element']=='H' or not(nei&free):continue
  nn=sorted(nei);v0=np.linalg.det(x0[nn[:3]]-x0[nn[3]]);v=np.linalg.det(x[nn[:3]]-x[nn[3]])
  if abs(v0)>1e-6 and v*v0<0:inversions.append({'id':a[i]['id'],'neighbors':[a[j]['id'] for j in nn],'hydrogen_neighbors':sum(a[j]['element']=='H' for j in nn)})
 partial=(out/'RESULT.json').read_text();assert '"checks": {' in partial
 actual=json.loads(partial.split(',\n  "checks": {')[0]+'\n}');maxgrad=float(np.max(abs(ga[sorted(free)])));maxmove=float(np.linalg.norm(x-x0,axis=1).max());checks={'optimizer_success':None,'freeH_gradient':maxgrad<=.1,'fixed_exact':bool(np.array_equal(x[fixed],x0[fixed])),'C_H_lengths':bondmax<=.06,'H_angles_improve':actual['after_angles']['H_angle_energy_kcal_mol']<actual['before_angles']['H_angle_energy_kcal_mol'],'H_angles_bounded':actual['after_angles']['max_freeH_angle_deviation_deg']<=25,'displacement':maxmove<=2.2,'no_new_clashes':not clashes,'no_inversions':not inversions}
 r={'source_id':source,'status':'geometry_recovered_original_admission_unavailable','inputs':record(p/'INPUTS.json'),'saved_coordinates':record(out/'coordinates_A.json'),'saved_forces':record(out/'forces.json'),'partial_result':record(out/'RESULT.json'),'actual_preserved_evidence':actual,'checks':checks,'optimizer_status':'not_saved_before_serialization_error','max_freeH_gradient_kcal_mol_A':maxgrad,'before_max_freeH_gradient_kcal_mol_A':float(np.max(abs(gb[sorted(free)]))),'max_displacement_A':maxmove,'max_CH_length_error_A':bondmax,'minimum_nonbonded_freeH_distance_A':mind,'new_clashes':clashes,'inversions':inversions,'free_H_ids':[a[i]['id'] for i in sorted(free)],'implementation':record(__file__),'new_molecular_calls':0};write_new(out/'RECOVERED.json',r);rows.append({'source':source,'checks':checks,'maxgrad':maxgrad,'inversions':inversions});print(source,json.dumps(rows[-1]))
write_new(base/'RECOVERY_COLLECTION.json',rows)
