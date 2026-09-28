"""Transfer reviewed real H repair to fixed full regions without molecular calls."""
from pathlib import Path
import sys,json,copy,math
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,xyz
from affordable_response import cap_jacobians
from metal_environment_lanm_reference import check_config

def prepare(source,base):
 olddir=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1'/source;c=check_config(read_json(olddir/'INPUTS.json'));adpath=ROOT/'workspaces/metal_environment_response_20260926/scaffold_H_repair_v1'/source/'REVIEWED_ADMISSION.json';ad=read_json(adpath);assert ad['admitted'] and ad['inputs']==record(olddir/'INPUTS.json');export=read_json(verify(ad['parent_export']));original=read_json(verify(export['atoms']));coords=read_json(verify(ad['saved_coordinates']));assert coords['source_ids']==[v['id'] for v in original];free=set(ad['free_H_ids']);full=copy.deepcopy(original);checks={};oldby={v['id']:v for v in original}
 for a,pos in zip(full,coords['xyz_A']):
  if a['id'] not in free:assert a['xyz_A']==pos
  else:assert a['element']=='H'
  a['xyz_A']=pos
 A={v['id']:np.array(v['xyz_A']) for v in full};B={sid:v.copy() for sid,v in A.items()};pert=c['perturbation'];ca,cb=[A[sid] for sid in pert['axis']];axis=(cb-ca)/np.linalg.norm(cb-ca);angle=math.radians(pert['right_hand_angle_degrees']);moving=pert['rotate_source_atoms']
 for sid in moving:
  v=A[sid]-cb;B[sid]=cb+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+axis*np.dot(axis,v)*(1-np.cos(angle))
 bonds=read_json(verify(export['bonds']));bondres=max(abs(np.linalg.norm(A[u]-A[v])-np.linalg.norm(B[u]-B[v])) for u,v in bonds);assert bondres<1e-12
 d=base/source;d.mkdir(parents=True,exist_ok=False);write_new(d/'physical_source_A.json',full);physB=[dict(v,xyz_A=B[v['id']].tolist()) for v in full];write_new(d/'physical_source_B.json',physB);oldenv=read_json(verify(c['environment_atoms']));env=copy.deepcopy(oldenv)
 for a in env:
  if a['id'] in A:a['xyz_A']=A[a['id']].tolist();assert np.array_equal(A[a['id']],B[a['id']])
  else:assert a['kind']=='spectator_formal_monopole'
 write_new(d/'environment_atoms.json',env);pc=d/'environment.pc';pc.write_text(str(len(env))+'\n'+''.join(' '.join(format(x,'.17g') for x in (v['charge_e'],*v['xyz_A']))+'\n' for v in env));boundary=read_json(verify(c['boundary_mapping']));boundary['original_boundary_mapping']=c['boundary_mapping'];boundary['coordinate_preparation']='reviewed_C_bound_H_v1'
 for row in boundary['ledger']:
  row['original_dipole_diagnostic']={k:row[k] for k in ('dipole_before_eA','retained_MM_dipole_after_eA')}
  prefix=row['partial_residue'];local=[v for v in full if v['id'].startswith(prefix)];remain=[v for v in env if v['id'].startswith(prefix)]
  row['dipole_before_eA']=sum((v['charge_e']*np.array(v['xyz_A']) for v in local),np.zeros(3)).tolist();row['retained_MM_dipole_after_eA']=sum((v['charge_e']*np.array(v['xyz_A']) for v in remain),np.zeros(3)).tolist()
 write_new(d/'boundary_mapping.json',boundary)
 maps={};configs={}
 for label,physical in [('A',A),('B',B)]:
  oldmap=read_json(verify(c['core_mapping'][label]));mapping=copy.deepcopy(oldmap)
  for a in mapping:
   if a['kind']=='source':a['xyz_A']=physical[a['source_id']].tolist()
   elif a['kind']=='cap':
    u,v=physical[a['retained_source_id']],physical[a['omitted_source_id']];length=a['length_A'];ja,jb=cap_jacobians(u,v,length);a.update(xyz_A=(u+length*(v-u)/np.linalg.norm(v-u)).tolist(),retained_xyz_A=u.tolist(),omitted_xyz_A=v.tolist(),jacobian_retained=ja.tolist(),jacobian_omitted=jb.tolist());assert a['xyz_A']==oldmap[a['qm_index']]['xyz_A']
  write_new(d/f'core_mapping_{label}.json',mapping);maps[label]=record(d/f'core_mapping_{label}.json');endpoints={}
  for metal,e0 in c['configurations'][label]['endpoints'].items():
   oldxyz=xyz(verify(e0['xyz']));path=d/f'{metal}_{label}.xyz';path.write_text(str(len(mapping))+'\n'+f'{source} reviewed C-bound-H, {label}, target{metal}; new preparation identity\n'+''.join(a[0]+' '+' '.join(format(t,'.17g') for t in m['xyz_A'])+'\n' for a,m in zip(oldxyz,mapping)));endpoints[metal]=dict(e0,xyz=record(path));assert [a[0] for a in xyz(path)]==[a[0] for a in oldxyz]
  configs[label]={'endpoints':endpoints,'pointcharges':record(pc),'core_mapping':maps[label]}
 new=copy.deepcopy(c);new.update(protocol_id=c['protocol_id']+'_CboundH_repaired_v2',original_inputs=record(olddir/'INPUTS.json'),source_state=record(adpath),reviewed_admission=record(adpath),physical_source=record(d/'physical_source_A.json'),physical_sources={'A':record(d/'physical_source_A.json'),'B':record(d/'physical_source_B.json')},core_mapping=maps,boundary_mapping=record(d/'boundary_mapping.json'),environment_atoms=record(d/'environment_atoms.json'),configurations=configs,implementation=record(__file__),preparation_plan=record(Path(__file__).with_name('PLAN.md')),energy_force_calls=0,coordinate_preparation='reviewed_C_bound_H_v1');write_new(d/'INPUTS.json',new);check_config(new)
 # Direct physical maps: source QM+field cover all physical atoms except four explicit zero-charge CAs.
 for label,physical in [('A',A),('B',B)]:
  for a in read_json(verify(maps[label])):
   if a['kind']=='source':assert np.array_equal(a['xyz_A'],physical[a['source_id']])
  for a in env:
   if a['id'] in physical:assert np.array_equal(a['xyz_A'],physical[a['id']])
 assert [v['charge_e'] for v in env]==[v['charge_e'] for v in oldenv];assert new['charge_accounting']==c['charge_accounting'];assert new['waters']==c['waters'];assert new['spectator_occupancy']==c['spectator_occupancy'];assert [v['id'] for v in env]==[v['id'] for v in oldenv]
 changed=[v['id'] for v in full if v['xyz_A']!=oldby[v['id']]['xyz_A']];assert set(changed)<=free
 result={'status':'passed','source_id':source,'inputs':record(d/'INPUTS.json'),'original_inputs':record(olddir/'INPUTS.json'),'reviewed_admission':record(adpath),'changed_real_H_atoms':len(changed),'max_B_bond_length_change_A':bondres,'checks':['exact unchanged heavy/water/exchangeableH coordinates and identities','only reviewed source C-bound H transferred','same original region charge inventory and spectators','B exact original torsion rule from repaired A','physical source matches all real QM and MM field coordinates','caps/Jacobians regenerated from unchanged source cuts','paired target XYZs and electronic states unchanged except intended species','A/B identical field values and ordering','native input contract and recursive pins pass'],'molecular_calls':0};write_new(d/'CHECK.json',result);return result
if __name__=='__main__':
 base=Path(sys.argv[1]).resolve();base.mkdir(parents=True,exist_ok=False);results=[prepare(source,base) for source in ['Hans8DQ2','Hans8FNR','Mex8FNS']];write_new(base/'CHECKS.json',results);print(json.dumps(results,indent=2))
