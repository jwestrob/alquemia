import sys,copy,json,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,record,write_new
from metal_environment_force_assembly import prepare_mapping
from metal_environment_lanm_reference import check_config,torsion_tangent
src=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Hans8DQ2/INPUTS.json'
out=ROOT/'workspaces/metal_environment_response_20260926/glu91_motion_v1';out.mkdir(exist_ok=False)
c=read_json(src);mp=prepare_mapping(src,'A','La');atoms=mp['physical_atoms'];ids=[a['id'] for a in atoms];ix={s:i for i,s in enumerate(ids)};x=np.array([a['xyz_A'] for a in atoms]);ex=read_json(verify(c['export']));bonds=read_json(verify(ex['bonds']));adj={s:set() for s in ids}
for a,b in bonds:adj[a].add(b);adj[b].add(a)
a,b='A/91//CB','A/91//CG';adj[a].remove(b);adj[b].remove(a);moving={b};todo=[b]
while todo:
 k=todo.pop()
 for n in adj[k]-moving:moving.add(n);todo.append(n)
assert a not in moving
adj[a].add(b);adj[b].add(a)
axis=x[ix[b]]-x[ix[a]];axis/=np.linalg.norm(axis);pivot=x[ix[b]];mi=[ix[s] for s in sorted(moving)];moving_nonzero=[s for s in sorted(moving) if np.linalg.norm(np.cross(axis,x[ix[s]]-pivot))>1e-12]
assert set(moving)<={m.get('source_id') for m in mp['core']}
assert not moving&{m['id'] for m in mp['field']}
for cap in mp['caps']:assert cap['retained_id'] not in moving and cap['omitted_id'] not in moving
cc=copy.deepcopy(c);cc['protocol_id']='nikasha_Hans8DQ2_EF3_CboundH_repaired_Glu91_pm1_v1';cc['origin_inputs']=record(src);cc['origin_configuration']='A';cc['implementation']=record(__file__);cc['preparation_plan']=record(Path(__file__).with_name('PLAN.md'));cc['perturbation']=dict(mode_id='Glu91_CB_CG_source_graph_pm1_v1',axis=[a,b],axis_unit=axis.tolist(),pivot_xyz_A=pivot.tolist(),rotate_source_atoms=moving_nonzero,changed_source_ids=moving_nonzero,right_hand_angle_degrees=2.,origin_angles_degrees={'A':-1.,'B':1.},type='whole graph carboxylate torsion',selection='development selected using original-H differential forces');cc['configurations']={};cc['physical_sources']={};cc['core_mapping']={}
basephysical=read_json(verify(c['physical_sources']['A']));tangents={};checks={}
def tangent(y):
 t=np.zeros_like(y);t[mi]=np.cross(axis,y[mi]-pivot);norm=float(np.linalg.norm(t));return dict(raw_physical_A_per_radian=t.tolist(),physical_unit_tangent=(t/norm).tolist(),physical_norm_A_per_radian=norm)
tangents['origin']=tangent(x)
for label,angle in [('A',-1.),('B',1.)]:
 y=x.copy();v=x[mi]-pivot;t=math.radians(angle);y[mi]=pivot+v*math.cos(t)+np.cross(axis,v)*math.sin(t)+np.outer(v@axis,axis)*(1-math.cos(t))
 bond_error=max(abs(np.linalg.norm(y[ix[u]]-y[ix[v]])-np.linalg.norm(x[ix[u]]-x[ix[v]])) for u,v in bonds)
 volerr=0.
 for sid,neighbors in adj.items():
  if len(neighbors)!=4:continue
  n=sorted(neighbors);v0=np.linalg.det(np.array([x[ix[z]]-x[ix[n[3]]] for z in n[:3]]));v1=np.linalg.det(np.array([y[ix[z]]-y[ix[n[3]]] for z in n[:3]]));volerr=max(volerr,abs(v1-v0));assert v0*v1>=0
 assert bond_error<1e-10 and volerr<1e-9
 mapping=copy.deepcopy(mp['core'])
 for m in mapping:
  if m.get('source_id') in moving:m['xyz_A']=y[ix[m['source_id']]].tolist()
 p=out/f'core_mapping_{label}.json';write_new(p,mapping);conf=copy.deepcopy(c['configurations']['A']);conf.update(core_mapping=record(p),geometry_id='Glu91_'+('minus' if angle<0 else 'plus')+'1deg',angle_from_origin_degrees=angle)
 for metal in ('La','Dy'):
  old=Path(verify(c['configurations']['A']['endpoints'][metal]['xyz'])).read_text().splitlines();lines=[old[0],conf['geometry_id']]
  for m,line in zip(mapping,old[2:]):lines.append(line.split()[0]+' '+' '.join(f'{z:.12f}' for z in m['xyz_A']))
  p=out/f'{metal}_{label}.xyz';p.write_text('\n'.join(lines)+'\n');conf['endpoints'][metal]['xyz']=record(p)
 cc['configurations'][label]=conf;cc['core_mapping'][label]=conf['core_mapping'];physical=copy.deepcopy(basephysical)
 for m in physical:
  if m['id'] in moving:m['xyz_A']=y[ix[m['id']]].tolist()
 p=out/f'physical_source_{label}.json';write_new(p,physical);cc['physical_sources'][label]=record(p)
 tangents[label]=tangent(y);tangents[label]['core_A_per_radian']=torsion_tangent(cc,label).tolist();checks[label]=dict(max_all_source_bond_change_A=bond_error,max_tetrahedral_signed_volume_change_A3=volerr,max_displacement_A=float(np.max(np.linalg.norm(y-x,axis=1))),changed_source_ids=moving_nonzero,field_unchanged=True,caps_unchanged=True)
cc['physical_source']=cc['physical_sources']['A'];cc['energy_force_calls']=0;cc['source_state_sensitivity']='modeled development displacements; no populations or affinity';p=out/'TANGENTS.json';write_new(p,dict(origin_inputs=record(src),physical_source_ids=ids,moving_source_ids=sorted(moving),states=tangents));cc['tangents']=record(p);check_config(cc);p=out/'INPUTS.json';write_new(p,cc)
for label in ('A','B'):
 for metal in ('La','Dy'):prepare_mapping(p,label,metal)
write_new(out/'CHECKS.json',dict(checks=checks,paired_source_maps_and_cap_jacobians=True,inputs=record(p),energy_force_calls=0));print(p)
