"""Source-defined physical load projections; analysis only, no molecular calls."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,InvalidArtifact
from metal_environment_force_assembly import prepare_mapping,target_source_id
AXES={'Asp85_CA_CB':('A/85//CA','A/85//CB'),'Glu91_CB_CG':('A/91//CB','A/91//CG'),'scaffold_N83_CA83':('A/83//N','A/83//CA')}

def design(inputs,plan,output):
 c=read_json(inputs);ex=read_json(verify(c['export']));bonds=read_json(verify(ex['bonds']));configs={};nearest=None;fixed_metal_direction=None
 for state in ('A','B'):
  mp=prepare_mapping(inputs,state,'La');atoms=mp['physical_atoms'];ids=[a['id'] for a in atoms];ix={sid:i for i,sid in enumerate(ids)};x=np.array([a['xyz_A'] for a in atoms]);mi=ix[target_source_id(next(a for a in mp['core'] if a['kind']=='metal'))];adj={sid:set() for sid in ids}
  for a,b in bonds:adj[a].add(b);adj[b].add(a)
  if nearest is None:
   qids={a['source_id'] for a in mp['core'] if a['kind']=='source'};distance,nearest=min((float(np.linalg.norm(x[ix[sid]]-x[mi])),sid) for sid in qids if atoms[ix[sid]]['element']=='O');fixed_metal_direction=(x[ix[nearest]]-x[mi])/distance
  modes={};t=np.zeros_like(x);t[mi]=fixed_metal_direction;modes['metal_toward_nearest_O']=dict(unit_tangent=t.tolist(),raw_tangent=t.tolist(),tangent_norm=1.,moving_source_ids=[ids[mi]],nearest_oxygen_source_id=nearest,metal_oxygen_distance_A=float(np.linalg.norm(x[ix[nearest]]-x[mi])),amplitude_units='angstrom')
  for name,(a,b) in AXES.items():
   graph={k:set(v) for k,v in adj.items()}
   if b not in graph[a]:raise InvalidArtifact('Requested source graph bond absent')
   graph[a].remove(b);graph[b].remove(a);move={b};todo=[b]
   while todo:
    now=todo.pop()
    for nxt in graph[now]-move:move.add(nxt);todo.append(nxt)
   if a in move:raise InvalidArtifact('Rotation cut is not graph bridge')
   moved=sorted(ix[sid] for sid in move);axis=x[ix[b]]-x[ix[a]];axis/=np.linalg.norm(axis);t=np.zeros_like(x);t[moved]=np.cross(axis,x[moved]-x[ix[b]]);norm=float(np.linalg.norm(t))
   if norm<=1e-12:raise InvalidArtifact('Zero physical displacement mode')
   modes[name]=dict(axis_source_ids=[a,b],unit_tangent=(t/norm).tolist(),raw_tangent=t.tolist(),tangent_norm=norm,moving_source_ids=sorted(move),moving_nonzero_count=int(np.count_nonzero(np.linalg.norm(t,axis=1)>1e-12)),RMS_A_per_radian=float(np.sqrt(np.mean(np.sum(t[moved]**2,axis=1)))),amplitude_units='radian')
  names=list(modes);us=np.array([np.array(modes[n]['unit_tangent']).ravel() for n in names]);gram=us@us.T
  ch=[i for i,a in enumerate(atoms) if a['element']=='H' and len(adj[a['id']])==1 and atoms[ix[next(iter(adj[a['id']]))]]['element']=='C'];h=[i for i,a in enumerate(atoms) if a['element']=='H' and i not in set(ch)];heavy=[i for i,a in enumerate(atoms) if a['element']!='H']
  configs[state]=dict(physical_atoms=atoms,modes=modes,mode_order=names,normalized_mode_gram=gram.tolist(),subsets={'C_bound_H':ch,'other_H':h,'heavy':heavy})
 write_new(output,dict(protocol='original_H_frozen_f_physical_projections_v1',inputs=record(inputs),plan=record(plan),implementation=record(__file__),configurations=configs))

def analyze(design_path,assembled,output):
 d=read_json(design_path);verify(d['inputs']);verify(d['plan']);verify(d['implementation']);rows={}
 for state,conf in d['configurations'].items():
  data={m:read_json(Path(assembled)/(m+'_'+state+'.json')) for m in ('La','Dy')};ids=[a['id'] for a in conf['physical_atoms']];x=np.array([a['xyz_A'] for a in conf['physical_atoms']]);gs={}
  for m,r in data.items():
   if [a['id'] for a in r['physical_atoms']]!=ids or not np.array_equal([a['xyz_A'] for a in r['physical_atoms']],x):raise InvalidArtifact('Assembly physical coordinates/order differ')
   if r['total_gradient_kcal_mol_A'] is None:raise InvalidArtifact('Complete finite gradient unavailable')
   gs[m]={'electronic':np.array(r['electronic']['gradient_kcal_mol_A']),'classical':np.array(r['classical_gradient_kcal_mol_A']),'total':np.array(r['total_gradient_kcal_mol_A'])}
  modes={}
  for name,mode in conf['modes'].items():
   u=np.array(mode['unit_tangent']);t=np.array(mode['raw_tangent']);comp={}
   for component in ('electronic','classical','total'):
    vals={m:float(np.sum(gs[m][component]*u)) for m in ('La','Dy')};delta=vals['Dy']-vals['La'];comp[component]=dict(gradient_load_kcal_mol_A=vals,common_mean_gradient_load_kcal_mol_A=(vals['Dy']+vals['La'])/2,Dy_minus_La_gradient_load_kcal_mol_A=delta,Dy_minus_La_force_load_kcal_mol_A=-delta,raw_generalized_gradient={m:float(np.sum(gs[m][component]*t)) for m in ('La','Dy')},raw_generalized_units='kcal/mol/'+mode['amplitude_units'])
   modes[name]=dict(components=comp,tangent_norm=mode['tangent_norm'],moving_count=len(mode['moving_source_ids']),moving_nonzero_count=mode.get('moving_nonzero_count',1),RMS_A_per_radian=mode.get('RMS_A_per_radian'),axis_source_ids=mode.get('axis_source_ids'),nearest_oxygen_source_id=mode.get('nearest_oxygen_source_id'),metal_oxygen_distance_A=mode.get('metal_oxygen_distance_A'))
  subsets={}
  for name,inds in conf['subsets'].items():
   sub={}
   for component in ('electronic','classical','total'):
    l=gs['La'][component][inds];y=gs['Dy'][component][inds];mean=(l+y)/2;delta=y-l
    sub[component]=dict(La_norm=float(np.linalg.norm(l)),Dy_norm=float(np.linalg.norm(y)),common_mean_norm=float(np.linalg.norm(mean)),differential_norm=float(np.linalg.norm(delta)),common_RMS_cartesian=float(np.sqrt(np.mean(mean**2))),differential_RMS_cartesian=float(np.sqrt(np.mean(delta**2))))
   subsets[name]=dict(atom_count=len(inds),units='kcal/mol/angstrom',components=sub)
  rows[state]=dict(modes=modes,subsets=subsets,normalized_mode_gram=conf['normalized_mode_gram'],mode_order=conf['mode_order'],assemblies={m:record(Path(assembled)/(m+'_'+state+'.json')) for m in ('La','Dy')})
 write_new(output,dict(design=record(design_path),rows=rows,full_hybrid_derivative_qualified=False,solvent_consistent=False,optimization_performed=False,affinity=None))

if __name__=='__main__':
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='operation',required=True);a=s.add_parser('design');a.add_argument('--inputs',required=True);a.add_argument('--plan',required=True);a.add_argument('--output',required=True);a=s.add_parser('analyze');a.add_argument('--design',required=True);a.add_argument('--assembled',required=True);a.add_argument('--output',required=True);a=p.parse_args()
 if a.operation=='design':design(a.inputs,a.plan,a.output)
 else:analyze(a.design,a.assembled,a.output)
