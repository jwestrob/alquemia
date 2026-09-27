"""Exact-source Thr159 partition diagnostic preparation; no energy Context."""
from pathlib import Path
import argparse,copy,hashlib,json,sys
import numpy as np

def rec(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):return json.loads(Path(p).read_text())
def ver(r):assert rec(r['path'])==r,r;return Path(r['path'])
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();sys.path.insert(0,str(a.repository.resolve()/'scripts'))
 from affordable_state import source_protein
 from affordable_response import cap_jacobians
 import openmm
 base=read(a.inputs);env={s:read(ver(v['atoms'])) for s,v in base['environments'].items()};oldmap=read(ver(base['core_mapping']));atoms,pos,q,physical,excluded=source_protein(ver(base['source']));byid={v['id']:i for i,v in enumerate(physical)};prefix='A/159/ /';ca,cb=[byid[prefix+n] for n in ('CA','CB')]
 # Use the actual topology, cut only the validated source CA-CB bond.
 from openmm import app
 pdb=app.PDBFile(str(ver(base['source'])));mod=app.Modeller(pdb.topology,pdb.positions);mod.delete([r for r in mod.topology.residues() if r.chain.id!='A' or r.name in ('PQQ','LA','CA','CE')]);top=list(mod.topology.atoms());assert len(top)==len(atoms)
 edges={tuple(sorted((u.index,v.index))) for u,v in mod.topology.bonds()};assert tuple(sorted((ca,cb))) in edges
 adj={i:set() for i in range(len(atoms))}
 for u,v in edges:
  if {u,v}!={ca,cb}:adj[u].add(v);adj[v].add(u)
 seen={cb};todo=[cb]
 while todo:
  i=todo.pop()
  for j in adj[i]-seen:seen.add(j);todo.append(j)
 assert {atoms[i].name for i in seen}=={'CB','HB','OG1','HG1','CG2','HG21','HG22','HG23'}
 assert all(atoms[i].residue.id=='159' and atoms[i].residue.chain.id=='A' for i in seen)
 srcids=[physical[i]['id'] for i in sorted(seen)];assert not any(m.get('source_id') in srcids for m in oldmap)
 allres=[i for i,v in enumerate(physical) if v['id'].startswith(prefix)];assert abs(sum(q[i] for i in allres))<1e-12
 oldA={v['id']:v for v in env['A']};assert all(abs(oldA[physical[i]['id']]['charge_e']-q[i])<1e-12 for i in allres)
 removed=seen|{ca};retained=[i for i in allres if i not in removed];recipients=[byid[prefix+n] for n in ('N','C')];increment=-sum(q[i] for i in retained)/2
 ledger={'fragment_id':'A:THR159','formal_charge':0,'original_residue_charge_e':float(sum(q[i] for i in allres)),'source_qm_indices':sorted(seen),'source_qm_ids':srcids,'removed_source_indices':sorted(removed),'removed_source_ids':[physical[i]['id'] for i in sorted(removed)],'recipients':recipients,'recipient_ids':[physical[i]['id'] for i in recipients],'each_increment_e':float(increment),'retained_MM_charge_before_e':float(sum(q[i] for i in retained))}
 pa=np.array(oldA[prefix+'CB']['xyz_A']);pb=np.array(oldA[prefix+'CA']['xyz_A']);cap=pa+1.09*(pb-pa)/np.linalg.norm(pb-pa);ja,jb=cap_jacobians(pa,pb,1.09)
 result={'protocol_id':'nikasha_Thr159_expanded_region_partition_preparation_v1','base_inputs':rec(a.inputs),'plan':rec(Path(__file__).with_name('PLAN.md')),'source':base['source'],'source_state':base['source_state'],'assembly':base['assembly'],'microstate':base['microstate'],'explicit_waters':base['explicit_waters'],'implementation':rec(__file__),'source_graph_implementation':rec(a.repository/'scripts/affordable_state.py'),'openmm_version':openmm.__version__,'forcefield':read(ver(base['boundary_mapping']))['forcefield'],'orca':base['orca'],'endpoints':{},'environments':{},'new_boundary_ledger':ledger,'physical_configurations':'Exactly original small-region A/B physical atoms and coordinates; extra synthetic cap at Thr CB-CA; B HG1 now QM.','molecular_evaluations':0,'complete_hybrid_status':'unsupported_missing_cross_repulsion_dispersion_and_boundary_mechanics','classification':None}
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);checks=[];znumber={'H':1,'C':6,'N':7,'O':8,'S':16,'Ca':20,'La':57};expandedmaps={}
 for label in ('A','B'):
  old={v['id']:v for v in env[label]};new=[]
  for v in env[label]:
   if v['id'] in ledger['removed_source_ids']:continue
   vv=copy.deepcopy(v)
   if vv['id'] in ledger['recipient_ids']:vv['charge_e']+=float(increment)
   new.append(vv)
  assert len(new)==9078;assert abs(sum(v['charge_e'] for v in new)-base['environments'][label]['charge_e'])<1e-10
  write(out/f'environment_{label}.json',new);pc=out/f'environment_{label}.pc';pc.write_text(str(len(new))+'\n'+''.join(' '.join(format(t,'.17g') for t in (v['charge_e'],*v['xyz_A']))+'\n' for v in new));result['environments'][label]={'pointcharges':rec(pc),'atoms':rec(out/f'environment_{label}.json'),'count':len(new),'charge_e':float(sum(v['charge_e'] for v in new))}
  mapping=copy.deepcopy(oldmap)
  for sid in srcids:
   v=old[sid];mapping.append({'qm_index':len(mapping),'element':v['element'],'kind':'protein_source','id':sid,'source_id':sid,'source_index':v['source_index'],'xyz_A':v['xyz_A']})
  mapping.append({'qm_index':len(mapping),'element':'H','kind':'cap','id':'cap/A:THR159','xyz_A':cap.tolist(),'retained_source_id':prefix+'CB','omitted_source_id':prefix+'CA','retained_xyz_A':pa.tolist(),'omitted_xyz_A':pb.tolist(),'length_A':1.09,'jacobian_retained':ja.tolist(),'jacobian_omitted':jb.tolist(),'mapping':'r_CB + 1.09 (r_CA-r_CB)/norm(r_CA-r_CB)'})
  assert len(mapping)==63;expandedmaps[label]=mapping;write(out/f'core_mapping_{label}.json',mapping)
  for metal,oldend in base['endpoints'].items():
   lines=ver(oldend['xyz']).read_text().splitlines();oldrows=[l.split() for l in lines[2:]];assert len(oldrows)==54
   rows=[(r[0],[float(v) for v in r[1:4]]) for r in oldrows]+[(r['element'],r['xyz_A']) for r in mapping[54:]]
   xyz=out/f'{metal}_{label}.xyz';xyz.write_text('63\n'+f'{metal} {label} Thr159 expanded partition, charge={oldend["charge"]}, mult=1\n'+''.join(el+' '+' '.join(format(v,'.17g') for v in xyz)+'\n' for el,xyz in rows))
   n_all=sum(znumber[el] for el,xyz in rows)-oldend['charge'];ecp=46 if metal=='La' else 0;assert (n_all-ecp)%2==0
   result['endpoints'][f'{metal}_{label}']={'metal':metal,'environment':label,'xyz':rec(xyz),'charge':oldend['charge'],'multiplicity':1,'physical_electrons':n_all,'native_r2scan3c_ecp_core_electrons':ecp,'native_r2scan3c_explicit_electrons':n_all-ecp,'core_mapping':rec(out/f'core_mapping_{label}.json'),'pointcharges':result['environments'][label]['pointcharges']}
  # All old MM atoms are retained in field, promoted to QM, or recorded omitted CA.
  assert {v['id'] for v in env[label]}=={v['id'] for v in new}|set(srcids)|{prefix+'CA'}
  for v in new:assert v['xyz_A']==old[v['id']]['xyz_A']
  for r in mapping[54:-1]:assert r['xyz_A']==old[r['id']]['xyz_A']
  assert mapping[:54]==oldmap
  checks.append({'environment':label,'old_core_preserved':True,'all_old_environment_atoms_accounted_for':True,'new_QM_atoms':8,'new_caps':1,'removed_MM_atoms':9,'new_environment_charge_e':result['environments'][label]['charge_e']})
 assert result['environments']['A']['pointcharges']['sha256']==result['environments']['B']['pointcharges']['sha256']
 changed=[i for i,(u,v) in enumerate(zip(expandedmaps['A'],expandedmaps['B'])) if u!=v];assert changed==[54+srcids.index(prefix+'HG1')]
 result['checks']=checks;result['changed_QM_indices_A_to_B']=changed;result['checks_passed']=True;write(out/'INPUTS.json',result);print(json.dumps({'inputs':str(out/'INPUTS.json'),'ledger':ledger,'endpoints':{k:{f:v[f] for f in ('charge','physical_electrons','native_r2scan3c_explicit_electrons')} for k,v in result['endpoints'].items()},'changed_indices':changed},indent=2))
if __name__=='__main__':main()
