"""Real prepared Hans source, state, geometry, boundary and link-map checks."""
import argparse,json,sys
from pathlib import Path
import numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();sys.path.insert(0,str(a.repository.resolve()/'scripts'))
 from affordable_common import read_json,verify,record,write_new
 c=read_json(a.inputs);ex=read_json(verify(c['export']));atoms=read_json(verify(ex['atoms']));by={v['id']:v for v in atoms};parent=read_json(verify(c['physical_source']))['atoms'];assert len(atoms)==1887
 assert all(v['id']==w['id'] and v['xyz_A']==w['xyz_A'] for v,w in zip(atoms,parent))
 maps={k:read_json(verify(v)) for k,v in c['core_mapping'].items()};field=read_json(verify(c['environment_atoms']));boundary=read_json(verify(c['boundary_mapping']));bonds=read_json(verify(ex['bonds']));adj={sid:set() for sid in by}
 for u,v in bonds:adj[u].add(v);adj[v].add(u)
 for cfg in c['configurations'].values():
  verify(cfg['pointcharges'])
  for metal,e in cfg['endpoints'].items():
   lines=verify(e['xyz']).read_text().splitlines();assert int(lines[0])==195 and len(lines)==197; xyz=np.array([[float(v) for v in l.split()[1:]] for l in lines[2:]]);assert e['charge']==-1;assert e['multiplicity']==({'La':1,'Dy':6}[metal]);assert (e['all_electron_count']-(e['multiplicity']-1))%2==0
   assert np.array_equal(xyz,np.array([v['xyz_A'] for v in read_json(verify(cfg['core_mapping']))]))
 assert c['configurations']['A']['pointcharges']==c['configurations']['B']['pointcharges'];assert len(field)==1696 and abs(sum(v['charge_e'] for v in field)-7)<1e-10
 selected=set(boundary['all_selected_QM_source_ids']);omitted=set(boundary['removed_MM1_CA_ids']);assert not selected.intersection(v['id'] for v in field);assert not omitted.intersection(v['id'] for v in field)
 assert set(by)==selected|omitted|{v['id'] for v in field if v['kind']!='spectator_formal_monopole'}
 assert sorted((v['site'],v['element'],v['charge_e']) for v in field if v['kind']=='spectator_formal_monopole')==[('EF1','La',3.),('EF2','La',3.),('EF4','Na',1.)]
 physA={sid:np.array(v['xyz_A']) for sid,v in by.items()};physB={sid:x.copy() for sid,x in physA.items()}
 for v in maps['B']:
  if v['kind']=='source':physB[v['id']]=np.array(v['xyz_A'])
 changed={sid for sid in by if not np.array_equal(physA[sid],physB[sid])};assert changed==set(c['perturbation']['changed_source_ids'])
 bond_residual=max(abs(np.linalg.norm(physA[u]-physA[v])-np.linalg.norm(physB[u]-physB[v])) for u,v in bonds);assert bond_residual<1e-12
 candidates=[]
 for sid in changed:
  excluded={sid}|adj[sid]|set().union(*(adj[n] for n in adj[sid]))
  for other in set(by)-excluded:
   candidates.append((float(np.linalg.norm(physB[sid]-physB[other])),sid,other,float(np.linalg.norm(physA[sid]-physA[other]))))
 closest=min(candidates)
 jac_error=0.;h=1e-5
 for v in maps['A']:
  if v['kind']!='cap':continue
  x=np.array(v['retained_xyz_A']);y=np.array(v['omitted_xyz_A']);length=v['length_A'];f=lambda x,y:x+length*(y-x)/np.linalg.norm(y-x)
  for which in ('retained','omitted'):
   J=np.array(v['jacobian_'+which]);numeric=np.zeros((3,3))
   for axis in range(3):
    d=np.eye(3)[axis]*h;numeric[:,axis]=(f(x+d,y)-f(x-d,y))/(2*h) if which=='retained' else (f(x,y+d)-f(x,y-d))/(2*h)
   jac_error=max(jac_error,float(np.max(abs(J-numeric))))
 assert jac_error<1e-8
 result={'status':'passed','inputs':record(a.inputs),'checks':['exact native source IDs and normalized coordinates','all four paired state/XYZ mappings and electron parity','fixed A/B field and QM/MM charge closure','all physical atoms retained in QM/MM/zero-charge boundary mapping','actual La/La/Na spectator occupancy','Asp85 exact five-atom rotation and covalent-length preservation','four cap analytic Jacobians versus coordinate finite differences'], 'max_covalent_length_change_A':bond_residual,'closest_moved_nonbonded_B_distance_A':closest[0],'closest_pair':closest[1:3],'same_pair_A_distance_A':closest[3],'max_displacement_A':max(float(np.linalg.norm(physB[sid]-physA[sid])) for sid in by),'max_cap_jacobian_residual':jac_error,'molecular_energy_force_calls':0,'physical_surface_qualification':'not tested','implementation':record(__file__)};write_new(a.output,result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
