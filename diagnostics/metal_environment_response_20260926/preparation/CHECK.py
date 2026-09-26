"""Real-fixture input/mapping checks; no scientific executable or energies."""
from pathlib import Path
import argparse, json, hashlib, sys
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--repository',type=Path,required=True);a=p.parse_args()
sys.path.insert(0,str(a.repository/'scripts'));from affordable_response import cap_jacobians
x=json.loads(a.inputs.read_text());results=[]
def verify(r):
 z=Path(r['path']);assert hashlib.sha256(z.read_bytes()).hexdigest()==r['sha256'];return z
for z in x['endpoints'].values():verify(z['xyz'])
env={z:json.loads(verify(r['atoms']).read_text()) for z,r in x['environments'].items()};maps=json.loads(verify(x['core_mapping']).read_text());changes=[i for i,(u,v) in enumerate(zip(env['A'],env['B'])) if u!=v];assert changes==[x['perturbation']['pointcharge_index']];assert [r['charge_e'] for r in env['A']]==[r['charge_e'] for r in env['B']];results.append({'check':'same_inventory_charge_only_declared_H_moves','status':'pass','changed_indices':changes})
qa=np.loadtxt(verify(x['endpoints']['Ca']['xyz']),skiprows=2,usecols=(1,2,3));ql=np.loadtxt(verify(x['endpoints']['La']['xyz']),skiprows=2,usecols=(1,2,3));assert np.array_equal(qa,ql);results.append({'check':'paired_coordinate_identity','status':'pass','max_delta_A':float(abs(qa-ql).max())})
assert not any(r.get('mapping_status')=='explicitly_unresolved' for r in maps);assert sum(r.get('mapping_status')=='carver_local_hydrogen' for r in maps)==3;results.append({'check':'complete_source_and_synthetic_atom_accounting','status':'pass'})
caperrs=[]
for m in maps:
 if m['kind']!='cap':continue
 ra=np.array(m['retained_xyz_A']);rb=np.array(m['omitted_xyz_A']);ja,jb=cap_jacobians(ra,rb,m['length_A']);h=1e-5
 def f(v,w):return v+m['length_A']*(w-v)/np.linalg.norm(w-v)
 for which,jac in [(0,ja),(1,jb)]:
  numeric=np.zeros((3,3))
  for i in range(3):
   d=np.eye(3)[i]*h;numeric[:,i]=(f(ra+d,rb)-f(ra-d,rb))/(2*h) if which==0 else (f(ra,rb+d)-f(ra,rb-d))/(2*h)
  caperrs.append(float(abs(numeric-jac).max()))
assert max(caperrs)<1e-8;results.append({'check':'physical_cap_Jacobians_selected_real_coordinate_FD','status':'pass','max_absolute_error':max(caperrs),'step_A':1e-5,'tolerance':1e-8})
heavy=np.array([r['xyz_A'] for r in maps if r['element']!='H']);before=np.array(x['perturbation']['hydrogen_A_A']);after=np.array(x['perturbation']['hydrogen_B_A']);near=[float(np.linalg.norm(heavy-h,axis=1).min()) for h in (before,after)];assert min(near)>1.2;results.append({'check':'modeled_H_to_QM_heavy_no_short_clash','status':'pass','before_after_A':near,'minimum_A':1.2})
result={'status':'pass','checks':results,'molecular_evaluations':0,'scientific_integration_tests':0,'notes':'Coordinate finite differences validate analytic cap Jacobians, not an electronic energy or force.'};assert not a.output.exists();a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
