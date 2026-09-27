"""Real source-map checks and six-cell dry-run only."""
import argparse,sys,tempfile,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();sys.path.insert(0,str(a.repository/'scripts'))
from affordable_common import read_json,verify,write_new
from metal_environment_reference import prepare,validate
from metal_environment_force_checks import check_pins,xyz_data
from affordable_response import cap_jacobians
x=read_json(a.inputs);check_pins(x);checks=[]
_,ca,cq=xyz_data(verify(x['endpoints']['Ca']['xyz']));_,la,lq=xyz_data(verify(x['endpoints']['La']['xyz']));assert np.array_equal(cq,lq) and ca[1:]==la[1:] and len(cq)==52
z={'H':1,'C':6,'N':7,'O':8,'S':16,'Ca':20,'La':57};n={m:sum(z[e] for e in elements)-x['endpoints'][m]['charge'] for m,elements in [('Ca',ca),('La',la)]};assert all(e%2==0 for e in n.values());checks.append({'check':'paired52geometry_state_parity','status':'pass','all_electron_counts':n,'La_explicit_electrons':n['La']-46})
e={k:read_json(verify(v['atoms'])) for k,v in x['environments'].items()};assert len(e['A'])==len(e['B'])==1880;assert [i for i,(u,v) in enumerate(zip(e['A'],e['B'])) if u!=v]==[1292];assert [v['charge_e'] for v in e['A']]==[v['charge_e'] for v in e['B']];assert abs(sum(v['charge_e'] for v in e['A'])+4)<1e-10;checks.append({'check':'only_declared_environment_H_moves_charge_closure','status':'pass'})
mapping=read_json(verify(x['core_mapping']));err=[]
for cap in mapping:
 if cap['kind']!='cap':continue
 r=np.array(cap['retained_xyz_A']);s=np.array(cap['omitted_xyz_A']);l=cap['length_A'];ja,jb=cap_jacobians(r,s,l)
 def f(u,v):return u+l*(v-u)/np.linalg.norm(v-u)
 for k in range(3):
  d=np.eye(3)[k]*1e-5;err.extend([float(abs((f(r+d,s)-f(r-d,s))/2e-5-ja[:,k]).max()),float(abs((f(r,s+d)-f(r,s-d))/2e-5-jb[:,k]).max())])
assert max(err)<1e-8;checks.append({'check':'all11cap_Jacobians','status':'pass','maximum_absolute_error':max(err),'tolerance':1e-8,'FD_step_A':1e-5})
with tempfile.TemporaryDirectory() as d:
 r=prepare(a.inputs,Path(__file__).with_name('PLAN.md'),Path(d)/'ready',workers=6,mpi_ranks=57);validate(verify(r['manifest']));checks.append({'check':'existing6cellprepare_dry_run','status':'pass'})
write_new(a.output,{'status':'pass','checks':checks,'molecular_evaluations':0});print(json.dumps(checks,indent=2))
