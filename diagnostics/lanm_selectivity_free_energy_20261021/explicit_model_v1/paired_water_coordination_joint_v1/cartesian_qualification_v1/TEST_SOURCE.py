#!/usr/bin/env python3
"""Existing saved W8 outputs plus synthetic W9 parser controls; no new chemistry."""
import copy,importlib.util,json,os,sys,time
from pathlib import Path
H=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('joint_cartesian_controls',H/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def main(out):
 t=time.monotonic();cpu=time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False);checks=[]
 m.require(all(os.environ.get(k)=='1' for k in m.THREADS),'One thread required');cfg=m.default_config();m.metadata(cfg)
 def ok(v,label):
  if not v:raise AssertionError(label)
  checks.append(label)
 def bad(fn,label):
  try:fn()
  except (ValueError,KeyError):checks.append(label)
  else:raise AssertionError(label)
 ok(cfg['candidate_result'] is None and cfg['candidate_review'] is None,'Candidate/review remain unbound')
 ok(m.ARGS==['system.xyz','Y','Y','0.00001'] and sum(6*n for n in m.ATOMS.values())==636,'Exact arguments and derived636 ceiling')
 for key,val in (('maximum_native_programs',5),('maximum_internal_FD_energies',637),('maximum_workers',5)):
  c=copy.deepcopy(cfg);c[key]=val;bad(lambda:m.metadata(c),'Changed ceiling rejected '+key)
 c=copy.deepcopy(cfg);c['executables']['testgrad']['sha256']='0'*64;bad(lambda:m.metadata(c),'Changed testgrad binary rejected')
 bad(lambda:m.candidate_gate(cfg),'No candidate gate opens before result binding')
 native=m.module(cfg['native_helper'],'cartesian_saved_control_native');helper=m.module(cfg['numerical_helper'],'cartesian_saved_control_parser');old=m.read(cfg['previous_transfer'])
 for gid in m.CELLS[:2]:
  row=old['rows'][gid];check=helper.gradient(m.verify(row['cartesian_check']['output']),25,True,native);ok(check==row['cartesian_check'],'Exact admitted17612 Cartesian replay '+gid);ok(m.compare(check,row['model'])['passed'],'Existing selected E/F comparison '+gid)
 def output(name,numer_error=0.,n=28):
  text='Tinker Version 26.2\n Total Potential Energy : 0.0 Kcal/mole\n Cartesian Gradient Breakdown over Individual Atoms :\n'
  for i in range(1,n+1):text+=f'Anlyt {i} 0.0 0.0 0.0 0.0\nNumer {i} {numer_error if i==1 else 0.0:.9g} 0.0 0.0 0.0\n'
  p=out/name;p.write_text(text);return p
 p=output('synthetic28.out');check=helper.gradient(p,28,True,native);value={'energy':0.,'forces':[[0.,0.,0.] for _ in range(28)]}
 ok(m.compare(check,value)['passed'],'Synthetic28 complete numerical arrays pass');ok(check['fd_total_separation_A']==1e-5 and check['fd_each_side_displacement_A']==5e-6,'Unchanged derivative separation interpretation')
 p=output('synthetic28_boundary.out',.001);check=helper.gradient(p,28,True,native);ok(m.compare(check,value)['passed'],'Cartesian1e-3 inclusive')
 p=output('synthetic28_failed.out',.001001);check=helper.gradient(p,28,True,native);bad(lambda:m.compare(check,value),'Cartesian beyond1e-3 fails')
 ok(not check['fd_passed'] and check['max_component_error_kcal_mol_A']>.001 and len(check['numerical_gradient'])==28,'Failed tolerance retains full parsed discrepancy')
 p=output('synthetic28_energy.out');check=helper.gradient(p,28,True,native);v=copy.deepcopy(value);v['energy']=.0001;ok(m.compare(check,v)['passed'],'Saved energy1e-4 inclusive');v['energy']=.000101;bad(lambda:m.compare(check,v),'Saved energy beyond1e-4 fails')
 v=copy.deepcopy(value);v['forces'][0][0]=.001;ok(m.compare(check,v)['passed'],'Saved analytic force1e-3 inclusive');v['forces'][0][0]=.001001;bad(lambda:m.compare(check,v),'Saved analytic force beyond1e-3 fails')
 p=output('synthetic27_missing.out',n=27);bad(lambda:helper.gradient(p,28,True,native),'Missing atom/component rejected')
 p=output('synthetic28_nonfinite.out');p.write_text(p.read_text().replace('Numer 1 0','Numer 1 NaN',1));bad(lambda:helper.gradient(p,28,True,native),'Nonfinite native marker rejected')
 p=output('synthetic28_order.out');p.write_text(p.read_text().replace('Numer 2 ','Numer 3 ',1));bad(lambda:helper.gradient(p,28,True,native),'Numerical atom order rejected')
 # Record-only full-program completion precedes tolerance admission, preserving
 # the parsed check in failure rows. This is metadata inspection, not execution.
 src=(H/'RUN.py').read_text();a=src.index("with lock:target['check']=check;result['internal_FD_energies_completed']+=6*ATOMS[gid]");b=src.index("admission=compare(check,row['selected_model'])",a)
 ok(a<b and src.count("result['internal_FD_energies_completed']+=6*ATOMS[gid]")==1,'Complete numerical-array count recorded once before tolerance admission')
 result={'status':'passed','implementation':m.pin(__file__),'source':m.pin(H/'RUN.py'),'declaration':m.pin(H/'DECLARATION.md'),'saved17612':cfg['previous_transfer'],'checks':checks,'checks_passed':len(checks),'native_program_invocations':0,'new_internal_FD_energies':0,'actual17676_outcomes_opened':False,'candidate_bound':False,'synthetic_outputs_explicit':True,'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1}
 m.save(out/'RESULT.json',result);print(json.dumps({'status':'passed','checks':len(checks),'result':m.pin(out/'RESULT.json')}))
if __name__=='__main__':main(sys.argv[1])
