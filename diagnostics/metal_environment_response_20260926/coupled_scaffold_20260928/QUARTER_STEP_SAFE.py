"""One declared quarter-step, same actual model; preserve failed coarse check."""
from pathlib import Path
import sys,os,time,resource,json
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record
from result_io import write_new,write_configuration
from metal_environment_component_checks import build_system,GROUPS
if len(sys.argv)==3 and sys.argv[1]=='--render':
 print(Path(sys.argv[2]).read_text(),end='');raise SystemExit(0)
if len(sys.argv)!=2:raise SystemExit('usage: QUARTER_STEP_SAFE.py NEW_OUTPUT_DIRECTORY or --render SAVED_RESULT_JSON')
if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('allocation required')
start=time.monotonic();b=ROOT/'workspaces/metal_environment_response_20260926/coupled_scaffold_classical_v1';out=Path(sys.argv[1]).resolve();out.mkdir(exist_ok=False);ledger=read_json(b/'LEDGER.json');p=read_json(b/'particles.json');art={k:read_json(b/k) for k in ['particles.json','exceptions.json','bonded_terms.json']};old=read_json(b/'RESULT.json');x=np.array(read_json(b/'origin_physical.json')['xyz_A']);oldplus=np.array(read_json(b/'plus_physical.json')['xyz_A']);moving=np.where(np.linalg.norm(oldplus-x,axis=1)>1e-12)[0];ids={v['id']:i for i,v in enumerate(p)};axis=x[ids['A/83//CA']]-x[ids['A/83//N']];axis/=np.linalg.norm(axis);pivot=x[ids['A/83//CA']];h=np.deg2rad(.0125);pool={}
for label,angle in [('minus',-h),('plus',h)]:
 xx=x.copy();v=x[moving]-pivot;xx[moving]=pivot+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+np.outer(v@axis,axis)*(1-np.cos(angle));pool[label]=xx;write_new(out/f'{label}_physical.json',{'source_ids':[v['id'] for v in p],'xyz_A':xx.tolist()})
results={};calls=0
for metal in ['La','Dy']:
 system,_=build_system(ledger,art,metal);integ=mm.VerletIntegrator(.001);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'));ee={}
 for label,xx in pool.items():
  ctx.setPositions(xx*.1);ee[label]={};ff={}
  for group,name in enumerate(GROUPS):
   st=ctx.getState(getEnergy=True,getForces=True,groups={group});calls+=1;ee[label][name]=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole);ff[name]=st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom)
  write_configuration(out/f'{metal}_{label}_receipt.json',energies=ee[label],forces=ff,energy_units='kcal/mol',force_units='kcal/mol/angstrom',quantity='force_not_gradient',source_ids=[v['id'] for v in p],geometry=record(out/f'{label}_physical.json'),ledger=record(b/'LEDGER.json'),metal=metal,configuration=label)
 checks=[]
 for row in old['results'][metal]['checks']:
  name=row['component'];fd=(ee['plus'][name]-ee['minus'][name])/(2*h);err=fd-row['analytic'];half=next(v for v in read_json(b/'halfstep_receipt_recovery_v2/RESULT.json')['results'][metal]['checks'] if v['component']==name);rich=(4*fd-half['halfstep_derivative'])/3;checks.append({'component':name,'analytic':row['analytic'],'quarterstep_derivative':fd,'quarterstep_residual':err,'same_tolerance':row['tolerance'],'quarterstep_passed':abs(err)<=row['tolerance'],'coarse_residual':row['residual'],'halfstep_residual':half['halfstep_residual'],'half_over_quarter_residual_ratio':half['halfstep_residual']/err if abs(err)>1e-14 else None,'coarse_over_quarter_residual_ratio':row['residual']/err if abs(err)>1e-14 else None,'Richardson_derivative':rich,'Richardson_discrepancy':rich-row['analytic'],'Richardson_admission':'not_defined_no_new_gate'})
 results[metal]={'energies_kcal_mol':ee,'checks':checks};del ctx;del integ
r={'results':results,'status':'passed' if all(v['quarterstep_passed'] for row in results.values() for v in row['checks']) else 'failed','coarse_result':record(b/'RESULT.json'),'implementation':record(__file__),'plan':record(Path(__file__).with_name('QUARTER_STEP_PLAN.md')),'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':time.monotonic()-start,'peakRSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'energy_force_queries':calls,'new_configurations':2,'QM_calls':0,'original_failures_preserved':True};write_new(out/'RESULT.json',r);print((out/'RESULT.json').read_text(),end='')
