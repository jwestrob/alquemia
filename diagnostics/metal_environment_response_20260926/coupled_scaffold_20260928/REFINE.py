"""One declared halfstep, same actual model; preserve failed coarse check."""
from pathlib import Path
import sys,os,time,resource,json
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new
from metal_environment_component_checks import build_system,GROUPS
if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('allocation required')
start=time.monotonic();b=ROOT/'workspaces/metal_environment_response_20260926/coupled_scaffold_classical_v1';out=b/'halfstep_v1';out.mkdir(exist_ok=False);ledger=read_json(b/'LEDGER.json');p=read_json(b/'particles.json');art={k:read_json(b/k) for k in ['particles.json','exceptions.json','bonded_terms.json']};old=read_json(b/'RESULT.json');x=np.array(read_json(b/'origin_physical.json')['xyz_A']);oldplus=np.array(read_json(b/'plus_physical.json')['xyz_A']);moving=np.where(np.linalg.norm(oldplus-x,axis=1)>1e-12)[0];ids={v['id']:i for i,v in enumerate(p)};axis=x[ids['A/83//CA']]-x[ids['A/83//N']];axis/=np.linalg.norm(axis);pivot=x[ids['A/83//CA']];h=np.deg2rad(.025);pool={}
for label,angle in [('minus',-h),('plus',h)]:
 xx=x.copy();v=x[moving]-pivot;xx[moving]=pivot+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+np.outer(v@axis,axis)*(1-np.cos(angle));pool[label]=xx;write_new(out/f'{label}_physical.json',{'source_ids':[v['id'] for v in p],'xyz_A':xx.tolist()})
results={};calls=0
for metal in ['La','Dy']:
 system,_=build_system(ledger,art,metal);integ=mm.VerletIntegrator(.001);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'));ee={}
 for label,xx in pool.items():
  ctx.setPositions(xx*.1);ee[label]={}
  for group,name in enumerate(GROUPS):st=ctx.getState(getEnergy=True,getForces=True,groups={group});calls+=1;ee[label][name]=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
 checks=[]
 for row in old['results'][metal]['checks']:
  name=row['component'];fd=(ee['plus'][name]-ee['minus'][name])/(2*h);err=fd-row['analytic'];rich=(4*fd-row['central_difference'])/3;checks.append({'component':name,'analytic':row['analytic'],'halfstep_derivative':fd,'halfstep_residual':err,'same_tolerance':row['tolerance'],'halfstep_passed':abs(err)<=row['tolerance'],'coarse_residual':row['residual'],'coarse_over_half_residual_ratio':row['residual']/err if abs(err)>1e-14 else None,'Richardson_derivative':rich,'Richardson_discrepancy':rich-row['analytic'],'Richardson_admission':'not_defined_no_new_gate'})
 results[metal]={'energies_kcal_mol':ee,'checks':checks};del ctx;del integ
r={'results':results,'status':'passed' if all(v['halfstep_passed'] for row in results.values() for v in row['checks']) else 'failed','coarse_result':record(b/'RESULT.json'),'implementation':record(__file__),'plan':record(Path(__file__).with_name('REFINEMENT_PLAN.md')),'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':time.monotonic()-start,'peakRSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'energy_force_queries':calls,'new_configurations':2,'QM_calls':0,'original_failures_preserved':True};write_new(out/'RESULT.json',r);print(json.dumps(r,indent=2))
