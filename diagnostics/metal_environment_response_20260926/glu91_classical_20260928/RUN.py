"""Evaluate only prescribed Glu91 configurations with exact saved classical ledger."""
from pathlib import Path
import sys,os,time,resource
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify
from metal_environment_component_checks import build_system,GROUPS
sys.path.insert(0,str(Path(__file__).parents[1]/'coupled_scaffold_20260928'))
from result_io import write_new
W=ROOT/'workspaces/metal_environment_response_20260926'

def prepare():
 c=read_json(W/'glu91_motion_v1/INPUTS.json');lp=W/'scaffold_origin_DQ2_v1/Hans8DQ2/LEDGER.json';ledger=read_json(lp);origin=read_json(verify(ledger['inputs']));assert c['origin_inputs']==ledger['inputs']
 for key in ['boundary_mapping','export','environment_atoms','spectator_occupancy','charge_accounting','source_state']:
  assert c[key]==origin[key],key
 for key in ['boundary_mapping','export','environment_atoms','source_state']:verify(c[key])
 verify(ledger['native_system']);verify(ledger['reused_builder'])
 artifacts={name:read_json(verify(pin)) for name,pin in ledger['artifacts'].items()};p=artifacts['particles.json'];ids={v['id']:i for i,v in enumerate(p)};x0=np.array([v['xyz_A'] for v in p]);env=read_json(verify(c['environment_atoms']));xs={};targetid=None;checks={}
 for label in ['A','B']:
  conf=c['configurations'][label];physical=read_json(verify(c['physical_sources'][label]));assert [v['id'] for v in physical]==[v['id'] for v in p[:len(physical)]]
  x=x0.copy()
  for v in physical:x[ids[v['id']]]=v['xyz_A']
  changed={p[i]['id'] for i in range(len(p)) if not np.array_equal(x[i],x0[i])};assert changed==set(c['perturbation']['changed_source_ids'])
  assert conf['pointcharges']==origin['configurations']['A']['pointcharges'];verify(conf['pointcharges'])
  for v in env:assert np.array_equal(x[ids[v['id']]],v['xyz_A'])
  mapping=read_json(verify(conf['core_mapping']));expected=[]
  for v in mapping:
   if v['kind']=='metal':pos=x[-1];site=v['source_site'];targetid=f"{site['chain']}/{site['resid']}/{site['insertion_code']}/{site['name']}"
   elif v['kind']=='source':pos=x[ids[v['source_id']]]
   else:
    r=x[ids[v['retained_source_id']]];m=x[ids[v['omitted_source_id']]];pos=r+v['length_A']*(m-r)/np.linalg.norm(m-r)
   assert np.max(abs(pos-np.array(v['xyz_A'])))<1e-9;expected.append(pos)
  for metal in ['La','Dy']:
   ep=conf['endpoints'][metal];old=origin['configurations']['A']['endpoints'][metal]
   for key in ['charge','multiplicity','physical_spin_2S','all_electron_count','oxidation_state']:assert ep[key]==old[key]
   lines=verify(ep['xyz']).read_text().splitlines();coords=np.array([[float(v) for v in line.split()[1:4]] for line in lines[2:] if line.strip()]);assert np.max(abs(coords-expected))<1e-9
  checks[label]={'changed_source_ids':sorted(changed),'core_physical_mapping':True,'external_field_unchanged':True,'states_unchanged':True};xs[label]=x
 return c,ledger,artifacts,p,xs,targetid,checks,lp

def main():
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('one shared CPU allocation required')
 out=W/'glu91_classical_20260928';out.mkdir(exist_ok=False);start=time.monotonic();c,ledger,artifacts,p,xs,targetid,checks,lp=prepare();rr={};calls=0
 physicalids=[v['id'] if v['id']!='target_EF3' else targetid for v in p]
 for metal in ['La','Dy']:
  system,description=build_system(ledger,artifacts,metal);integ=mm.VerletIntegrator(.001);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'))
  for label,x in xs.items():
   ctx.setPositions(x*.1);energies={};gradients={}
   for group,name in enumerate(GROUPS):
    state=ctx.getState(getEnergy=True,getForces=True,groups={group});calls+=1;energies[name]=state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole);gradients[name]=-np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom))
   result={'source_ids':physicalids,'coordinates_A':x,'quantity':'gradient','units':'kcal/mol/angstrom','components_kcal_mol':energies,'component_gradients_kcal_mol_A':gradients,'gradient_kcal_mol_A':sum(gradients.values()),'total_kcal_mol':sum(energies.values()),'metal':metal,'configuration':label,'target_source_id':targetid,'QM_MM_Coulomb_included':False,'C4_induction_included':False,'ledger':record(lp),'inputs':record(W/'glu91_motion_v1/INPUTS.json'),'endpoint_state':c['configurations'][label]['endpoints'][metal],'builder_description':description}
   write_new(out/f'{metal}_{label}.json',result);rr[metal+'_'+label]=result
  del ctx;del integ
 common={label:{name:bool(abs(rr['La_'+label]['components_kcal_mol'][name]-rr['Dy_'+label]['components_kcal_mol'][name])<1e-9 and np.max(abs(rr['La_'+label]['component_gradients_kcal_mol_A'][name]-rr['Dy_'+label]['component_gradients_kcal_mol_A'][name]))<1e-9) for name in GROUPS[:3]} for label in xs}
 response={metal:{k:rr[metal+'_B']['components_kcal_mol'][k]-rr[metal+'_A']['components_kcal_mol'][k] for k in GROUPS} for metal in ['La','Dy']}
 result={'status':'complete','inputs':record(W/'glu91_motion_v1/INPUTS.json'),'ledger':record(lp),'endpoints':{key:record(out/f'{key}.json') for key in rr},'B_minus_A_components_kcal_mol':response,'B_minus_A_total_kcal_mol':{m:sum(v.values()) for m,v in response.items()},'La_response_minus_Dy_response_kcal_mol':sum(response['La'].values())-sum(response['Dy'].values()),'Dy_minus_La_per_geometry_kcal_mol':{label:rr['Dy_'+label]['total_kcal_mol']-rr['La_'+label]['total_kcal_mol'] for label in xs},'common_component_checks':common,'preparation_checks':checks,'component_queries':calls,'QM_calls':0,'optimization_steps':0,'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':time.monotonic()-start,'peakRSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'implementation':record(__file__),'plan':record(Path(__file__).with_name('PLAN.md')),'full_hybrid_qualified':False}
 write_new(out/'RESULT.json',result);print((out/'RESULT.json').read_text(),end='')
if __name__=='__main__':main()
