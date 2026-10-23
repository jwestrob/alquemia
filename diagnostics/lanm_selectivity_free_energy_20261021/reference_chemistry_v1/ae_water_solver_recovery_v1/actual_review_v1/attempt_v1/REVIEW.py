#!/usr/bin/env python3
"""Saved-only actual17648 admission and matched four-cell force comparison."""
import datetime as dt
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
BASE = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1'
WORK = BASE/'ae_water_solver_recovery_v1'


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def save(p,v):
    with Path(p).open('x') as f:
        json.dump(v,f,indent=2,allow_nan=False);f.write('\n')


def main():
    start,cpu=time.monotonic(),time.process_time()
    out=WORK/'actual_review_v1';out.mkdir(parents=True,exist_ok=False)
    comp=WORK/'paired_comparison_v1';comp.mkdir(parents=True,exist_ok=False)
    checks=[];pins={}
    def check(ok,why):
        if not ok:raise ValueError(why)
        checks.append(why)
    def verify(p):
        check(pin(p['path'])==p,'Pin '+p['path']);pins[p['path']]=p;return Path(p['path'])
    def read(p):return json.loads(verify(p).read_text())
    def module(p,name):
        s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
    path=WORK/'run_v1/execution/RESULT.json';result=read(pin(path))
    check(result['status']=='recovery_admitted' and result['job_id']=='17648' and all(result['checks'].values()),'Actual complete recovery admitted')
    wrapper=module(result['wrapper'],'saved_recovery_wrapper');release=read(result['release']);manifest=read(result['manifest'])
    for key in ('config','implementation','independent_review','preparation'):verify(release[key])
    verify(release['root_release']['declaration'])
    check(release['manifest']==result['manifest'] and release['implementation']==result['wrapper'] and release['native_invocations']==1 and release['MPI_ranks']==4 and release['CPUs_per_task']==1,'Exact finite released recovery')
    cfg,native,panel,cell=wrapper.load(verify(release['config']))
    raw=read(result['native_result']);directory=Path(result['native_result']['path']).parent
    begin=json.loads((directory/'START.json').read_text());pins[str(directory/'START.json')]=pin(directory/'START.json')
    check(begin['cell']==cell and begin['job_id']=='17648' and begin['ranks']==4 and begin['memory_MiB']==16000 and begin['MaxCore_MiB_per_rank']==3000,'Actual four-rank/default grant and exact cell')
    check(set(begin['thread_environment'].values())=={'1'} and begin['command']==[panel['orca']['path'],'endpoint.inp','--bind-to none'],'Actual native command and one thread/rank')
    for k in ('input','baseline_at_actual_resources','release','manifest','native_helper','panel'):verify(begin[k])
    check(begin['baseline_at_actual_resources']['sha256']==cfg['original_input']['sha256'],'Original baseline deck remains exact at actual resources')
    base,changed=wrapper.recovery_deck(native,panel,cell,4,3000)
    check((directory/'endpoint.inp').read_text()==changed and changed.replace(wrapper.INSERT,wrapper.ORIGINAL)==base,'Sole actual QuadRegionStart input delta')
    check((directory/'geometry.xyz').read_bytes()==verify(cell['geometry']['xyz']).read_bytes(),'Unchanged actual geometry')
    replay=native.collect_one(panel,cell,directory,raw['returncode'])
    check(all(raw[k]==v for k,v in replay.items()) and all(replay['checks'].values()),'Exact original native energy/gradient replay')
    recovery=wrapper.admit(cfg,native,panel,cell,directory,raw)
    check(all(result[k]==v for k,v in recovery.items()),'Exact additional recovery admission replay')
    text=verify(raw['output']).read_text();verify(raw['engrad']);verify(raw['input'])
    # Independent macro-line extraction: micro residuals never determine admission.
    macros=[]
    for i,line in enumerate(text.splitlines(),1):
        t=line.split()
        if 'MAcro)' in line and len(t)>=4 and t[0].isdigit():
            macros.append({'line':i,'iteration':int(t[0]),'energy_hartree':float(t[1]),'norm':float(t[2]),'mode':'TRAH' if '(TRAH' in line else 'NR'})
    check(len(macros)==40 and macros[-1]['iteration']==39 and macros[-1]['norm']==4.543915e-7 and macros[-1]['mode']=='NR','Independent actual final macro residual')
    check(macros[-1]['norm']<=2e-6 and 'SCF CONVERGED AFTER 534 CYCLES' in text and 'ORCA TERMINATED NORMALLY' in text,'Native convergence and active unchanged threshold')
    check('1.000e-06' in next(l for l in text.splitlines() if 'NR start threshold (gradient norm)' in l),'Executed solver entry threshold')
    check((directory/'endpoint.err').stat().st_size==0,'Empty native stderr')
    spin=float(re.findall(r'Expectation value of <S\*\*2>\s*:\s*([-+0-9.eE]+)',text)[-1])
    population=text.split('MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS',1)[1].split('MULLIKEN REDUCED',1)[0].strip()
    original=read(cfg['original_result'])
    check(original['status']=='failed_collection' and not original.get('energy_hartree') and not original.get('forces_kcal_mol_A'),'Original17260 failure remains unadmitted')
    for k in ('physical_state_qualified','method_accuracy_qualified','force_accuracy_qualified'):check(result[k] is False,'Retained false '+k)
    rows={};evidence={}
    for recipe,metal in [('ae_pbe0','La'),('ae_pbe0','Dy'),('ecp_pbe0','La'),('ecp_pbe0','Dy')]:
        cid=recipe+'__water_central_'+metal;c=next(v for v in panel['cells'] if v['cell_id']==cid)
        rp=Path(result['native_result']['path']) if recipe=='ae_pbe0' and metal=='Dy' else BASE/'pilot_v1/tasks'/cid/'execution/RESULT.json'
        r=read(pin(rp));parsed=native.collect_one(panel,c,rp.parent,r['returncode'])
        check(r['status']=='complete' and all(parsed['checks'].values()) and all(r[k]==v for k,v in parsed.items()),'Exact matched comparator raw replay '+cid)
        for k in ('output','engrad','input','geometry'):verify(r[k])
        check(c['geometry']['positions_A']==cell['geometry']['positions_A'] and c['charge']==3 and c['geometry']['atoms']==4,'Identical matched comparator coordinates/charge '+cid)
        check(c['physical_multiplicity']==(1 if metal=='La' else 6) and c['multiplicity']==(6 if recipe=='ae_pbe0' and metal=='Dy' else 1),'Preserved physical/effective spin '+cid)
        check(len(r['forces_kcal_mol_A'])==4 and all(math.isfinite(v) for f in r['forces_kcal_mol_A'] for v in f),'Finite full comparator gradient '+cid)
        rt=Path(r['output']['path']).read_text();scf=wrapper.scf_evidence(rt)
        check(scf['native_scf_converged'] and scf['orbital_residual_within_TolG'],'Comparator native orbital threshold '+cid)
        corrections=[float(v) for v in re.findall(r'Exchange energy change after final integration\s*:\s*([-+0-9.eE]+)',rt)]
        dispersion=[float(v) for v in re.findall(r'Dispersion correction\s+([-+0-9.eE]+)',rt)]
        rows[cid]=r
        evidence[cid]={'result':pin(rp),'job_id':r['job_id'],'recipe':recipe,'state':{k:c[k] for k in ('charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')},'energy_hartree_for_lineage_only':r['energy_hartree'],'input':r['input'],'output':r['output'],'engrad':r['engrad'],'geometry':c['geometry'],'native_checks':r['checks'],'SCF':scf,'final_COSX_energy_changes_hartree':corrections,'dispersion_energy_hartree':dispersion,'force_accuracy_qualified':False,'physical_state_qualified':False}
    x=cell['geometry']['positions_A'];delta=[x[1][i]-x[0][i] for i in range(3)];distance=math.sqrt(sum(v*v for v in delta));u=[v/distance for v in delta]
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    def subtract(a,b):return [[x-y for x,y in zip(aa,bb)] for aa,bb in zip(a,b)]
    def rigid(f):
        center=[sum(x[i][j] for i in (1,2,3))/3 for j in range(3)];net=[sum(f[i][j] for i in (1,2,3)) for j in range(3)];torque=[0.,0.,0.]
        for i in (1,2,3):
            a,b,c=[x[i][j]-center[j] for j in range(3)];p,q,r=f[i]
            torque[0]+=b*r-c*q;torque[1]+=c*p-a*r;torque[2]+=a*q-b*p
        return {'ion_force':f[0],'water_net_force':net,'water_centroid_torque_over_1A':torque}
    def metrics(f):
        b=rigid(f);a=b['ion_force'];other=b['water_net_force']+b['water_centroid_torque_over_1A']
        return {'full4_force_difference_kcal_mol_A':f,'atom_vector_RMS_kcal_mol_A':math.sqrt(sum(v*v for row in f for v in row)/4),'max_atomic_component_kcal_mol_A':max(abs(v) for row in f for v in row),'body_response_difference':b,'category_balanced_rigid_RMS_kcal_mol_A':math.sqrt(sum(v*v for v in a)/6+sum(v*v for v in other)/12),'max_body_component_kcal_mol_A':max(abs(v) for v in a+other)}
    slopes={};forces={}
    for cid,r in rows.items():
        forces[cid]=r['forces_kcal_mol_A'];body=rigid(forces[cid]);slope=-dot(body['water_net_force'],u)
        slopes[cid]={'outward_rigid_water_dE_dR_kcal_mol_A':slope,'metal_force_outward_projection_kcal_mol_A':dot(forces[cid][0],u),'translation_closure_kcal_mol_A':slope-dot(forces[cid][0],u),'rigid_body_response':body}
    method_difference={m:metrics(subtract(forces['ae_pbe0__water_central_'+m],forces['ecp_pbe0__water_central_'+m])) for m in ('La','Dy')}
    paired={recipe:subtract(forces[recipe+'__water_central_Dy'],forces[recipe+'__water_central_La']) for recipe in ('ae_pbe0','ecp_pbe0')}
    paired_slopes={recipe:slopes[recipe+'__water_central_Dy']['outward_rigid_water_dE_dR_kcal_mol_A']-slopes[recipe+'__water_central_La']['outward_rigid_water_dE_dR_kcal_mol_A'] for recipe in ('ae_pbe0','ecp_pbe0')}
    comparison={'status':'saved_matched_force_comparison_complete','implementation':pin(__file__),'recovery_result':pin(path),'source_cells':evidence,'radial_convention':'u=(O-metal)/distance; translate all water atoms outward; dE/dR=-sum_water F dot u','radial_unit_vector':u,'metal_O_distance_A':distance,'radial_slopes':slopes,'AE_minus_ECP_by_metal':method_difference,'paired_Dy_minus_La_force_kcal_mol_A':paired,'paired_radial_slopes':paired_slopes,'paired_slope_AE_minus_ECP_kcal_mol_A':paired_slopes['ae_pbe0']-paired_slopes['ecp_pbe0'],'paired_force_AE_minus_ECP':metrics(subtract(paired['ae_pbe0'],paired['ecp_pbe0'])),'comparison_is_accuracy_test':False,'balanced_binding_energy_available':False,'raw_different_element_or_recipe_energy_contrasts_used':False,'confounding':['frozen-f ECP versus explicit-electron scalar-ZORA Hamiltonians','different metal and ligand basis/auxiliary sets','effective singlet versus physical sextet for Dy','unqualified electronic root and spin-orbit omission','analytic AE force has no independent displacement qualification; final COSX correction retained'],'new_native_calls':0,'new_native_programs':0}
    save(comp/'RESULT.json',comparison)
    blob=(WORK/'run_v1/completion_event.json').read_bytes();(out/'completion_event.snapshot.json').write_bytes(blob);sched=json.loads(blob)['records']['17648']
    check(sched['state']=='COMPLETED' and sched['exit_code']=='0:0' and int(sched['allocated_cpus'])==4,'Actual scheduler successful completion')
    review={'status':'actual_recovery_raw_admission_passed','utc':dt.datetime.now(dt.timezone.utc).isoformat(),'implementation':pin(__file__),'actual_result':pin(path),'native_result':result['native_result'],'native_output':result['native_output'],'release':result['release'],'manifest':result['manifest'],'original_failure':cfg['original_result'],'completion_snapshot':pin(out/'completion_event.snapshot.json'),'comparison':pin(comp/'RESULT.json'),'checks':checks,'check_count':len(checks),'verified_pins':list(pins.values()),'SCF':result['SCF_evidence'],'independent_macro_rows':macros,'spin_squared':spin,'spin_squared_ideal':8.75,'Mulliken_population_raw':population,'native_cost':{'scheduler_elapsed_seconds':int(sched['elapsed_seconds']),'allocated_CPUs':int(sched['allocated_cpus']),'allocated_CPU_seconds':int(sched['elapsed_seconds'])*int(sched['allocated_cpus']),'wrapper_native_wall_seconds':result['elapsed_seconds'],'wrapper_wall_times_four_CPUs':result['allocated_cpu_seconds'],'measured_process_CPU_seconds':None,'measured_CPU_note':'Not present in bound completion/native receipts; allocated CPU time is not measured utilization.'},'new_native_calls':0,'local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'physical_state_qualified':False,'method_accuracy_qualified':False,'force_accuracy_qualified':False,'proactive_terminal_review_not_event_ACK':True}
    save(out/'REVIEW.json',review)
    print(json.dumps({'review':pin(out/'REVIEW.json'),'comparison':pin(comp/'RESULT.json'),'checks':len(checks),'paired_slopes':paired_slopes,'paired_slope_difference':comparison['paired_slope_AE_minus_ECP_kcal_mol_A'],'cost':review['native_cost'],'wall':review['wall_seconds'],'CPU':review['CPU_seconds']}))


if __name__=='__main__':main()
