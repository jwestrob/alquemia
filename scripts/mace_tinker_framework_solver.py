"""Declared twelve-call CPU framework solver control; no La/Ca score."""
from __future__ import annotations
import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import time
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_hybrid import rotation
from mace_tinker_capability import validate as validate_parent

PROTOCOL='Tinker26_2_AMOEBA2018_GK_framework_solver_control_v1'
VARIANTS=('all_standard','frozen_standard','frozen_tight','frozen_rigid')
TOL=dict(energy_kcal_mol=.01,component_kcal_mol=1e-8,frozen_dipole_eA=1e-12)
ENERGY_NAMES=('total','permanent','polarization','solvation_with_nonpolar',
    'bond','angle','stretch_bend','urey','angle_angle','opbend','opdist','improper','imptors',
    'torsion','pitors','strtor','angtor','tortor','vdw','repel','dispersion','charge',
    'charge_dipole','dipole','charge_transfer','reaction_field','metal','restraint','extra')
NONPOLAR=(1.4,.0334,.103,.11,.0135,1.7025,1.3275,.033428,1.,.75,1.056)


def parse_parameters(text):
    rows={}; result={}
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if line.strip().startswith('GK_FLAGS'):
            result['GK_flags']=line.strip()[len('GK_FLAGS'):].split();continue
        if f[0]=='PREP':
            i=int(f[1])
            if i in rows:raise InvalidArtifact('duplicate parameter row')
            if len(f)!=18:raise InvalidArtifact('parameter column count')
            rows[i]=dict(atomic_number=int(f[2]),xyz_A=list(map(float,f[3:6])),
                radius_A=float(f[6]),descreen_radius_A=float(f[7]),scale=float(f[8]),neck=float(f[9]),
                charge_e=float(f[10]),polarizability_A3=float(f[11]),damping_A_half=float(f[12]),
                response_allowed=f[13]=='T',cavity_radius_A=float(f[14]),dispersion_radius_A=float(f[15]),
                dispersion_epsilon_kcal=float(f[16]),dispersion_coefficient=float(f[17]))
        elif f[0]=='ALQUEMIA_PARAMETERS':result['inventory']=list(map(int,f[1:]))
        elif f[0]=='GLOBAL':result['global']=list(map(float,f[1:]))
        elif f[0]=='NONPOLAR':result['nonpolar']=list(map(float,f[1:]))
    if not {'inventory','global','nonpolar','GK_flags'}<=result.keys():
        raise InvalidArtifact('missing native parameters')
    n=result['inventory'][0]
    if set(rows)!=set(range(1,n+1)):raise InvalidArtifact('incomplete native parameter atoms')
    result['atoms']=[rows[i] for i in range(1,n+1)]
    # JSON encoding also rejects nonfinite nested values.
    cache_key(result)
    return result


def check_parameters(params,task,mapping):
    n=len(mapping['system_atom_ids']); frozen=set(task['frozen_indices'])
    if params['inventory'][:4]!=[n,n,n,len(frozen)]:raise InvalidArtifact('native state inventory differs')
    expected_global=[1.,78.3,2.455,.30,.09,task['poleps'],100.,1e12]
    if params['global']!=expected_global:raise InvalidArtifact('unexpected native solver settings')
    if params['nonpolar']!=list(NONPOLAR) or params['GK_flags']!=['GRYCUK','T','T']:
        raise InvalidArtifact('unexpected solvent constants or flags')
    physical={a['id']:a for a in mapping['physical_atoms']};z={'H':1,'C':6,'N':7,'O':8,'S':16}
    source=np.array([physical[i]['xyz_A'] for i in mapping['system_atom_ids']])
    expected=source@np.array(task['rotation']).T+np.array(task['translation_A'])
    for i,(atom,pid,coords) in enumerate(zip(params['atoms'],mapping['system_atom_ids'],expected),1):
        if atom['atomic_number']!=z[physical[pid]['element']]:raise InvalidArtifact('native element changed')
        if atom['xyz_A']!=list(coords):raise InvalidArtifact('native geometry changed')
        if atom['response_allowed']!=(i not in frozen):raise InvalidArtifact('source mask not restored')
        if atom['radius_A']<=0 or atom['descreen_radius_A']<=0 or atom['scale']<0:
            raise InvalidArtifact('invalid native radius/scale')
    return True


def native_call(executable,task,mode,log,threads):
    cmd=[str(verify(executable)),str(verify(task['xyz'])),str(verify(task['mask'])),mode]
    env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_MAX_ACTIVE_LEVELS':'1',
         'OMP_PROC_BIND':'spread','OMP_PLACES':'cores','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    cpu=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
    with Path(log).open('x') as out:
        p=subprocess.run(cmd,cwd=Path(task['xyz']['path']).parent,env=env,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(command=cmd,returncode=p.returncode,wall_seconds=time.monotonic()-start,
        child_CPU_seconds=after.ru_utime+after.ru_stime-cpu.ru_utime-cpu.ru_stime,
        peak_child_RSS_KiB=after.ru_maxrss,threads=threads,log=record(log))


def prepare(parent_manifest,software,plan,output):
    parent=validate_parent(parent_manifest);sw=read_json(software)
    if sw['returncode']!=0:raise InvalidArtifact('solver frontend unavailable')
    verify(sw['executable']);verify(sw['library'])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    impl=root/'implementation';impl.mkdir()
    for name in ('affordable_common.py','mace_hybrid.py','mace_amoeba_capability.py','mace_tinker_capability.py','mace_tinker_framework_solver.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    tasks=[];start=time.monotonic();cpu=time.process_time()
    for old in parent['tasks']:
        mapping=read_json(verify(old['mapping']))
        lines=verify(old['xyz']).read_text().splitlines();source_rows=[s.split() for s in lines[1:]]
        for variant in VARIANTS:
            name=old['case_id']+'_'+variant;directory=root/'tasks'/name;directory.mkdir(parents=True)
            mat=rotation() if variant=='frozen_rigid' else np.eye(3)
            shift=np.array([17.,-23.,9.]) if variant=='frozen_rigid' else np.zeros(3)
            coords=np.array([[float(v) for v in row[2:5]] for row in source_rows])@mat.T+shift
            output_lines=[lines[0]]
            for row,pos in zip(source_rows,coords):
                output_lines.append(' '.join(row[:2]+[format(float(v),'.17g') for v in pos]+row[5:]))
            (directory/'framework.xyz').write_text('\n'.join(output_lines)+'\n')
            frozen=[] if variant=='all_standard' else mapping['frozen_indices']
            (directory/'framework.freeze').write_text(str(len(frozen))+'\n'+''.join(str(i)+'\n' for i in frozen))
            eps=1e-7 if variant in ('frozen_tight','frozen_rigid') else 1e-5
            key=['parameters '+str(verify(parent['parameters'])),'multipoleterm ONLY','polarizeterm',
                 'solvateterm','solvate GK','dielectric 1.0','gk-radius SOLUTE','gkc 2.455',
                 'polarization MUTUAL','polar-iter 100','polar-eps '+str(eps)]
            (directory/'framework.key').write_text('\n'.join(key)+'\n')
            task=dict(task_id=name,case_id=old['case_id'],variant=variant,poleps=eps,
                rotation=mat.tolist(),translation_A=shift.tolist(),frozen_indices=frozen,mapping=old['mapping'],
                xyz=record(directory/'framework.xyz'),mask=record(directory/'framework.freeze'),key=record(directory/'framework.key'))
            receipt=native_call(sw['executable'],task,'preflight',directory/'preflight.log',1)
            write_new(directory/'preflight_receipt.json',receipt)
            text=verify(receipt['log']).read_text()
            if receipt['returncode'] or 'ALQUEMIA_PREFLIGHT_COMPLETE' not in text:
                raise InvalidArtifact(f'parameter preflight failed: {name}')
            params=parse_parameters(text);check_parameters(params,task,mapping)
            write_new(directory/'parameters.json',params);task['parameters']=record(directory/'parameters.json')
            task['preflight_receipt']=record(directory/'preflight_receipt.json')
            tasks.append(task)
    config=dict(protocol=PROTOCOL,parent=record(parent_manifest),software=record(software),executable=sw['executable'],
        native_parameters=parent['parameters'],plan=record(plan),tolerances=TOL,threads=64,
        requested_energy_calls=12,scientific_energy_calls_so_far=0,new_DFT_calls=0,new_MACE_calls=0,
        implementation=[record(p) for p in sorted(impl.glob('*.py'))],tasks=tasks,
        preparation_receipt=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu))
    for task in tasks:
        task['cache_key']=cache_key(dict(task=task,protocol=PROTOCOL,executable=config['executable'],
            native_parameters=config['native_parameters'],implementation=config['implementation'],tolerances=TOL,threads=64))
    write_new(root/'manifest.json',config)
    return config


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['tolerances']!=TOL or len(m['tasks'])!=12:
        raise InvalidArtifact('undeclared solver configuration')
    for pin in [m['parent'],m['software'],m['executable'],m['native_parameters'],m['plan']]+m['implementation']:verify(pin)
    for task in m['tasks']:
        for k in ('mapping','xyz','mask','key','parameters','preflight_receipt'):verify(task[k])
        body={k:v for k,v in task.items() if k!='cache_key'}
        expected=cache_key(dict(task=body,protocol=PROTOCOL,executable=m['executable'],
            native_parameters=m['native_parameters'],implementation=m['implementation'],tolerances=TOL,threads=64))
        if expected!=task['cache_key']:raise InvalidArtifact('changed scientific cache key')
    return m


def parse_energy(text,task):
    params=parse_parameters(text);frozen=set(task['frozen_indices']);response={};energies=None;timing=None
    if 'ALQUEMIA_ENERGY_COMPLETE' not in text or 'not Converged' in text or 'ALQUEMIA_ENERGY_BEGIN' not in text:
        raise InvalidArtifact('incomplete/nonconverged native energy')
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='ENERGIES':
            if energies is not None or len(f)!=len(ENERGY_NAMES)+1:raise InvalidArtifact('native energy column count')
            energies=dict(zip(ENERGY_NAMES,map(float,f[1:])))
        elif f[0]=='RESPONSE':
            index=int(f[1])
            if index in response or len(f)!=15:raise InvalidArtifact('native response column count/duplicate')
            response[index]=list(map(float,f[2:]))
        elif f[0]=='KERNEL_TIMING':timing=list(map(float,f[1:]))
    summaries=re.findall(r'Induced Dipoles\s*:\s*Iterations\s+(\d+)\s+RMS Change\s+([-+0-9.Ee]+)',text)
    if not summaries or any(float(e)+5e-11>task['poleps'] for _,e in summaries):
        raise InvalidArtifact('native residual does not prove requested convergence')
    if energies is None or timing is None or set(response)!=set(range(1,params['inventory'][0]+1)):
        raise InvalidArtifact('missing native energy/response data')
    cache_key(dict(energies=energies,response=response,timing=timing))
    if any(row[0]<=0 for row in response.values()):raise InvalidArtifact('nonpositive Born radius')
    residual=energies['total']-math.fsum(energies[k] for k in ENERGY_NAMES[1:4])
    disabled={k:energies[k] for k in ENERGY_NAMES[4:]}
    max_frozen=max((max(abs(v) for v in response[i][1:]) for i in frozen),default=0.)
    return dict(status='computed_framework_control',energies_kcal_mol=energies,
        parameters=params,response_eA={str(i):v for i,v in response.items()},
        response_layout='Born_radius_A; vacuum_d[3]; vacuum_p[3]; solvent_d[3]; solvent_p[3]',
        component_residual_kcal_mol=residual,disabled_components_zero=all(v==0 for v in disabled.values()),
        max_frozen_dipole_eA=max_frozen,solver_summaries=[dict(iterations=int(i),reported_RMS_Debye=float(e)) for i,e in summaries],
        kernel_wall_seconds=timing[0],kernel_CPU_seconds=timing[1],numerical_score=None,full_model_qualified=False)


def execute(path):
    m=validate(path);root=Path(path).resolve().parent
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=64:raise InvalidArtifact('declared64CPU allocation required')
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for task in m['tasks']:
            directory=root/'tasks'/task['task_id'];out=directory/'result.json'
            if out.exists():
                old=read_json(out)
                if old['cache_key']!=task['cache_key']:raise InvalidArtifact('incompatible previous task')
                verify(old['receipt']['log']);continue
            receipt=native_call(m['executable'],task,'run',directory/'energy.log',64)
            text=verify(receipt['log']).read_text();r=dict(task_id=task['task_id'],cache_key=task['cache_key'],
                status='failed',receipt=receipt,job_id=os.environ['SLURM_JOB_ID'],
                energy_call_attempted='ALQUEMIA_ENERGY_BEGIN' in text)
            try:
                if receipt['returncode']:raise InvalidArtifact('native process failed')
                r.update(parse_energy(text,task))
                expected=read_json(verify(task['parameters']));actual=r['parameters']
                expected['inventory'][-1]=64
                if actual!=expected:raise InvalidArtifact('runtime native parameters differ from preflight')
            except Exception as error:r.update(status='failed',failure_reason=f'{type(error).__name__}: {error}')
            write_new(out,r);print(json.dumps({k:r.get(k) for k in ('task_id','status','failure_reason')}),flush=True)


def collect(path,output):
    m=validate(path);root=Path(path).resolve().parent;rows={};checks=[]
    for task in m['tasks']:
        p=root/'tasks'/task['task_id']/'result.json'
        if not p.exists():rows[task['task_id']]=dict(status='missing');continue
        r=read_json(p)
        if r['cache_key']!=task['cache_key']:raise InvalidArtifact('task result key differs')
        verify(r['receipt']['log']);rows[task['task_id']]={k:v for k,v in r.items() if k not in ('parameters','response_eA')}
        if r['status']=='computed_framework_control':
            checks.extend([dict(name=task['task_id']+'_component',pass_=abs(r['component_residual_kcal_mol'])<=TOL['component_kcal_mol'] and r['disabled_components_zero']),
                           dict(name=task['task_id']+'_source_freeze',pass_=r['max_frozen_dipole_eA']<=TOL['frozen_dipole_eA'])])
    for case in sorted(set(t['case_id'] for t in m['tasks'])):
        a,b,c=(rows[case+'_'+v] for v in ('frozen_standard','frozen_tight','frozen_rigid'))
        if all(r['status']=='computed_framework_control' for r in (a,b,c)):
            for label,left,right in (('refinement',a,b),('rigid',b,c)):
                components={k:right['energies_kcal_mol'][k]-left['energies_kcal_mol'][k] for k in ENERGY_NAMES[:4]}
                checks.append(dict(name=case+'_'+label,component_changes_kcal_mol=components,
                    pass_=abs(components['total'])<=TOL['energy_kcal_mol']))
    result=dict(manifest=record(path),protocol=PROTOCOL,cases=rows,checks=checks,
        complete=all(r['status']=='computed_framework_control' for r in rows.values()),
        gates_pass=len(checks)==30 and all(c['pass_'] for c in checks),
        attempted_energy_calls=sum(r.get('energy_call_attempted',False) for r in rows.values()),
        native_solver_wall_seconds=sum(r.get('kernel_wall_seconds',0) for r in rows.values()),
        numerical_score=None,full_model_qualified=False,baseline_changed=False)
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare')
    for name in ('parent-manifest','software','plan','output'):s.add_argument('--'+name,required=True)
    for name in ('dry-run','execute','collect'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name=='collect':s.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':print(json.dumps({'tasks':len(prepare(a.parent_manifest,a.software,a.plan,a.output)['tasks'])}))
    elif a.command=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_energy_calls':0}))
    elif a.command=='execute':execute(a.manifest)
    else:
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','gates_pass','attempted_energy_calls','native_solver_wall_seconds')}))


if __name__=='__main__':main()
