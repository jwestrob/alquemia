"""Verify native vacuum/GK field-energy identities using real framework controls."""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_tinker_framework_solver import (ENERGY_NAMES,check_parameters,native_call,
    parse_parameters,validate as validate_parent)

PROTOCOL='Tinker26_2_native_direct_field_accounting_v1'
TOL=dict(accounting_kcal_mol=1e-7,rigid_field_e_A2=1e-8)


def build(parent_software,source,output):
    parent=read_json(parent_software)
    if parent['returncode']!=0:raise InvalidArtifact('parent software failed')
    library=verify(parent['library']);verify(parent['executable'])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    target=root/Path(source).name;shutil.copyfile(source,target)
    cmd=[parent['command'][0],'-O2','-g','-fno-fast-math','-fopenmp',
         '-I'+str(library.parent/'mod'),str(target),str(library),
         '-lfftw3_threads','-lfftw3','-o',str(root/'native_field_accounting')]
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (root/'build.log').open('x') as f:
        p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result=dict(command=cmd,returncode=p.returncode,source=record(target),library=parent['library'],
        parent_software=record(parent_software),compiler=parent['compiler'],log=record(root/'build.log'),
        wall_seconds=time.monotonic()-start,
        child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        peak_child_RSS_KiB=after.ru_maxrss,new_energy_calls=0)
    if not p.returncode:result['executable']=record(root/'native_field_accounting')
    write_new(root/'receipt.json',result)
    return result


def task_key(task,m):
    return cache_key(dict(task={k:v for k,v in task.items() if k!='cache_key'},
        protocol=m['protocol'],software=m['software'],implementation=m['implementation'],
        parent=m['parent'],native_parameters=m['native_parameters'],tolerances=m['tolerances'],threads=m['threads']))


def prepare(parent_manifest,software,plan,output):
    start=time.monotonic();cpu=time.process_time()
    parent=validate_parent(parent_manifest);sw=read_json(software)
    if sw['returncode']:raise InvalidArtifact('accounting frontend unavailable')
    for k in ('executable','source','library','parent_software','log'):verify(sw[k])
    if sw['library']!=read_json(verify(parent['software']))['library']:
        raise InvalidArtifact('different native library')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    impl=root/'implementation';impl.mkdir()
    for pin in parent['implementation']:
        p=verify(pin);shutil.copyfile(p,impl/p.name)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    m=dict(protocol=PROTOCOL,parent=record(parent_manifest),software=record(software),
        executable=sw['executable'],native_parameters=parent['native_parameters'],plan=record(plan),
        implementation=[record(p) for p in sorted(impl.glob('*.py'))],threads=64,tolerances=TOL,
        new_energy_calls=6,direct_field_queries=12,new_DFT_calls=0,new_MACE_calls=0,tasks=[])
    for old in parent['tasks']:
        result=Path(parent_manifest).resolve().parent/'tasks'/old['task_id']/'result.json'
        previous=read_json(result)
        if previous['status']!='computed_framework_control' or previous['cache_key']!=old['cache_key']:
            raise InvalidArtifact('parent result not complete/compatible')
        verify(previous['receipt']['log'])
        modes=['fields']+(['static'] if old['variant'] in ('all_standard','frozen_rigid') else [])
        for mode in modes:
            task=copy.deepcopy(old);task.pop('cache_key')
            task.update(task_id=old['task_id']+'_'+mode,parent_task_id=old['task_id'],mode=mode,
                        parent_result=record(result))
            directory=root/'tasks'/task['task_id'];directory.mkdir(parents=True)
            for k in ('xyz','mask','key'):
                src=verify(old[k]);target=directory/src.name
                shutil.copyfile(src,target);task[k]=record(target)
            task['cache_key']=task_key(task,m);m['tasks'].append(task)
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m)
    return m


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['tolerances']!=TOL or m['threads']!=64 or len(m['tasks'])!=18:
        raise InvalidArtifact('changed accounting configuration')
    if [sum(t['mode']==v for t in m['tasks']) for v in ('fields','static')]!=[12,6]:
        raise InvalidArtifact('changed calculation inventory')
    for pin in [m['parent'],m['software'],m['executable'],m['native_parameters'],m['plan']]+m['implementation']:verify(pin)
    sw=read_json(verify(m['software']))
    for k in ('source','library','parent_software','log'):verify(sw[k])
    parent=validate_parent(verify(m['parent']))
    byid={t['task_id']:t for t in parent['tasks']}
    for t in m['tasks']:
        if t['cache_key']!=task_key(t,m):raise InvalidArtifact('changed task cache key')
        for k in ('xyz','mask','key','mapping','parameters','preflight_receipt','parent_result'):verify(t[k])
        old=byid[t['parent_task_id']]
        if any(t[k]['sha256']!=old[k]['sha256'] for k in ('xyz','mask','key','mapping','parameters')):
            raise InvalidArtifact('parent scientific input changed')
        r=read_json(verify(t['parent_result']))
        if r['cache_key']!=old['cache_key']:raise InvalidArtifact('parent result cache mismatch')
        verify(r['receipt']['log'])
    return m


def parse(text,task):
    mode=task['mode'];params=parse_parameters(text)
    if 'ALQUEMIA_ACCOUNTING_COMPLETE '+mode not in text or 'ALQUEMIA_'+mode.upper()+'_BEGIN' not in text:
        raise InvalidArtifact('incomplete accounting execution')
    if 'Induced Dipoles :  Iterations' in text or 'RMS Change' in text:
        raise InvalidArtifact('unexpected induced solve')
    rows={};responses={};energies=None;constants=None;timing=None;state=None
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0] in ('FIELD','RESPONSE'):
            target=rows if f[0]=='FIELD' else responses
            i=int(f[1])
            if i in target or len(f)!=15:raise InvalidArtifact('duplicate/malformed field or response')
            target[i]=list(map(float,f[2:]))
        elif f[0]=='ENERGIES':
            if energies is not None or len(f)!=len(ENERGY_NAMES)+1:raise InvalidArtifact('energy column count')
            energies=dict(zip(ENERGY_NAMES,map(float,f[1:])))
        elif f[0]=='CONSTANTS':constants=list(map(float,f[1:]))
        elif f[0]=='KERNEL_TIMING':timing=list(map(float,f[1:]))
        elif f[0]=='FINAL_STATE':state=[f[1]=='T',int(f[2])]
    n=params['inventory'][0];ids=set(range(1,n+1))
    if set(responses)!=ids or constants is None or timing is None or state is None:
        raise InvalidArtifact('missing response/constants/timing/state')
    if any(r[0]<=0 or any(v!=0 for v in r[1:]) for r in responses.values()):
        raise InvalidArtifact('not the constrained no-response/query state')
    if mode=='fields':
        if set(rows)!=ids or energies is not None or state!=[True,n-len(task['frozen_indices'])]:
            raise InvalidArtifact('invalid field-query state')
        if any(rows[i][0]!=responses[i][0] for i in ids):raise InvalidArtifact('Born radius mismatch')
    elif energies is None or rows or state!=[False,0] or energies['polarization']!=0:
        raise InvalidArtifact('invalid no-response energy')
    if energies is not None:
        if any(energies[k]!=0 for k in ENERGY_NAMES[4:]):raise InvalidArtifact('extra energy component')
        if abs(energies['total']-math.fsum(energies[k] for k in ENERGY_NAMES[1:4]))>TOL['accounting_kcal_mol']:
            raise InvalidArtifact('energy component sum mismatch')
    result=dict(status='computed_accounting_control',parameters=params,fields_e_A2=rows,
        fields_layout='Born_radius_A; vacuum_d[3]; vacuum_p[3]; solvent_d[3]; solvent_p[3]',
        response_eA=responses,energies_kcal_mol=energies,electric=constants[0],dielec=constants[1],
        kernel_wall_seconds=timing[0],kernel_CPU_seconds=timing[1],final_state=state,
        numerical_score=None,full_model_qualified=False)
    cache_key(result)
    return result


def accepted(task,directory):
    for p in sorted(Path(directory).glob('attempt_*/result.json')):
        r=read_json(p)
        if r['cache_key']!=task['cache_key']:raise InvalidArtifact('incompatible cached task')
        verify(r['receipt']['log'])
        if r['status']=='computed_accounting_control':return r,p
    return None


def check_runtime_parameters(params,task,threads):
    # The actual serialized parent coordinates are the invariant. Recomputing
    # R*x with another BLAS/CPU can differ in last bits despite identical input.
    expected=read_json(verify(task['parameters']));expected['inventory'][-1]=threads
    if params!=expected:raise InvalidArtifact('runtime parameters differ from pinned parent')
    lines=verify(task['xyz']).read_text().splitlines()[1:]
    xyz=[[float(v) for v in line.split()[2:5]] for line in lines]
    if xyz!=[a['xyz_A'] for a in params['atoms']]:
        raise InvalidArtifact('runtime coordinates differ from actual serialized input')


def recover(path,output):
    """Reparse the original outputs; zero native executions or changed tolerances."""
    m=validate(path);root=Path(path).resolve().parent;out=Path(output).resolve()
    out.mkdir(parents=True,exist_ok=False);results={}
    for task in m['tasks']:
        if accepted(task,root/'tasks'/task['task_id']):continue
        attempts=sorted((root/'tasks'/task['task_id']).glob('attempt_*/result.json'))
        if not attempts:raise InvalidArtifact('no actual output to recover')
        source=attempts[-1];old=read_json(source)
        if old['failure_reason']!='InvalidArtifact: native geometry changed' or old['receipt']['returncode']!=0:
            raise InvalidArtifact('not the recognized frozen-coordinate validation defect')
        r=parse(verify(old['receipt']['log']).read_text(),task)
        check_runtime_parameters(r['parameters'],task,m['threads'])
        r.update(task_id=task['task_id'],cache_key=task['cache_key'],receipt=old['receipt'],
            job_id=old['job_id'],energy_call_attempted=old['energy_call_attempted'],
            field_query_attempted=old['field_query_attempted'],original_failed_result=record(source),
            recovery_implementation=record(__file__),native_calls_during_recovery=0,
            recovery_reason='Validate identical frozen serialized coordinates and full parameter record instead of regenerating rotation')
        target=out/(task['task_id']+'.json');write_new(target,r);results[task['task_id']]=record(target)
    result=dict(manifest=record(path),implementation=record(__file__),results=results,new_native_calls=0)
    write_new(out/'recovery.json',result);return result


def execute(path):
    m=validate(path);root=Path(path).resolve().parent
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=m['threads']:raise InvalidArtifact('64CPU allocation required')
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for task in m['tasks']:
            directory=root/'tasks'/task['task_id']
            if accepted(task,directory):continue
            attempt=directory/f'attempt_{len(list(directory.glob("attempt_*")))+1:04d}'
            attempt.mkdir()
            receipt=native_call(m['executable'],task,task['mode'],attempt/'native.log',m['threads'])
            text=verify(receipt['log']).read_text()
            r=dict(task_id=task['task_id'],cache_key=task['cache_key'],receipt=receipt,
                job_id=os.environ['SLURM_JOB_ID'],status='failed',
                energy_call_attempted='ALQUEMIA_STATIC_BEGIN' in text,
                field_query_attempted='ALQUEMIA_FIELDS_BEGIN' in text)
            try:
                if receipt['returncode']:raise InvalidArtifact('native process failed')
                r.update(parse(text,task))
                check_runtime_parameters(r['parameters'],task,m['threads'])
            except Exception as error:r.update(status='failed',failure_reason=f'{type(error).__name__}: {error}')
            write_new(attempt/'result.json',r)
            print(json.dumps({k:r.get(k) for k in ('task_id','status','failure_reason')}),flush=True)


def physical_identity(params):
    p=copy.deepcopy(params)
    p['inventory']=p['inventory'][:3];p['global']=p['global'][:5]+p['global'][7:]
    for atom in p['atoms']:atom.pop('response_allowed')
    return p


def collect(path,output,recovery=None):
    m=validate(path);root=Path(path).resolve().parent;computed={};rows={};checks=[];pairs=[]
    recovered=read_json(recovery) if recovery else None
    if recovered:
        verify(recovered['implementation'])
        if recovered['manifest']!=record(path) or recovered['new_native_calls']!=0:
            raise InvalidArtifact('incompatible recovery')
    for task in m['tasks']:
        item=accepted(task,root/'tasks'/task['task_id'])
        if item is None and recovered and task['task_id'] in recovered['results']:
            p=verify(recovered['results'][task['task_id']]);r=read_json(p)
            verify(r['original_failed_result']);verify(r['receipt']['log']);verify(r['recovery_implementation'])
            if r['cache_key']!=task['cache_key'] or r['status']!='computed_accounting_control':
                raise InvalidArtifact('incompatible recovered task')
            check_runtime_parameters(r['parameters'],task,m['threads']);item=(r,p)
        if item:
            r,p=item;computed[task['task_id']]=r
            rows[task['task_id']]=dict(status=r['status'],result=record(p),kernel_wall_seconds=r['kernel_wall_seconds'])
        else:rows[task['task_id']]=dict(status='missing_or_failed')
    for t in (t for t in m['tasks'] if t['mode']=='fields'):
        static_id=t['case_id']+('_frozen_rigid' if t['variant']=='frozen_rigid' else '_all_standard')+'_static'
        if t['task_id'] not in computed or static_id not in computed:continue
        f=computed[t['task_id']];static=computed[static_id];old=read_json(verify(t['parent_result']))
        if physical_identity(f['parameters'])!=physical_identity(static['parameters']):
            raise InvalidArtifact('incompatible static reference geometry/parameters')
        n=f['parameters']['inventory'][0]
        fields=np.array([f['fields_e_A2'][str(i)] for i in range(1,n+1)])
        mu=np.array([old['response_eA'][str(i)] for i in range(1,n+1)])
        noresponse=np.array([static['response_eA'][str(i)] for i in range(1,n+1)])
        if np.max(abs(mu[:,0]-fields[:,0]))>1e-12 or np.max(abs(mu[:,0]-noresponse[:,0]))>1e-12:
            raise InvalidArtifact('query/parent/reference cavity changed')
        frozen=np.array(t['frozen_indices'],dtype=int)-1
        if len(frozen) and np.any(mu[frozen,1:]!=0):raise InvalidArtifact('old source dipoles not frozen')
        factor=-.5*f['electric']/f['dielec']
        vac=factor*math.fsum((mu[:,1:4]*fields[:,4:7]).ravel())
        gk=factor*math.fsum((mu[:,7:10]*fields[:,10:13]).ravel())
        vac_dual=factor*math.fsum((mu[:,4:7]*fields[:,1:4]).ravel())
        gk_dual=factor*math.fsum((mu[:,10:13]*fields[:,7:10]).ravel())
        e=old['energies_kcal_mol'];s=static['energies_kcal_mol']
        residual=dict(vacuum=vac-e['polarization'],solvent=gk-(e['polarization']+e['solvation_with_nonpolar']-s['solvation_with_nonpolar']),
                      permanent=e['permanent']-s['permanent'])
        row=dict(task_id=t['parent_task_id'],vacuum_contraction_kcal_mol=vac,solvent_contraction_kcal_mol=gk,
            native_mutual_components_kcal_mol=e,native_no_response_components_kcal_mol=s,
            accounting_residuals_kcal_mol=residual,
            d_p_reciprocity_difference_kcal_mol=dict(vacuum=vac_dual-vac,solvent=gk_dual-gk),
            pass_=all(abs(v)<=TOL['accounting_kcal_mol'] for v in residual.values()))
        pairs.append(row);checks.append(dict(name=t['parent_task_id']+'_accounting',pass_=row['pass_']))
    for case in sorted({t['case_id'] for t in m['tasks']}):
        a=case+'_frozen_tight_fields';b=case+'_frozen_rigid_fields'
        if a not in computed or b not in computed:continue
        t=next(t for t in m['tasks'] if t['task_id']==b);n=computed[a]['parameters']['inventory'][0]
        av=np.array([computed[a]['fields_e_A2'][str(i)][1:] for i in range(1,n+1)]).reshape(n,4,3)
        bv=np.array([computed[b]['fields_e_A2'][str(i)][1:] for i in range(1,n+1)]).reshape(n,4,3)
        delta=float(np.max(abs(av@np.array(t['rotation']).T-bv)))
        checks.append(dict(name=case+'_rigid_field',max_component_error_e_A2=delta,pass_=delta<=TOL['rigid_field_e_A2']))
    attempts=[read_json(p) for p in root.glob('tasks/*/attempt_*/result.json')]
    result=dict(protocol=PROTOCOL,manifest=record(path),recovery=record(recovery) if recovery else None,complete=len(computed)==18,
        gates_pass=len(checks)==15 and all(c['pass_'] for c in checks),checks=checks,comparisons=pairs,tasks=rows,
        actual_energy_calls=sum(r['energy_call_attempted'] for r in attempts),
        actual_field_queries=sum(r['field_query_attempted'] for r in attempts),
        native_process_wall_seconds=sum(r['receipt']['wall_seconds'] for r in attempts),
        native_process_CPU_seconds=sum(r['receipt']['child_CPU_seconds'] for r in attempts),
        native_kernel_wall_seconds=sum(r.get('kernel_wall_seconds',0) for r in attempts),
        numerical_score=None,full_model_qualified=False,baseline_changed=False)
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('build')
    for name in ('parent-software','source','output'):s.add_argument('--'+name,required=True)
    s=sub.add_parser('prepare')
    for name in ('parent-manifest','software','plan','output'):s.add_argument('--'+name,required=True)
    for name in ('dry-run','execute','collect','recover'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name in ('collect','recover'):s.add_argument('--output',required=True)
        if name=='collect':s.add_argument('--recovery')
    a=p.parse_args()
    if a.command=='build':print(json.dumps(build(a.parent_software,a.source,a.output)))
    elif a.command=='prepare':print(json.dumps({'tasks':len(prepare(a.parent_manifest,a.software,a.plan,a.output)['tasks'])}))
    elif a.command=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_energy_calls':0}))
    elif a.command=='execute':execute(a.manifest)
    elif a.command=='recover':print(json.dumps({'recovered_tasks':len(recover(a.manifest,a.output)['results']),'new_native_calls':0}))
    else:
        r=collect(a.manifest,a.output,a.recovery);print(json.dumps({k:r[k] for k in ('complete','gates_pass','actual_energy_calls','actual_field_queries')}))


if __name__=='__main__':main()
