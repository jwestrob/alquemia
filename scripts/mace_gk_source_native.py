"""Native source-only GK component verification on immutable physical cavities."""
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
from mace_density_multipoles import parse_moments
from mace_tinker_framework_solver import ENERGY_NAMES

PROTOCOL='native_GK_source_only_component_validation_v1'


def frontend(text):
    old="  if (trim(operation)/='static'.and.trim(operation)/='solve'.and.trim(operation)/='identity') error stop 'invalid mode'"
    new="  if (trim(operation)/='source_self'.and.trim(operation)/='empty') error stop 'invalid component mode'"
    marker="  write(*,*) 'ALQUEMIA_PARAMETERS',n,npole,npolar,frozen_count,omp_get_max_threads()"
    if text.count(old)!=1 or text.count(marker)!=1:raise InvalidArtifact('qualified frontend source differs')
    # Existing static branch disables use_polar and starts every induced array at zero.
    zero="  do i=1,n\n    if (.not.frozen(i).or.trim(operation)=='empty') pole(:,i)=0d0\n  end do\n"
    return text.replace(old,new).replace(marker,zero+marker)


def build(parent,output):
    p=read_json(parent);root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    text=frontend(verify(p['frontend']).read_text());f=root/'source_component.f90';f.write_text(text)
    lib=verify(p['library']);cmd=[p['command'][0],'-O2','-g','-fno-fast-math','-fopenmp','-I'+str(lib.parent/'mod'),str(f),str(verify(p['derived_routine'])),str(lib),'-lfftw3_threads','-lfftw3','-o',str(root/'source_component')]
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (root/'build.log').open('x') as out:r=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    res=dict(parent=record(parent),frontend=record(f),command=cmd,returncode=r.returncode,log=record(root/'build.log'),library=p['library'],derived_routine=p['derived_routine'],wall_seconds=time.monotonic()-start,child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
    if not r.returncode:res['executable']=record(root/'source_component')
    write_new(root/'receipt.json',res)
    if r.returncode:raise InvalidArtifact('isolated component frontend failed to compile')
    return res


def task_key(t,m):return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},protocol=m['protocol'],software=m['software'],analysis=m['analysis'],plan=m['plan'],implementation=m['implementation']))


def prepare(analysis,software,plan,output):
    a=read_json(analysis);c=read_json(verify(a['config']));sw=read_json(software)
    if not a['checks_pass'] or sw['returncode']:raise InvalidArtifact('analytical/native prerequisites failed')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    parent=read_json(verify(read_json(verify(next(iter(c['collections'].values()))))['manifest']))
    for n,pin in parent['implementation'].items():shutil.copyfile(verify(pin),impl/n)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    m=dict(protocol=PROTOCOL,analysis=record(analysis),software=record(software),plan=record(plan),tasks=[],threads=64,requested_energy_calls=3*len(c['collections']),implementation={p.name:record(p) for p in impl.glob('*.py')})
    for case,pin in c['collections'].items():
        rr=read_json(verify(pin));pm=read_json(verify(rr['manifest']))
        if read_json(verify(pm['software']))['library']!=sw['library']:raise InvalidArtifact('parent native library differs')
        for state in ('Ca','La','empty'):
            metal='Ca' if state=='empty' else state
            old=next(t for t in pm['static_tasks'] if t['case_id']==case and t['state']==metal and t['variant']=='primary')
            tid=case+'_'+state;d=root/'tasks'/tid;d.mkdir(parents=True)
            t=dict(task_id=tid,case_id=case,state=state,mode='empty' if state=='empty' else 'source_self',frozen_indices=old['frozen_indices'],parent_result=rr['tasks'][old['task_id']]['result'])
            for name in ('xyz','key','mask','overrides'):
                dst=d/Path(old[name]['path']).name;shutil.copyfile(verify(old[name]),dst);t[name]=record(dst)
            t['cache_key']=task_key(t,m);m['tasks'].append(t)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path);a=read_json(verify(m['analysis']));sw=read_json(verify(m['software']));p=read_json(verify(sw['parent']))
    if m['protocol']!=PROTOCOL or m['threads']!=64 or len(m['tasks'])!=m['requested_energy_calls']:raise InvalidArtifact('component scope differs')
    for pin in [m['plan'],*m['implementation'].values(),sw['frontend'],sw['executable'],sw['library'],sw['derived_routine']]:verify(pin)
    if frontend(verify(p['frontend']).read_text())!=verify(sw['frontend']).read_text():raise InvalidArtifact('component frontend transformation differs')
    expected={(case,state) for case in a['cases'] for state in ('Ca','La','empty')}
    if len(m['tasks'])!=len(expected) or {(t['case_id'],t['state']) for t in m['tasks']}!=expected:raise InvalidArtifact('case/state inventory differs')
    for t in m['tasks']:
        if t['cache_key']!=task_key(t,m):raise InvalidArtifact('changed component task')
        for name in ('xyz','key','mask','overrides','parent_result'):verify(t[name])
    return m


def parse(text,t):
    if 'ALQUEMIA_DENSITY_GK_COMPLETE '+t['mode'] not in ' '.join(text.split()):raise InvalidArtifact('incomplete native component')
    if any(x in text for x in ('ALQUEMIA_FIELDS_BEGIN','ALQUEMIA_RESPONSE_BEGIN')):raise InvalidArtifact('unexpected native field or response operation')
    marker='ALQUEMIA_MOMENTS_COMPLETE'
    if marker not in text:raise InvalidArtifact('missing actual native moment export')
    data=parse_moments(text[:text.index(marker)]+marker,None)
    if not data['pass_']:raise InvalidArtifact('native moment rotation failed')
    old=read_json(verify(t['parent_result']));params=copy.deepcopy(old['parameters']);frozen=set(t['frozen_indices'])
    for i,atom in enumerate(params['atoms'],1):
        if t['state']=='empty' or i not in frozen:atom['charge_e']=0.
    if data['parameters']!=params:raise InvalidArtifact('native geometry, cavity or parameters differ')
    for i,(new,prior) in enumerate(zip(data['atoms'],old['moments']),1):
        expect=prior['global_'] if i in frozen and t['state']!='empty' else [0.]*13
        if new['global_']!=expect:raise InvalidArtifact('native source/environment moments differ')
    energies=None;responses={};timing=None;active=None
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='ENERGIES':energies=dict(zip(ENERGY_NAMES,map(float,f[1:])))
        elif f[0]=='RESPONSE':
            i=int(f[1])
            if i in responses or len(f)!=15:raise InvalidArtifact('duplicate or invalid native response row')
            responses[i]=list(map(float,f[2:]))
        elif f[0]=='KERNEL_TIMING':timing=list(map(float,f[1:]))
        elif f[0]=='ACTIVE_TERMS':active=f[1:]
    n=params['inventory'][0]
    if set(responses)!=set(range(1,n+1)) or not energies or timing is None or active!=[t['mode'],'T','T','T']:raise InvalidArtifact('missing static component evidence')
    if energies['polarization'] or any(energies[k] for k in ENERGY_NAMES[4:]) or abs(energies['total']-energies['permanent']-energies['solvation_with_nonpolar'])>1e-7:raise InvalidArtifact('energy components differ')
    err=max(abs(responses[i][0]-old['fields'][str(i)][0]) for i in responses)
    if err>1e-10 or any(x for row in responses.values() for x in row[1:]):raise InvalidArtifact('cavity changed or induced response nonzero')
    result=dict(status='computed',energies=energies,maximum_Born_difference_A=err,kernel_wall_seconds=timing[0],kernel_CPU_seconds=timing[1])
    cache_key(result);return result


def accepted(t,root):
    for p in sorted((root/'tasks'/t['task_id']).glob('attempt_*/result.json')):
        r=read_json(p)
        if r['cache_key']!=t['cache_key']:raise InvalidArtifact('incompatible component cache')
        if r.get('receipt'):verify(r['receipt']['log'])
        if r['status']=='computed':return r,p
    return None


def execute(path,retry=False):
    m=validate(path);root=Path(path).resolve().parent;sw=read_json(verify(m['software']))
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=64:raise InvalidArtifact('declared64CPU allocation required')
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for t in m['tasks']:
            if accepted(t,root):continue
            d=root/'tasks'/t['task_id'];attempts=sorted(d.glob('attempt_*'))
            if attempts and not retry:continue
            out=d/f'attempt_{len(attempts)+1:04d}';out.mkdir();r=dict(task_id=t['task_id'],cache_key=t['cache_key'],status='failed',receipt=None)
            try:
                cmd=[str(verify(sw['executable'])),str(verify(t['xyz'])),str(verify(t['mask'])),str(verify(t['overrides'])),t['mode']]
                start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
                env={**os.environ,'OMP_NUM_THREADS':'64','OMP_MAX_ACTIVE_LEVELS':'1','OMP_PROC_BIND':'spread','OMP_PLACES':'cores','OPENBLAS_NUM_THREADS':'1'}
                with (out/'native.log').open('x') as f:proc=subprocess.run(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env)
                after=resource.getrusage(resource.RUSAGE_CHILDREN)
                r['receipt']=dict(command=cmd,returncode=proc.returncode,wall_seconds=time.monotonic()-start,child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,peak_RSS_KiB=after.ru_maxrss,job_id=os.environ['SLURM_JOB_ID'],log=record(out/'native.log'))
                if proc.returncode:raise InvalidArtifact('native execution failed')
                r.update(parse(verify(r['receipt']['log']).read_text(),t))
            except Exception as exc:r['failure_reason']=str(exc)
            write_new(out/'result.json',r);print(json.dumps({k:r.get(k) for k in ('task_id','status','failure_reason')}),flush=True)


def collect(path,output):
    m=validate(path);root=Path(path).resolve().parent;a=read_json(verify(m['analysis']));rows={};cases={};checks=[]
    for t in m['tasks']:
        pair=accepted(t,root);rows[t['task_id']]=dict(status='unavailable') if pair is None else dict(status='computed',result=record(pair[1]))
    for case,expected in a['cases'].items():
        if any(rows[case+'_'+s]['status']=='unavailable' for s in ('Ca','La','empty')):cases[case]=dict(status='unavailable');continue
        results={s:read_json(verify(rows[case+'_'+s]['result'])) for s in ('Ca','La','empty')}
        values={s:results[s]['energies']['solvation_with_nonpolar']-results['empty']['energies']['solvation_with_nonpolar'] for s in ('Ca','La')}
        values['R']=values['Ca']-values['La'];errors={s:values[s]-expected['endpoints'][s]['source_self_kcal'] for s in ('Ca','La')}
        errors['R']=values['R']-expected['contrasts']['source_self_kcal'];checks.append(dict(case_id=case,errors_kcal=errors,pass_=max(abs(x) for x in errors.values())<=1e-7));cases[case]=dict(status='computed',source_self_kcal=values)
    attempts=[read_json(p) for p in root.glob('tasks/*/attempt_*/result.json')]
    r=dict(protocol=PROTOCOL,manifest=record(path),complete=all(x['status']=='computed' for x in rows.values()),cases=cases,tasks=rows,checks=checks,checks_pass=len(checks)==len(a['cases']) and all(x['pass_'] for x in checks),native_attempts=len(attempts),native_energy_calls=sum('ALQUEMIA_ENERGY_BEGIN' in verify(x['receipt']['log']).read_text() for x in attempts if x.get('receipt')),new_score=None)
    write_new(output,r);return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('build');q.add_argument('--parent',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('prepare')
    for k in ('analysis','software','plan','output'):q.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op=='execute':q.add_argument('--retry-failed',action='store_true')
        if op=='collect':q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='build':r=build(a.parent,a.output);print(r['returncode'])
    elif a.command=='prepare':r=prepare(a.analysis,a.software,a.plan,a.output);print({'tasks':len(r['tasks'])})
    elif a.command=='dry-run':r=validate(a.manifest);print({'tasks':len(r['tasks']),'status':'pass'})
    elif a.command=='execute':execute(a.manifest,a.retry_failed)
    else:r=collect(a.manifest,a.output);print({k:r[k] for k in ('complete','checks_pass','native_energy_calls')})
