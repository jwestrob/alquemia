"""Native field-input replay; no new quantum source or biological score."""
from __future__ import annotations
import argparse
import copy
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
from mace_native_field_accounting import check_runtime_parameters,physical_identity,validate as validate_parent
from mace_tinker_framework_solver import native_call,parse_parameters

PROTOCOL='Tinker26_2_supplied_direct_fields_native_response_replay_v1'
TOL=dict(dipole_eA=1e-9,energy_kcal_mol=1e-7,rigid_dipole_eA=1e-8)


def derive(text):
    start=text.index('      subroutine induce0c\n');end=text.index('      end\n',start)+len('      end\n')
    original=text[start:end]
    if original.count('      call dfield0d (field,fieldp,fields,fieldps)')!=1:
        raise InvalidArtifact('unexpected native direct-field call')
    out=original.replace('      subroutine induce0c\n','      subroutine alquemia_induce0c (source_fields)\n',1)
    out=out.replace('      implicit none\n','      implicit none\n      real*8 source_fields(3,n,4)\n',1)
    out=out.replace('      call dfield0d (field,fieldp,fields,fieldps)',
        '      field(:,:) = source_fields(:,:,1)\n      fieldp(:,:) = source_fields(:,:,2)\n'
        '      fields(:,:) = source_fields(:,:,3)\n      fieldps(:,:) = source_fields(:,:,4)')
    return original,out


def build(parent_software,native_source,frontend,output):
    parent=read_json(parent_software);library=verify(parent['library'])
    if parent['returncode']:raise InvalidArtifact('parent software unavailable')
    text=Path(native_source).read_text();original,derived=derive(text)
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(native_source,root/'original_induce.f')
    (root/'original_induce0c.f').write_text(original);(root/'supplied_induce0c.f').write_text(derived)
    shutil.copyfile(frontend,root/Path(frontend).name)
    compiler=parent['command'][0]
    cmd=[compiler,'-O2','-g','-fno-fast-math','-fopenmp','-I'+str(library.parent/'mod'),
        str(root/Path(frontend).name),str(root/'supplied_induce0c.f'),str(library),
        '-lfftw3_threads','-lfftw3','-o',str(root/'native_field_input')]
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (root/'build.log').open('x') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt=dict(command=cmd,returncode=r.returncode,library=parent['library'],parent_software=record(parent_software),
        original_full_source=record(root/'original_induce.f'),original_routine=record(root/'original_induce0c.f'),
        derived_routine=record(root/'supplied_induce0c.f'),frontend=record(root/Path(frontend).name),
        compiler=parent['compiler'],log=record(root/'build.log'),transformation_implementation=record(__file__),
        wall_seconds=time.monotonic()-start,child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        peak_child_RSS_KiB=after.ru_maxrss,new_response_solves=0)
    # Freeze the transformation implementation as part of the build product.
    shutil.copyfile(__file__,root/Path(__file__).name)
    receipt['transformation_implementation']=record(root/Path(__file__).name)
    if not r.returncode:receipt['executable']=record(root/'native_field_input')
    write_new(root/'receipt.json',receipt);return receipt


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},protocol=PROTOCOL,
        software=m['software'],implementation=m['implementation'],native_parameters=m['native_parameters'],tolerances=TOL,threads=64))


def prepare(parent_manifest,parent_report,software,plan,output):
    start=time.monotonic();cpu=time.process_time()
    parent=validate_parent(parent_manifest);report=read_json(parent_report);sw=read_json(software)
    if not report['complete'] or not report['gates_pass'] or report['manifest']!=record(parent_manifest):
        raise InvalidArtifact('native accounting prerequisite not qualified')
    if sw['returncode']:raise InvalidArtifact('field input frontend unavailable')
    original,derived=derive(verify(sw['original_full_source']).read_text())
    if original!=verify(sw['original_routine']).read_text() or derived!=verify(sw['derived_routine']).read_text():
        raise InvalidArtifact('native routine transformation differs')
    if sw['library']!=read_json(verify(parent['software']))['library']:raise InvalidArtifact('native library changed')
    verify(sw['executable']);verify(sw['frontend'])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for pin in parent['implementation']:
        p=verify(pin);shutil.copyfile(p,impl/p.name)
    for name in ('mace_native_field_accounting.py','mace_native_field_input.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    m=dict(protocol=PROTOCOL,parent=record(parent_manifest),parent_report=record(parent_report),software=record(software),
        executable=sw['executable'],native_parameters=parent['native_parameters'],plan=record(plan),threads=64,
        implementation=[record(p) for p in sorted(impl.glob('*.py'))],tolerances=TOL,tasks=[],
        new_response_solves=6,new_energy_calls=0,new_DFT_calls=0,new_MACE_calls=0)
    for old in parent['tasks']:
        if old['mode']!='fields' or old['variant'] not in ('frozen_tight','frozen_rigid'):continue
        t=copy.deepcopy(old);t.pop('cache_key');t['task_id']=old['parent_task_id']+'_supplied'
        directory=root/'tasks'/t['task_id'];directory.mkdir(parents=True)
        t['field_result']=report['tasks'][old['task_id']]['result']
        static_id=old['case_id']+('_frozen_rigid' if old['variant']=='frozen_rigid' else '_all_standard')+'_static'
        t['static_result']=report['tasks'][static_id]['result']
        f=read_json(verify(t['field_result']));s=read_json(verify(t['static_result']))
        if physical_identity(f['parameters'])!=physical_identity(s['parameters']):raise InvalidArtifact('incompatible static state')
        n=f['parameters']['inventory'][0];lines=[str(n)]
        for i in range(1,n+1):lines.append(str(i)+' '+' '.join(format(v,'.17g') for v in f['fields_e_A2'][str(i)]))
        field_path=directory/'framework.fields';field_path.write_text('\n'.join(lines)+'\n');t['supplied_fields']=record(field_path)
        for k in ('xyz','mask','key'):
            source=verify(t[k]);target=directory/source.name;shutil.copyfile(source,target);t[k]=record(target)
        t['cache_key']=key(t,m);m['tasks'].append(t)
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);return m


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['tolerances']!=TOL or len(m['tasks'])!=6 or m['threads']!=64:
        raise InvalidArtifact('changed response configuration')
    for pin in [m['parent'],m['parent_report'],m['software'],m['executable'],m['native_parameters'],m['plan']]+m['implementation']:verify(pin)
    sw=read_json(verify(m['software']))
    for k in ('library','original_full_source','original_routine','derived_routine','frontend','transformation_implementation','log'):verify(sw[k])
    original,derived=derive(verify(sw['original_full_source']).read_text())
    if original!=verify(sw['original_routine']).read_text() or derived!=verify(sw['derived_routine']).read_text():
        raise InvalidArtifact('unapproved native routine edit')
    parent=validate_parent(verify(m['parent']));old={t['parent_task_id']:t for t in parent['tasks'] if t['mode']=='fields'}
    for t in m['tasks']:
        if t['cache_key']!=key(t,m):raise InvalidArtifact('changed response cache key')
        for k in ('xyz','mask','key','mapping','parameters','field_result','static_result','supplied_fields','parent_result'):verify(t[k])
        if any(t[k]['sha256']!=old[t['parent_task_id']][k]['sha256'] for k in ('xyz','mask','key','parameters','mapping')):
            raise InvalidArtifact('parent scientific input changed')
        for k in ('field_result','static_result','parent_result'):verify(read_json(verify(t[k]))['receipt']['log'])
    return m


def parse(text,t):
    if 'ALQUEMIA_RESPONSE_COMPLETE' not in text or 'not Converged' in text:
        raise InvalidArtifact('incomplete/nonconverged native response')
    params=parse_parameters(text);check_runtime_parameters(params,t,64)
    rows={};timing=None;constants=None
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='RESPONSE':
            i=int(f[1])
            if i in rows or len(f)!=15:raise InvalidArtifact('duplicate/malformed response')
            rows[i]=list(map(float,f[2:]))
        elif f[0]=='KERNEL_TIMING':timing=list(map(float,f[1:]))
        elif f[0]=='CONSTANTS':constants=list(map(float,f[1:]))
    summaries=re.findall(r'Induced Dipoles\s*:\s*Iterations\s+(\d+)\s+RMS Change\s+([-+0-9.Ee]+)',text)
    if not summaries or any(float(e)+5e-11>t['poleps'] for _,e in summaries):
        raise InvalidArtifact('requested native convergence not demonstrated')
    n=params['inventory'][0]
    if set(rows)!=set(range(1,n+1)) or timing is None or constants is None:raise InvalidArtifact('missing response evidence')
    if any(rows[i][0]<=0 for i in rows) or any(v!=0 for i in t['frozen_indices'] for v in rows[i][1:]):
        raise InvalidArtifact('invalid cavity or source dipole changed')
    result=dict(status='computed_response_replay',parameters=params,response_eA=rows,
        electric=constants[0],dielec=constants[1],kernel_wall_seconds=timing[0],kernel_CPU_seconds=timing[1],
        solver_summaries=[dict(iterations=int(i),reported_RMS_Debye=float(e)) for i,e in summaries],
        numerical_score=None,full_model_qualified=False)
    cache_key(result);return result


def accepted(t,directory):
    for p in sorted(Path(directory).glob('attempt_*/result.json')):
        r=read_json(p)
        if r['cache_key']!=t['cache_key']:raise InvalidArtifact('incompatible previous attempt')
        verify(r['receipt']['log'])
        if r['status']=='computed_response_replay':return r,p
    return None


def execute(path):
    m=validate(path);root=Path(path).resolve().parent
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=64:raise InvalidArtifact('64CPU allocation required')
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for t in m['tasks']:
            directory=root/'tasks'/t['task_id']
            if accepted(t,directory):continue
            attempt=directory/f'attempt_{len(list(directory.glob("attempt_*")))+1:04d}';attempt.mkdir()
            receipt=native_call(m['executable'],t,str(verify(t['supplied_fields'])),attempt/'native.log',64)
            text=verify(receipt['log']).read_text();r=dict(task_id=t['task_id'],cache_key=t['cache_key'],status='failed',
                receipt=receipt,job_id=os.environ['SLURM_JOB_ID'],response_solve_attempted='ALQUEMIA_RESPONSE_BEGIN' in text)
            try:
                if receipt['returncode']:raise InvalidArtifact('native process failed')
                r.update(parse(text,t))
            except Exception as error:r.update(status='failed',failure_reason=f'{type(error).__name__}: {error}')
            write_new(attempt/'result.json',r);print(json.dumps({k:r.get(k) for k in ('task_id','status','failure_reason')}),flush=True)


def collect(path,output):
    m=validate(path);root=Path(path).resolve().parent;rows={};checks=[];responses={}
    for t in m['tasks']:
        item=accepted(t,root/'tasks'/t['task_id'])
        if not item:rows[t['task_id']]=dict(status='missing_or_failed');continue
        r,p=item;n=r['parameters']['inventory'][0];mu=np.array([r['response_eA'][str(i)] for i in range(1,n+1)])
        old=read_json(verify(t['parent_result']));f=read_json(verify(t['field_result']));s=read_json(verify(t['static_result']))
        original=np.array([old['response_eA'][str(i)] for i in range(1,n+1)])
        fields=np.array([f['fields_e_A2'][str(i)] for i in range(1,n+1)])
        if np.max(abs(mu[:,0]-fields[:,0]))>1e-12:raise InvalidArtifact('replay cavity changed')
        if r['electric']!=f['electric'] or r['dielec']!=f['dielec']:raise InvalidArtifact('native conversion constant changed')
        c=-.5*r['electric']/r['dielec']
        vac=c*math.fsum((mu[:,1:4]*fields[:,4:7]).ravel());solv=c*math.fsum((mu[:,7:10]*fields[:,10:13]).ravel())
        total=s['energies_kcal_mol']['total']+solv
        errors=dict(max_dipole_eA=float(np.max(abs(mu[:,1:]-original[:,1:]))),
            vacuum_kcal_mol=vac-old['energies_kcal_mol']['polarization'],total_kcal_mol=total-old['energies_kcal_mol']['total'])
        pass_=errors['max_dipole_eA']<=TOL['dipole_eA'] and max(abs(errors[k]) for k in ('vacuum_kcal_mol','total_kcal_mol'))<=TOL['energy_kcal_mol']
        rows[t['task_id']]=dict(status=r['status'],result=record(p),errors=errors,
            vacuum_contraction_kcal_mol=vac,solvent_contraction_kcal_mol=solv,total_kcal_mol=total,
            static_components_kcal_mol=s['energies_kcal_mol'],parent_components_kcal_mol=old['energies_kcal_mol'],
            solver_summaries=r['solver_summaries'])
        checks.append(dict(name=t['task_id']+'_replay',pass_=pass_));responses[t['task_id']]=mu[:,1:].reshape(n,4,3)
    for case in sorted({t['case_id'] for t in m['tasks']}):
        a=case+'_frozen_tight_supplied';b=case+'_frozen_rigid_supplied'
        if a not in responses or b not in responses:continue
        t=next(t for t in m['tasks'] if t['task_id']==b)
        d=float(np.max(abs(responses[a]@np.array(t['rotation']).T-responses[b])))
        e=rows[b]['total_kcal_mol']-rows[a]['total_kcal_mol']
        checks.append(dict(name=case+'_rigid',max_dipole_eA=d,energy_difference_kcal_mol=e,
            pass_=d<=TOL['rigid_dipole_eA'] and abs(e)<=TOL['energy_kcal_mol']))
    attempts=[read_json(p) for p in root.glob('tasks/*/attempt_*/result.json')]
    r=dict(protocol=PROTOCOL,manifest=record(path),complete=len(responses)==6,gates_pass=len(checks)==9 and all(c['pass_'] for c in checks),
        tasks=rows,checks=checks,actual_response_solves=sum(a['response_solve_attempted'] for a in attempts),
        native_kernel_wall_seconds=sum(a.get('kernel_wall_seconds',0) for a in attempts),
        native_process_wall_seconds=sum(a['receipt']['wall_seconds'] for a in attempts),
        native_process_CPU_seconds=sum(a['receipt']['child_CPU_seconds'] for a in attempts),
        numerical_score=None,full_model_qualified=False,baseline_changed=False)
    write_new(output,r);return r


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('build')
    for name in ('parent-software','native-source','frontend','output'):s.add_argument('--'+name,required=True)
    s=sub.add_parser('prepare')
    for name in ('parent-manifest','parent-report','software','plan','output'):s.add_argument('--'+name,required=True)
    for name in ('dry-run','execute','collect'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name=='collect':s.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='build':
        r=build(a.parent_software,a.native_source,a.frontend,a.output);print(json.dumps(r));raise SystemExit(r['returncode'])
    elif a.command=='prepare':print(json.dumps({'tasks':len(prepare(a.parent_manifest,a.parent_report,a.software,a.plan,a.output)['tasks'])}))
    elif a.command=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_response_solves':0}))
    elif a.command=='execute':execute(a.manifest)
    else:
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','gates_pass','actual_response_solves')}))


if __name__=='__main__':main()
