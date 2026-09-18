"""Native operator replay for a frozen AMOEBA reaction-field transfer."""
from pathlib import Path
import argparse
import fcntl
import json
import math
import os
import resource
import shutil
import subprocess
import time
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new
from mace_density_multipoles import parse_moments

PROTOCOL = 'native_AMOEBA_GK_frozen_response_functional_check_v1'
CASES = ('GGR_2FW0', 'GGR_2FVY')
TOL = dict(input_native=1e-12, born_A=1e-10, static_kcal=1e-7,
           reaction_field=1e-10, cross_kcal=1e-7, stationary_kcal=1e-5,
           residual_Debye=1e-7, reciprocity_kcal=1e-7)

TAIL = r'''
  allocate(expected_born(n),saved(13,n),mud(3,n),mup(3,n),drive(3,n,2))
  allocate(vd(3,n),vp(3,n),sd(3,n),sp(3,n),fd(3,n),fp(3,n),fs(3,n),fps(3,n))
  open(newunit=mask_unit,file=trim(field_path),status='old',action='read',iostat=ios)
  if (ios/=0) error stop 'missing frozen response'
  read(mask_unit,*,iostat=ios) count_input
  if (ios/=0.or.count_input/=n) error stop 'frozen response count differs'
  do i=1,n
    read(mask_unit,*,iostat=ios) index,expected_born(i),mud(:,i),mup(:,i),drive(:,i,1),drive(:,i,2)
    if (ios/=0.or.index/=i) error stop 'frozen response order differs'
    if (.not.all(ieee_is_finite([expected_born(i),mud(:,i),mup(:,i),drive(:,i,1),drive(:,i,2)]))) &
      error stop 'nonfinite response input'
    if (frozen(i).and.any([mud(:,i),mup(:,i)]/=0d0)) error stop 'source response is not frozen'
  end do
  read(mask_unit,*,iostat=ios) index
  if (ios>=0) error stop 'extra frozen response data'
  close(mask_unit)
  call system_clock(start_clock,clock_rate)
  call cpu_time(begin_cpu)
  call replica(0d0)
  call born
  if (any(abs(rborn-expected_born)>1d-10)) error stop 'frozen cavity differs'
  switch_mode='MPOLE'
  call switch(switch_mode)
  write(*,*) 'MUTUAL_CUTOFF_SQUARED',off2
  if (off2/=mpolecut**2) error stop 'native field cutoff differs'
  uind=mud; uinp=mup; uinds=mud; uinps=mup
  write(*,*) 'ALQUEMIA_MUTUAL_OPERATOR_BEGIN'
  flush(6)
  call ufield0d(vd,vp,sd,sp)
  do i=1,n
    write(*,*) 'MUTUAL',i,rborn(i),polarity(i),vd(:,i),vp(:,i),sd(:,i),sp(:,i)
  end do
  write(*,*) 'ALQUEMIA_PERMANENT_OPERATOR_BEGIN'
  call dfield0d(fd,fp,fs,fps)
  do i=1,n
    write(*,*) 'REACTION',i,fs(:,i)-fd(:,i),fps(:,i)-fp(:,i)
  end do
  saved=rpole
  uind=0d0; uinp=0d0; uinds=0d0; uinps=0d0
  use_polar=.false.
  call enp(ecav,edisp)
  write(*,*) 'NONPOLAR_ENERGY',ecav,edisp
  do k=1,3
    rpole=saved
    if (k==2) rpole(2:4,:)=saved(2:4,:)+0.5d0*(mud+mup)
    if (k==3) then
      rpole=0d0
      rpole(2:4,:)=0.5d0*(mud-mup)
    end if
    es=0d0
    write(*,*) 'ALQUEMIA_STATIC_GK_BEGIN',k
    flush(6)
    call egk
    write(*,*) 'STATIC_GK',k,es
  end do
  call cpu_time(end_cpu)
  call system_clock(end_clock)
  write(*,*) 'KERNEL_TIMING',real(end_clock-start_clock,8)/real(clock_rate,8),end_cpu-begin_cpu
  write(*,*) 'ALQUEMIA_FROZEN_RESPONSE_COMPLETE',n
  call final
end program alquemia_tinker_density_gk
'''


def frontend(text):
    marker = '  call system_clock(start_clock,clock_rate)'
    if text.count(marker) != 1:
        raise InvalidArtifact('unexpected qualified frontend timing marker')
    text = text[:text.index(marker)]
    text = text.replace('  use chgpot\n', '  use chgpot\n  use units, only: debye\n', 1)
    text = text.replace('  use limits\n', '  use limits\n  use shunt, only: off2\n', 1)
    text += "  write(*,*) 'NATIVE_DEBYE',debye\n"
    text = text.replace("if (command_argument_count()<4.or.command_argument_count()>5)",
                        "if (command_argument_count()/=5)")
    old = "  if (trim(operation)/='static'.and.trim(operation)/='solve'.and.trim(operation)/='identity') error stop 'invalid mode'"
    if text.count(old) != 1:
        raise InvalidArtifact('unexpected native frontend modes')
    text = text.replace(old, "  if (trim(operation)/='replay') error stop 'invalid mode'")
    text = text.replace("  if (trim(operation)=='solve') then", "  if (trim(operation)=='replay') then", 1)
    declaration = ('  real*8, allocatable :: saved(:,:),mud(:,:),mup(:,:),drive(:,:,:)\n'
                   '  real*8, allocatable :: vd(:,:),vp(:,:),sd(:,:),sp(:,:)\n'
                   '  real*8 :: ecav,edisp\n')
    text = text.replace('  logical, allocatable :: frozen(:)\n',
                        '  logical, allocatable :: frozen(:)\n' + declaration, 1)
    return text + TAIL


def build(parent, output):
    p = read_json(parent); root = Path(output).resolve(); root.mkdir(parents=True, exist_ok=False)
    lib = verify(p['library']); f = root/'frozen_response.f90'
    f.write_text(frontend(verify(p['frontend']).read_text()))
    shutil.copyfile(__file__, root/Path(__file__).name)
    cmd = [p['command'][0], '-O2', '-g', '-fno-fast-math', '-fopenmp', '-I'+str(lib.parent/'mod'),
           str(f), str(lib), '-lfftw3_threads', '-lfftw3', '-o', str(root/'frozen_response')]
    start = time.monotonic(); before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with (root/'build.log').open('x') as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    r = dict(parent=record(parent), frontend=record(f), library=p['library'], command=cmd,
             returncode=proc.returncode, log=record(root/'build.log'), implementation=record(root/Path(__file__).name),
             wall_seconds=time.monotonic()-start,
             CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
    if not proc.returncode: r['executable'] = record(root/'frozen_response')
    write_new(root/'receipt.json', r)
    if proc.returncode: raise InvalidArtifact('native replay frontend compilation failed')
    return r


def task_key(t, m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},
                         protocol=m['protocol'], software=m['software'], plan=m['plan'],
                         implementation=m['implementation'], tolerances=m['tolerances']))


def source_arrays(s, r):
    n = s['parameters']['inventory'][0]
    if s['parameters'] != r['parameters'] or s['moments'] != r['moments']:
        raise InvalidArtifact('static/response physical inputs differ')
    mu = np.array([r['response'][str(i)] for i in range(1,n+1)])
    rows = verify(r['supplied_fields']).read_text().splitlines()
    if int(rows[0]) != n or len(rows) != n+1:
        raise InvalidArtifact('supplied field inventory differs')
    f = np.array([list(map(float, row.split())) for row in rows[1:]])
    if not np.array_equal(f[:,0], np.arange(1,n+1)) or np.max(abs(f[:,1]-mu[:,0])) > 1e-12:
        raise InvalidArtifact('supplied field mapping/cavity differs')
    return mu[:,0], mu[:,7:10], mu[:,10:13], f[:,8:11], f[:,11:14]


def prepare(collection, software, plan, output):
    co = read_json(collection); pm = read_json(verify(co['manifest'])); sw = read_json(software)
    if not co['complete'] or not co['numerical_pass'] or sw['returncode']:
        raise InvalidArtifact('qualified native archive/software required')
    if read_json(verify(pm['software']))['library'] != sw['library']:
        raise InvalidArtifact('native library changed')
    root = Path(output).resolve(); root.mkdir(parents=True, exist_ok=False); impl = root/'implementation'; impl.mkdir()
    for name, pin in pm['implementation'].items(): shutil.copyfile(verify(pin), impl/name)
    shutil.copyfile(__file__, impl/Path(__file__).name)
    m = dict(protocol=PROTOCOL, collection=record(collection), software=record(software), plan=record(plan),
             implementation={p.name:record(p) for p in impl.glob('*.py')}, tolerances=TOL, threads=64,
             requested_replays=4, new_response_solves=0, new_DFT_calls=0, new_MACE_calls=0, tasks=[])
    for case in CASES:
        for metal in ('Ca','La'):
            st = next(t for t in pm['static_tasks'] if t['case_id']==case and t['state']==metal and t['variant']=='primary')
            rt = next(t for t in pm['response_tasks'] if t['case_id']==case and t['state']==metal and t['variant']=='primary' and not t['task_id'].endswith('_standard'))
            sp = co['tasks'][st['task_id']]['result']; rp = co['tasks'][rt['task_id']]['result']
            s = read_json(verify(sp)); r = read_json(verify(rp)); arrays = source_arrays(s,r)
            d = root/'tasks'/(case+'_'+metal); d.mkdir(parents=True)
            t = dict(task_id=d.name, case_id=case, metal=metal, static=sp, response=rp,
                     frozen_indices=st['frozen_indices'])
            for name in ('xyz','key','mask','overrides'):
                target = d/Path(st[name]['path']).name; shutil.copyfile(verify(st[name]),target); t[name]=record(target)
            rows = np.column_stack((np.arange(1,len(arrays[0])+1), *arrays))
            target = d/'frozen_response.dat'
            target.write_text(str(len(rows))+'\n'+''.join(str(int(row[0]))+' '+' '.join(format(v,'.17g') for v in row[1:])+'\n' for row in rows))
            t['input'] = record(target); t['cache_key']=task_key(t,m); m['tasks'].append(t)
    write_new(root/'manifest.json',m); return validate(root/'manifest.json')


def validate(path):
    m=read_json(path); sw=read_json(verify(m['software'])); p=read_json(verify(sw['parent']))
    if m['protocol']!=PROTOCOL or m['tolerances']!=TOL or m['threads']!=64 or len(m['tasks'])!=4:
        raise InvalidArtifact('native replay scope changed')
    for pin in (m['collection'],m['plan'],*m['implementation'].values(),sw['library'],sw['executable'],sw['frontend'],sw['implementation']): verify(pin)
    if frontend(verify(p['frontend']).read_text())!=verify(sw['frontend']).read_text():
        raise InvalidArtifact('native wrapper changed')
    if {(t['case_id'],t['metal']) for t in m['tasks']}!={(c,s) for c in CASES for s in ('Ca','La')}:
        raise InvalidArtifact('native replay input inventory changed')
    for t in m['tasks']:
        if t['cache_key']!=task_key(t,m): raise InvalidArtifact('native replay cache changed')
        for name in ('static','response','xyz','key','mask','overrides','input'): verify(t[name])
        s=read_json(verify(t['static'])); r=read_json(verify(t['response']))
        for a in (s,r): verify(a['receipt']['log'])
        arrays=source_arrays(s,r); expected=np.column_stack((np.arange(1,len(arrays[0])+1),*arrays))
        lines=verify(t['input']).read_text().splitlines()
        if int(lines[0])!=len(expected) or not np.array_equal(np.loadtxt(lines[1:]),expected):
            raise InvalidArtifact('frozen response differs from real archive')
    return m


def parse(text,t):
    if 'ALQUEMIA_FROZEN_RESPONSE_COMPLETE' not in text:
        raise InvalidArtifact('native replay incomplete')
    s=read_json(verify(t['static'])); r=read_json(verify(t['response'])); born,md,mp,fd,fp=source_arrays(s,r)
    marker='ALQUEMIA_MOMENTS_COMPLETE'
    moments=parse_moments(text[:text.index(marker)]+marker,None)
    if not moments['pass_'] or moments['parameters']!=s['parameters']:
        raise InvalidArtifact('native physical system changed')
    if any(a['global_']!=b['global_'] for a,b in zip(moments['atoms'],s['moments'])):
        raise InvalidArtifact('native global multipoles changed')
    mutual={}; reaction={}; energies={}; nonpolar=None; timing=None
    for line in text.splitlines():
        f=line.split()
        if not f: continue
        if f[0] in ('MUTUAL','REACTION','STATIC_GK'):
            target={'MUTUAL':mutual,'REACTION':reaction,'STATIC_GK':energies}[f[0]]; i=int(f[1])
            if i in target: raise InvalidArtifact('duplicate native row')
            target[i]=list(map(float,f[2:]))
        elif f[0]=='NONPOLAR_ENERGY': nonpolar=list(map(float,f[1:]))
        elif f[0]=='KERNEL_TIMING': timing=list(map(float,f[1:]))
    n=len(born)
    if set(mutual)!=set(range(1,n+1)) or set(reaction)!=set(mutual) or set(energies)!={1,2,3} or nonpolar is None or timing is None:
        raise InvalidArtifact('missing native operator/component evidence')
    u=np.array([mutual[i] for i in range(1,n+1)]); rf=np.array([reaction[i] for i in range(1,n+1)])
    if u.shape!=(n,14) or rf.shape!=(n,6): raise InvalidArtifact('native vector shape differs')
    c=s['electric']/s['dielec']; active=u[:,1]>0; active[np.array(t['frozen_indices'])-1]=False
    if np.any(md[~active]) or np.any(mp[~active]): raise InvalidArtifact('nonresponsive dipole changed')
    ad=np.zeros_like(md); ap=np.zeros_like(mp)
    ad[active]=md[active]/u[active,1,None]-u[active,8:11]
    ap[active]=mp[active]/u[active,1,None]-u[active,11:14]
    dot=lambda a,b: math.fsum((a[active]*b[active]).ravel())
    induction=-.5*c*dot(md,fp)
    stationary=.5*c*(dot(mp,ad)-dot(mp,fd)-dot(md,fp))
    cross=energies[2][0]-energies[3][0]
    field_cross=energies[1][0]-.5*c*(dot(md+mp,rf[:,:3])+dot(mp,u[:,8:11]-u[:,2:5]))
    # Native Tinker debye constant; separately printed below for runtime verification.
    debye_line=next(line.split() for line in text.splitlines() if line.split()[:1]==['NATIVE_DEBYE'])
    debye=float(debye_line[1])
    residuals=[debye*math.sqrt(math.fsum((((a[active]-f[active])*u[active,1,None])**2).ravel())/n)
               for a,f in ((ad,fd),(ap,fp))]
    errors=dict(born_A=float(np.max(abs(u[:,0]-born))),
                static_kcal=energies[1][0]-(s['energies_kcal_mol']['solvation_with_nonpolar']-sum(nonpolar)),
                reaction_field=float(np.max(abs(rf[:,:3]-rf[:,3:]))),
                cross_kcal=cross-field_cross, stationary_kcal=stationary-induction,
                residual_Debye=max(residuals), reciprocity_kcal=.5*c*(dot(mp,ad)-dot(md,ap)))
    return dict(status='computed',checks=[dict(name=k,error=v,tolerance=TOL[k],pass_=abs(v)<=TOL[k]) for k,v in errors.items()],
                stationary_kcal=stationary, induction_kcal=induction, GK_cross_kcal=cross,
                GK_cross_field_kcal=field_cross, static_GK_kcal={str(k):v[0] for k,v in energies.items()},
                nonpolar_kcal=nonpolar, residual_RMS_Debye=residuals,
                kernel_wall_seconds=timing[0],kernel_CPU_seconds=timing[1],new_score=None)


def execute(path,output):
    m=validate(path); sw=read_json(verify(m['software'])); root=Path(path).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=64:
        raise InvalidArtifact('64-CPU allocation required')
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for t in m['tasks']:
            d=root/'tasks'/t['task_id']; out=d/'attempt_0001'; out.mkdir()
            cmd=[str(verify(sw['executable'])),str(verify(t['xyz'])),str(verify(t['mask'])),str(verify(t['overrides'])),'replay',str(verify(t['input']))]
            start=time.monotonic(); before=resource.getrusage(resource.RUSAGE_CHILDREN)
            env={**os.environ,'OMP_NUM_THREADS':'64','OMP_MAX_ACTIVE_LEVELS':'1','OMP_PROC_BIND':'spread','OMP_PLACES':'cores','OPENBLAS_NUM_THREADS':'1'}
            with (out/'native.log').open('x') as log:
                proc=subprocess.run(cmd,cwd=d,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env)
            after=resource.getrusage(resource.RUSAGE_CHILDREN)
            receipt=dict(command=cmd,returncode=proc.returncode,wall_seconds=time.monotonic()-start,
                         CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                         peak_RSS_KiB=after.ru_maxrss,log=record(out/'native.log'),slurm_job_id=os.environ['SLURM_JOB_ID'])
            r=dict(task_id=t['task_id'],cache_key=t['cache_key'],manifest=record(path),receipt=receipt,status='failed')
            try:
                if proc.returncode: raise InvalidArtifact('native replay failed')
                r.update(parse(verify(receipt['log']).read_text(),t))
            except Exception as exc: r['failure_reason']=str(exc)
            write_new(out/'result.json',r); print(t['task_id'],r['status'],r.get('failure_reason',''),flush=True)
    return collect(path,output)


def collect(path,output):
    m=validate(path); root=Path(path).resolve().parent; rows={}
    for t in m['tasks']:
        p=root/'tasks'/t['task_id']/'attempt_0001/result.json'; r=read_json(p)
        if r['cache_key']!=t['cache_key'] or r['manifest']!=record(path): raise InvalidArtifact('native replay receipt mismatch')
        verify(r['receipt']['log']); rows[t['task_id']]=dict(result=record(p),**r)
    complete=all(r['status']=='computed' for r in rows.values())
    result=dict(protocol=PROTOCOL,manifest=record(path),rows=rows,complete=complete,
                checks_pass=complete and all(c['pass_'] for r in rows.values() for c in r['checks']),
                new_response_solves=0,new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False)
    write_new(output,result); return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('build'); p.add_argument('--parent',required=True); p.add_argument('--output',required=True)
    p=sub.add_parser('prepare')
    for k in ('collection','software','plan','output'): p.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        p=sub.add_parser(op); p.add_argument('--manifest',required=True)
        if op!='dry-run': p.add_argument('--output',required=True)
    a=parser.parse_args()
    if a.command=='build': r=build(a.parent,a.output)
    elif a.command=='prepare': r=prepare(a.collection,a.software,a.plan,a.output)
    elif a.command=='dry-run': r=validate(a.manifest)
    elif a.command=='execute': r=execute(a.manifest,a.output)
    else: r=collect(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','tasks','implementation')}))
